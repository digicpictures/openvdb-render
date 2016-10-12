#include "volume_sampler.h"

#include <thread>

#include <maya/MHWGeometry.h>
#include <maya/MShaderManager.h>
#include <maya/MTextureManager.h>

#include <openvdb/tools/Dense.h>
#include <openvdb/tools/GridTransformer.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/ValueTransformer.h>
#include <openvdb/tree/LeafManager.h>
#include <openvdb/Exceptions.h>

// Disable "decorated name length exceeded, name was truncated" warning.
#pragma warning(disable: 4503)

using namespace MHWRender;

//==============================================================================
class ProgressBar
{
public:
    ProgressBar(const MString& msg, unsigned int max)
    {
        // Display a progress bar if Maya is running in UI mode
        fShowProgress = (MGlobal::mayaState() == MGlobal::kInteractive);
        reset(msg, max);
    }
    void reset(const MString& msg, unsigned int max)
    {
        MStatus status;
        beginProgress(msg, max);
    }
    ~ProgressBar()
    {
        endProgress();
    }
    void stepProgress() const
    {
        if (fShowProgress) {
            MGlobal::executeCommand("progressBar -e -s 1 $gMainProgressBar");
        }
    }
    bool isCancelled() const
    {
        int isCancelled = 0;
        if (fShowProgress) {
            MGlobal::executeCommand("progressBar -q -ic $gMainProgressBar", isCancelled);
        }
        if (isCancelled) {
            MStatus status;
            const MString interruptMsg = "Interrupted by user";
            MGlobal::displayInfo(interruptMsg);
            return true;
        }
        return false;
    }
private:
    // Forbidden and not implemented.
    ProgressBar(const ProgressBar&);
    const ProgressBar& operator=(const ProgressBar&);
    void beginProgress(const MString& msg, unsigned int max) const
    {
        if (fShowProgress) {
            MString maxValue, progressBarCmd;
            // Progress from 0 to max
            if (max <= 0) {
                max = 1;
            }
            maxValue += max;
            // Clear previous isCancelled flag
            MGlobal::executeCommand("progressBar -e -bp -ii 1 $gMainProgressBar");
            MGlobal::executeCommand("progressBar -e -ep $gMainProgressBar");
            // Initialize the progress bar
            progressBarCmd.format("progressBar -e -bp -ii 1 -st \"^1s\" -max ^2s $gMainProgressBar",
                msg, maxValue);
            MGlobal::executeCommand(progressBarCmd);
        }
    }
    void endProgress() const
    {
        if (fShowProgress) {
            MGlobal::executeCommand("progressBar -e -ep $gMainProgressBar");
        }
    }
    bool fShowProgress;  // whether to show the progress bar
};

//==============================================================================
template<class TreeType>
class MinMaxVoxel
{
public:
    typedef openvdb::tree::LeafManager<TreeType> LeafArray;
    typedef typename TreeType::ValueType ValueType;

    MinMaxVoxel(LeafArray& leafs)
        : mLeafArray(leafs), mMin(std::numeric_limits<ValueType>::max()), mMax(-mMin) {}
    MinMaxVoxel(const MinMaxVoxel<TreeType>& rhs, tbb::split)
        : mLeafArray(rhs.mLeafArray), mMin(std::numeric_limits<ValueType>::max()), mMax(-mMin) {}

    void runParallel() { tbb::parallel_reduce(mLeafArray.getRange(), *this); }
    void runSerial() { (*this)(mLeafArray.getRange()); }

    const ValueType& minVoxel() const { return mMin; }
    const ValueType& maxVoxel() const { return mMax; }

    inline void operator()(const tbb::blocked_range<size_t>& range)
    {
        typename TreeType::LeafNodeType::ValueOnCIter iter;

        for (size_t n = range.begin(); n < range.end(); ++n) {
            iter = mLeafArray.leaf(n).cbeginValueOn();
            for (; iter; ++iter) {
                const ValueType value = iter.getValue();
                mMin = std::min(mMin, value);
                mMax = std::max(mMax, value);
            }
        }
    }

    inline void join(const MinMaxVoxel<TreeType>& rhs)
    {
        mMin = std::min(mMin, rhs.mMin);
        mMax = std::max(mMax, rhs.mMax);
    }

private:
    LeafArray& mLeafArray;
    ValueType mMin, mMax;
};

openvdb::Vec3d operator/(const openvdb::Coord& lhs, const openvdb::Coord& rhs)
{
    return{ double(lhs.x()) / rhs.x(), double(lhs.y()) / rhs.y(), double(lhs.z()) / rhs.z() };
}

//==============================================================================
MTexture* volumeTextureFromGrid(const openvdb::tools::MultiResGrid<openvdb::FloatTree>& grid, const openvdb::Coord& texture_extents, MTextureManager* texture_manager)
{
    //ProgressBar pb("Resampling VDB grid...", 100);
    //pb.stepProgress();

    using namespace openvdb;

    // Calculate the bounding box of the finest grid.
    CoordBBox fine_bbox = grid.grid(0)->evalActiveVoxelBoundingBox();

    // BBox splitting treats both ends of the bbox as inclusive, so subtract 1 from the slice counts to get inclusive bbox.
    const auto coarse_bbox = openvdb::CoordBBox(openvdb::Coord(), texture_extents - openvdb::Coord(1, 1, 1));

    // Create coarse -> fine index transform funcition.
    const auto coarse_voxel_size = fine_bbox.extents() / texture_extents;
    const auto index_transform = [&fine_bbox, &coarse_voxel_size](const Coord& coarse_index) {
        return fine_bbox.min().asVec3d() + coarse_index.asVec3d() * coarse_voxel_size;
    };

    // Calculate LOD level.
    const auto max_component = std::max(std::max(coarse_voxel_size.x(), coarse_voxel_size.y()), coarse_voxel_size.z());
    const auto lod_level = std::log2(max_component);

    // Get min and max voxel values and create value transform function.
    // TODO: leave values as is and use shader parameter for the value range.
    float min_voxel_value, max_voxel_value;
    tree::LeafManager<const FloatTree> leafs(grid.grid(0)->tree());
    {
        MinMaxVoxel<const FloatTree> minmax(leafs);
        minmax.runParallel();
        min_voxel_value = minmax.minVoxel();
        max_voxel_value = minmax.maxVoxel();
    }
    const auto value_transform = [min_voxel_value, max_voxel_value](float value) {
        return (value - min_voxel_value) / (max_voxel_value - min_voxel_value);
    };

    // Sample the multires grid on the sparse lattice using the LOD level.

    std::vector<float> samples(texture_extents.x() * texture_extents.y() * texture_extents.z());
    const auto stride = openvdb::Vec3i(1, texture_extents.x(), texture_extents.x() * texture_extents.y());

    tbb::parallel_for(coarse_bbox, [&grid, &index_transform, lod_level, &value_transform, &stride, &samples](const CoordBBox& bbox) {
        constexpr int SAMPLING_ORDER = 1;
        typedef int32_t index_t;
        for (index_t z = bbox.min().z(); z <= bbox.max().z(); ++z) {
            const index_t base_z = z * stride.z();
            for (index_t y = bbox.min().y(); y <= bbox.max().y(); ++y) {
                const index_t base_y = base_z + y * stride.y();
                for (index_t x = bbox.min().x(); x <= bbox.max().x(); ++x) {
                    const auto sample_pos = index_transform(openvdb::Coord(x, y, z));
                    const auto sample_value = grid.sampleValue<SAMPLING_ORDER>(sample_pos, lod_level);
                    const auto out_index = base_y + x;
                    samples[out_index] = value_transform(sample_value);
                }
            }
        }
    });


    // Create texture from the samples.

    MHWRender::MTextureDescription texture_desc;
    texture_desc.fWidth = texture_extents.x();
    texture_desc.fHeight = texture_extents.y();
    texture_desc.fDepth = texture_extents.z();
    texture_desc.fBytesPerRow = 4 * texture_desc.fWidth;
    texture_desc.fBytesPerSlice = texture_desc.fBytesPerRow * texture_desc.fHeight;
    texture_desc.fMipmaps = 1;
    texture_desc.fArraySlices = 1;
    texture_desc.fFormat = kR32_FLOAT;
    texture_desc.fTextureType = kVolumeTexture;
    texture_desc.fEnvMapType = kEnvNone;

    return texture_manager->acquireTexture("", texture_desc, samples.data(), true);
}
