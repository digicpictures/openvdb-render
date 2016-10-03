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

using namespace MHWRender;

//==============================================================================
// CLASS ProgressBar
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

MTexture* volumeTextureFromGrid(const openvdb::FloatGrid* grid, const openvdb::Coord& slice_counts, MTextureManager* texture_manager)
{
    ProgressBar pb("Resampling VDB grid...", 100);
    pb.stepProgress();

#if 1
    using namespace openvdb;

    //Timer timer;

    CoordBBox grid_bbox;
    grid->tree().evalActiveVoxelBoundingBox(grid_bbox);

    // Scale and filter grid to slice_counts size.
    //timer.start("Scaling and filtering grid...");
    FloatGrid::Ptr filtered_grid = FloatGrid::create(grid->background());
    const auto scale_factor = slice_counts / grid_bbox.extents();
    tools::GridTransformer(math::scale<Mat4R>(scale_factor)).transformGrid<tools::BoxSampler, FloatGrid>(*grid, *filtered_grid);
    CoordBBox filtered_grid_bbox;
    filtered_grid->tree().evalActiveVoxelBoundingBox(filtered_grid_bbox);
    //timer.end();

    // Transform grid values to [0, 1].
    //timer.start("Rescaling grid values to [0, 1]...");
    float min_voxel_value, max_voxel_value;
    tree::LeafManager<const FloatTree> leafs(filtered_grid->tree());
    {
        MinMaxVoxel<const FloatTree> minmax(leafs);
        minmax.runParallel();
        min_voxel_value = minmax.minVoxel();
        max_voxel_value = minmax.maxVoxel();
    }
    tools::foreach(filtered_grid->beginValueAll(), [min_voxel_value, max_voxel_value](const FloatGrid::ValueAllIter& iter) {
        iter.setValue((*iter - min_voxel_value) / (max_voxel_value - min_voxel_value));
    });
    //timer.end();

    // Create dense grid.
    //timer.start("Creating dense grid...");
    using FloatDenseGrid = tools::Dense<float, tools::LayoutXYZ>;
    FloatDenseGrid dense_grid(filtered_grid_bbox);
    tools::copyToDense(*filtered_grid, dense_grid);
    //timer.end();

    // Create output texture.
    //out_volume.genTexture(filtered_grid_bbox.extents(), dense_grid.data());
    const auto extents = filtered_grid_bbox.extents();
#else
    const auto extents = slice_counts;
#endif

    MHWRender::MTextureDescription texture_desc;
    texture_desc.fWidth = extents.x();
    texture_desc.fHeight = extents.y();
    texture_desc.fDepth = extents.z();
    texture_desc.fBytesPerRow = 4 * texture_desc.fWidth;
    texture_desc.fBytesPerSlice = texture_desc.fBytesPerRow * texture_desc.fHeight;
    //texture_desc.fBytesPerRow = 0;
    //texture_desc.fBytesPerSlice = 0;
    texture_desc.fMipmaps = 1;
    texture_desc.fFormat = kR32_FLOAT;
    texture_desc.fTextureType = kVolumeTexture;
    texture_desc.fEnvMapType = kEnvNone;

    return texture_manager->acquireTexture("", texture_desc, dense_grid.data(), false);

#if 0
    std::vector<float> v(extents.x() * extents.y() * extents.z());
    for (int k = 0; k < extents.z(); ++k) {
        for (int j = 0; j < extents.y(); ++j) {
            for (int i = 0; i < extents.x(); ++i) {
                int idx = i + extents.x() * j + extents.x() * extents.y() * k;
                v[idx] = (float(i) / extents.x()) * (float(j) / extents.y()) * (float(k) / extents.z());
            }
        }
    }
    return texture_manager->acquireTexture("", texture_desc, v.data(), false);
#endif
}
