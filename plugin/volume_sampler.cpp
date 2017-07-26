#include "volume_sampler.h"

#include <algorithm>

#include <openvdb/tools/Dense.h>
#include <openvdb/tools/GridTransformer.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/ValueTransformer.h>
#include <openvdb/tree/LeafManager.h>
#include <openvdb/Types.h>

#include <tbb/parallel_for.h>

#include "progress_bar.h"
#include "vdb_maya_utils.hpp"

template <typename SampleType>
bool VolumeSampler<SampleType>::sampleGrid(const openvdb::FloatGrid& grid, const openvdb::Coord& texture_extents, const VolumeBufferHandleType& output)
{
    assert(output);

    // Calculate number of LOD levels.
    const auto grid_extents = getIndexSpaceBoundingBox(grid).extents().asVec3d();
    const auto num_levels = size_t(openvdb::math::Ceil(std::log2(maxComponentValue(grid_extents))));

    if (num_levels > 1)
    {
        // Create and sample MultiResGrid.
        openvdb::tools::MultiResGrid<openvdb::FloatTree> multires(num_levels, grid);
        return sampleMultiResGrid(multires, texture_extents, output);
    }
    else
    {
        // Use box filter.
        return sampleGridWithBoxFilter(grid, texture_extents, output);
    }
}

namespace {

    template <typename T>
    struct ValueRange
    {
        T min, max;
        ValueRange() noexcept : min(std::numeric_limits<T>::max()), max(std::numeric_limits<T>::min()) {}
        ValueRange(T min_, T max_) noexcept : min(min_), max(max_) {}
        void update(T value) { min = std::min(min, value); max = std::max(max, value); }
        void update(const ValueRange<T>& value_range) { min = std::min(min, value_range.min); max = std::max(max, value_range.max); }
        ValueRange merge(const ValueRange<T>& value_range) { return{ std::min(min, value_range.min), std::max(max, value_range.max) }; }
    };
    typedef ValueRange<float> FloatRange;

    template <typename RealType>
    void setHeader(const FloatRange& value_range, const openvdb::BBoxd& bbox, VolumeBufferHeader<RealType>& output)
    {
        output.value_range[0] = static_cast<RealType>(value_range.min);
        output.value_range[1] = static_cast<RealType>(value_range.max);
        const auto& extents = bbox.extents();
        output.size[0] = static_cast<RealType>(extents.x());
        output.size[1] = static_cast<RealType>(extents.y());
        output.size[2] = static_cast<RealType>(extents.z());
        const auto& origin = bbox.min();
        output.origin[0] = static_cast<RealType>(origin.x());
        output.origin[1] = static_cast<RealType>(origin.y());
        output.origin[2] = static_cast<RealType>(origin.z());
    }

    template <typename SamplingFunc, typename VoxelType>
    bool sampleVolume(const openvdb::Coord& extents, SamplingFunc sampling_func, ProgressBar *progress_bar,
                      VoxelType *out_voxel_array, FloatRange& out_value_range)
    {
        const auto domain = openvdb::CoordBBox(openvdb::Coord(), extents - openvdb::Coord(1, 1, 1));
        const auto num_voxels = domain.volume();

        // Initialize progress bar.
        if (progress_bar)
            progress_bar->setMaxProgress(num_voxels);

        // Sample on a lattice.
        typedef tbb::enumerable_thread_specific<FloatRange> PerThreadRange;
        PerThreadRange per_thread_ranges;
        const auto stride = openvdb::Vec3i(1, extents.x(), extents.x() * extents.y());
        tbb::atomic<bool> cancelled;
        cancelled = false;
        tbb::parallel_for(domain, [&sampling_func, &stride, progress_bar, &cancelled,
            &per_thread_ranges, out_voxel_array](const openvdb::CoordBBox& bbox) {
            const auto local_extents = bbox.extents();
            const auto progress_step = local_extents.x() * local_extents.y();

            // Loop through local bbox.
            PerThreadRange::reference this_thread_range = per_thread_ranges.local();
            for (auto z = bbox.min().z(); z <= bbox.max().z(); ++z) {
                for (auto y = bbox.min().y(); y <= bbox.max().y(); ++y) {
                    for (auto x = bbox.min().x(); x <= bbox.max().x(); ++x) {
                        if (cancelled)
                            return;
                        const auto domain_index = openvdb::Vec3i(x, y, z);
                        const auto linear_index = domain_index.dot(stride);
                        const auto sample_value = sampling_func(domain_index);
                        out_voxel_array[linear_index] = sample_value;
                        this_thread_range.update(sample_value);
                    }
                }

                // Report progress.
                if (progress_bar)
                {
                    if (progress_bar->isCancelled())
                        cancelled = true;
                    progress_bar->addProgress(progress_step);
                }
            }
        });
        if (cancelled)
            return false;

        // Merge per-thread value ranges.
        out_value_range = FloatRange();
        for (const FloatRange& per_thread_range : per_thread_ranges) {
            out_value_range.update(per_thread_range);
        }

        // Remap sample values to [0, 1].
        typedef tbb::blocked_range<size_t> tbb_range;
        auto buffer = out_voxel_array;
        tbb::parallel_for(tbb_range(0, num_voxels),
            [buffer, &out_value_range](const tbb_range& range) {
            for (auto i = range.begin(); i < range.end(); ++i) {
                buffer[i] = unlerp(out_value_range.min, out_value_range.max, buffer[i]);
            }
        });

        if (progress_bar)
            progress_bar->setProgress(100);

        return true;
    }

} // unnamed namespace

template <typename SampleType>
bool VolumeSampler<SampleType>::sampleGridWithBoxFilter(const openvdb::FloatGrid& grid, const openvdb::Coord& texture_extents, const VolumeBufferHandleType& output)
{
    assert(output);

    // Return false if the grid bbox is empty.
    const auto grid_bbox_is = getIndexSpaceBoundingBox(grid);
    if (grid_bbox_is.empty())
        return false;

    // Set up sampling func.
    const auto sampler = openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler>(grid);
    const auto bbox_world = grid.transform().indexToWorld(grid_bbox_is);
    const auto domain_extents = texture_extents.asVec3d();
    auto sampling_func = [&sampler, &bbox_world, &domain_extents](const openvdb::Vec3d& domain_index) -> SampleType {
        const auto sample_pos_ws = (domain_index + 0.5) / domain_extents * bbox_world.extents() + bbox_world.min();
        return sampler.wsSample(sample_pos_ws);
    };

    FloatRange value_range;
    const auto res = sampleVolume(texture_extents, sampling_func, m_progress_bar, output.voxel_array, value_range);
    setHeader<SampleType>(value_range, bbox_world, *output.header);
    return res;
}

template <typename SampleType>
bool VolumeSampler<SampleType>::sampleMultiResGrid(const openvdb::tools::MultiResGrid<openvdb::FloatTree>& multires, const openvdb::Coord& texture_extents, const VolumeBufferHandleType& output)
{
    assert(output);

    // Return false if the grid bbox is empty.
    const auto grid_bbox_is = getIndexSpaceBoundingBox(*multires.grid(0));
    if (grid_bbox_is.empty())
        return false;

    // Calculate LOD level.
    const auto grid_extents = grid_bbox_is.extents().asVec3d();
    const auto coarse_voxel_size = grid_extents / texture_extents.asVec3d();
    const auto max_levels = openvdb::math::Ceil(std::log2(maxComponentValue(grid_extents)));
    const auto lod_level = clamp(std::log2(maxComponentValue(coarse_voxel_size)), 0, max_levels);

    // Set up sampling func.
    const auto bbox_world = multires.grid(0)->transform().indexToWorld(grid_bbox_is);
    const auto domain_extents = texture_extents.asVec3d();
    auto sampling_func = [&multires, lod_level, &bbox_world, &domain_extents](const openvdb::Vec3d& domain_index) -> SampleType {
        const auto sample_pos_ws = (domain_index + 0.5) / domain_extents * bbox_world.extents() + bbox_world.min();
        const auto sample_pos_is = multires.transform().worldToIndex(sample_pos_ws);
        return multires.sampleValue<1>(sample_pos_is, lod_level);
    };

    FloatRange value_range;
    const auto res = sampleVolume(texture_extents, sampling_func, m_progress_bar, output.voxel_array, value_range);
    setHeader<SampleType>(value_range, bbox_world, *output.header);
    return res;
}

template class VolumeSampler<float>;
template class VolumeSampler<half>;
