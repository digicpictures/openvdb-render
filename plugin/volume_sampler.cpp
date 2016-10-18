#include "volume_sampler.h"

#include <algorithm>
#include <atomic>
#include <random>
#include <thread>

#include <boost/algorithm/clamp.hpp>

#include <maya/MHWGeometry.h>
#include <maya/MShaderManager.h>
#include <maya/MTextureManager.h>

#include <openvdb/tools/Dense.h>
#include <openvdb/tools/GridTransformer.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/ValueTransformer.h>
#include <openvdb/tree/LeafManager.h>
#include <openvdb/Exceptions.h>

#include <tbb/parallel_for.h>

#include "util.h"

// Disable "decorated name length exceeded, name was truncated" warning.
#pragma warning(disable: 4503)

using namespace MHWRender;
using namespace openvdb;

namespace {

    template <typename T>
    using identity_t = T;

    template <typename T>
    T unlerp(identity_t<T> a, identity_t<T> b, T x)
    {
        return (x - a) / (b - a);
    }

    template <typename VecT>
    typename VecT::value_type maxComponentValue(const VecT& v)
    {
        return std::max(std::max(v.x(), v.y()), v.z());
    }

    MHWRender::MTextureManager* getTextureManager()
    {
        return MHWRender::MRenderer::theRenderer()->getTextureManager();
    }

} // unnamed namespace

VolumeTexture VolumeSampler::sampleGridWithBoxFilter(const openvdb::FloatGrid& grid, const openvdb::Coord& texture_extents)
{
    const auto sampler = openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler>(grid);
    const auto world_bbox = grid.transform().indexToWorld(getIndexSpaceBoundingBox(&grid));
    const auto domain_extents = texture_extents.asVec3d();
    auto sampling_func = [&sampler, &world_bbox, &domain_extents](const openvdb::Vec3d& domain_index) {
        const auto sample_pos_ws = (domain_index + 0.5) / domain_extents * world_bbox.extents() + world_bbox.min();
        return sampler.wsSample(sample_pos_ws);
    };
    return sampleVolume(texture_extents, sampling_func);
}

VolumeTexture VolumeSampler::sampleGridWithMipmapFilter(const openvdb::FloatGrid& grid, const openvdb::Coord& texture_extents)
{
    // Calculate LOD level.
    const auto index_bbox = getIndexSpaceBoundingBox(&grid);
    const auto grid_extents = index_bbox.extents().asVec3d();
    const auto coarse_voxel_size = grid_extents / texture_extents.asVec3d();
    const auto max_levels = openvdb::math::Ceil(std::log2(maxComponentValue(grid_extents)));
    const auto lod_level = boost::algorithm::clamp(std::log2(maxComponentValue(coarse_voxel_size)), 0, max_levels);
    const auto num_levels = size_t(openvdb::math::Ceil(lod_level)) + 1;

    if (num_levels == 1) {
        // No need for mult res grid.
        return sampleGridWithBoxFilter(grid, texture_extents);
    }

    // Create and sample multi res grid.
    openvdb::tools::MultiResGrid<openvdb::FloatTree> multires(num_levels, grid);

    const auto world_bbox = grid.transform().indexToWorld(index_bbox);
    const auto domain_extents = texture_extents.asVec3d();
    auto sampling_func = [&multires, lod_level, &world_bbox, &domain_extents](const openvdb::Vec3d& domain_index) {
        const auto sample_pos_ws = (domain_index + 0.5) / domain_extents * world_bbox.extents() + world_bbox.min();
        const auto sample_pos_is = multires.transform().worldToIndex(sample_pos_ws);
        return multires.sampleValue<1>(sample_pos_is, lod_level);
    };
    return sampleVolume(texture_extents, sampling_func);
}

template <typename SamplingFunc>
VolumeTexture VolumeSampler::sampleVolume(const openvdb::Coord& extents, SamplingFunc sampling_func)
{
    const auto domain = openvdb::CoordBBox(openvdb::Coord(), extents - openvdb::Coord(1, 1, 1));
    m_buffer.resize(domain.volume());

    // Sample on a lattice.
    using PerThreadRange = tbb::enumerable_thread_specific<FloatRange>;
    PerThreadRange per_thread_ranges;
    const auto stride = openvdb::Vec3i(1, extents.x(), extents.x() * extents.y());
    tbb::parallel_for(domain, [&sampling_func, &stride, &per_thread_ranges, output = m_buffer.data()](const CoordBBox& bbox) {
        PerThreadRange::reference this_thread_range = per_thread_ranges.local();
        for (auto z = bbox.min().z(); z <= bbox.max().z(); ++z) {
            for (auto y = bbox.min().y(); y <= bbox.max().y(); ++y) {
                for (auto x = bbox.min().x(); x <= bbox.max().x(); ++x) {
                    const auto domain_index = openvdb::Vec3i(x, y, z);
                    const auto linear_index = domain_index.dot(stride);
                    const auto sample_value = sampling_func(domain_index);
                    output[linear_index] = sample_value;
                    this_thread_range.update(sample_value);
                }
            }
        }
    });
    FloatRange value_range;
    for (const FloatRange& per_thread_range : per_thread_ranges) {
        value_range.update(per_thread_range);
    }

    // Remap sample values to [0, 1].
    using tbb_range = tbb::blocked_range<size_t>;
    tbb::parallel_for(tbb_range(0, m_buffer.size()),
                      [buffer = m_buffer.data(), &value_range](const tbb_range& range) {
        for (auto i = range.begin(); i < range.end(); ++i) {
            buffer[i] = unlerp(value_range.min, value_range.max, buffer[i]);
        }
    });

    return { extents, m_buffer.data(), value_range };
}

MHWRender::MTexture*
VolumeTexture::acquireVolumeTexture(const openvdb::Coord& texture_extents, const float* pixel_data)
{
    MHWRender::MTextureDescription texture_desc;
    texture_desc.fWidth = texture_extents.x();
    texture_desc.fHeight = texture_extents.y();
    texture_desc.fDepth = texture_extents.z();
    texture_desc.fBytesPerRow = 4 * texture_desc.fWidth;
    texture_desc.fBytesPerSlice = texture_desc.fBytesPerRow * texture_desc.fHeight;
    texture_desc.fMipmaps = 0;
    texture_desc.fArraySlices = 1;
    texture_desc.fFormat = kR32_FLOAT;
    texture_desc.fTextureType = kVolumeTexture;
    texture_desc.fEnvMapType = kEnvNone;
    return getTextureManager()->acquireTexture("", texture_desc, pixel_data, true);
}
