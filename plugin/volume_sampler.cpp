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

} // unnamed namespace

VolumeTexture
VolumeSampler::sampleMultiResGrid(const openvdb::FloatGrid& grid, const openvdb::Coord& texture_extents)
{
    // Calculate LOD level.
    const auto index_bbox = getIndexSpaceBoundingBox(&grid);
    const auto grid_extents = index_bbox.extents().asVec3d();
    const auto coarse_voxel_size = grid_extents / texture_extents.asVec3d();
    const auto max_levels = openvdb::math::Ceil(std::log2(maxComponent(grid_extents)));
    const auto lod_level = boost::algorithm::clamp(std::log2(maxComponent(coarse_voxel_size)), 0, max_levels);
    const auto num_levels = size_t(openvdb::math::Ceil(lod_level)) + 1;

    if (num_levels == 1) {
        // No need for mult res grid.
        return sampleGrid(grid, texture_extents);
    }

    // Create multi res grid.
    openvdb::tools::MultiResGrid<openvdb::FloatTree> multires(num_levels, grid);

    // Sample the grid on a uniform lattice.
    const auto world_bbox = grid.transform().indexToWorld(index_bbox);
    const auto domain = openvdb::CoordBBox(openvdb::Coord(), texture_extents - openvdb::Coord(1, 1, 1));
    m_buffer.resize(domain.volume());
    std::atomic<float> min_value = 0.0, max_value = 0.0;
    samplingLoop(m_buffer.data(), domain, [&multires, lod_level, domain_extents = texture_extents.asVec3d(), &world_bbox, &min_value, &max_value](const openvdb::Coord& index, float& output) {
        const auto sample_pos_ws = (index.asVec3d() + 0.5) / domain_extents * world_bbox.extents() + world_bbox.min();
        const auto sample_pos_is = multires.transform().worldToIndex(sample_pos_ws);
        const auto sample_value = multires.sampleValue<1>(sample_pos_is, lod_level);
        min_value.store(std::min(sample_value, min_value.load(std::memory_order::memory_order_relaxed)), std::memory_order::memory_order_relaxed);
        max_value.store(std::max(sample_value, max_value.load(std::memory_order::memory_order_relaxed)), std::memory_order::memory_order_relaxed);
        output = sample_value;
    });

    samplingLoop(m_buffer.data(), domain, [min_value = min_value.load(), max_value = max_value.load()](const openvdb::Coord, float& output){
        output = unlerp(min_value, max_value, output);
    });

    return { texture_extents, m_buffer.data(), ValueRange(min_value, max_value), m_texture_manager };
}

VolumeTexture
VolumeSampler::sampleGrid(const openvdb::FloatGrid& grid, const openvdb::Coord& texture_extents)
{
    const auto world_bbox = grid.transform().indexToWorld(getIndexSpaceBoundingBox(&grid));

    const auto domain = openvdb::CoordBBox(openvdb::Coord(), texture_extents - openvdb::Coord(1, 1, 1));
    m_buffer.resize(domain.volume());

    auto sampler = openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler>(grid);

    std::atomic<float> min_value = 0.0, max_value = 0.0;
    samplingLoop(m_buffer.data(), domain,
        [domain_extents = texture_extents.asVec3d(), &world_bbox, &sampler, &min_value, &max_value](const openvdb::Coord& index, float& output) {
        const auto sample_pos_ws = (index.asVec3d() + 0.5) / domain_extents * world_bbox.extents() + world_bbox.min();
        float sample_value = sampler.wsSample(sample_pos_ws);

        min_value.store(std::min(sample_value, min_value.load(std::memory_order::memory_order_relaxed)), std::memory_order::memory_order_relaxed);
        max_value.store(std::max(sample_value, max_value.load(std::memory_order::memory_order_relaxed)), std::memory_order::memory_order_relaxed);
        output = sample_value;
    });

    samplingLoop(m_buffer.data(), domain, [min_value = min_value.load(), max_value = max_value.load()](const openvdb::Coord, float& output){
        output = unlerp(min_value, max_value, output);
    });

    return { texture_extents, m_buffer.data(), ValueRange(min_value, max_value), m_texture_manager };
}

void
VolumeSampler::samplingLoop(float* output, const openvdb::CoordBBox& domain, std::function<void(openvdb::Coord,float&)> sampling_func)
{

    // Iterate through coarse voxels and pass in fine coords to sampling func.
    const auto extents = domain.extents();
    const auto stride = openvdb::Vec3i(1, extents.x(), extents.x() * extents.y());
    tbb::parallel_for(domain, [&stride, &sampling_func, output](const CoordBBox& bbox) {
        typedef int32_t index_t;
        for (index_t z = bbox.min().z(); z <= bbox.max().z(); ++z) {
            const index_t base_z = z * stride.z();
            for (index_t y = bbox.min().y(); y <= bbox.max().y(); ++y) {
                const index_t base_y = base_z + y * stride.y();
                for (index_t x = bbox.min().x(); x <= bbox.max().x(); ++x) {
                    const auto out_index = base_y + x;
                    sampling_func(openvdb::Coord(x, y, z), output[out_index]);
                }
            }
        }
    });
}

MHWRender::MTexture*
VolumeTexture::acquireVolumeTexture(const openvdb::Coord& texture_extents, const float* pixel_data, MHWRender::MTextureManager* texture_manager)
{
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
    return texture_manager->acquireTexture("", texture_desc, pixel_data, true);
}
