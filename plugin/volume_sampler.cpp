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
VolumeSampler::sampleMultiResGrid(const openvdb::tools::MultiResGrid<openvdb::FloatTree>& multires, const openvdb::Coord& texture_extents)
{
    const auto& grid = multires.grid(0);
    const auto index_bbox = getIndexSpaceBoundingBox(&grid);
    const auto world_bbox = grid->transform().indexToWorld(index_bbox);

    const auto domain = openvdb::CoordBBox(openvdb::Coord(), texture_extents - openvdb::Coord(1, 1, 1));
    m_buffer.resize(domain.volume());

    // Calculate LOD level.
    const auto coarse_voxel_size = index_bbox.extents().asVec3d() / texture_extents.asVec3d();
    const auto max_component = std::max(std::max(coarse_voxel_size.x(), coarse_voxel_size.y()), coarse_voxel_size.z());
    const auto lod_level = boost::algorithm::clamp(std::log2(max_component), 0, double(multires.coarsestLevel()));

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

    auto sampler = openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::PointSampler>(grid);

    std::atomic<float> min_value = 0.0, max_value = 0.0;
    samplingLoop(m_buffer.data(), domain,
        [&texture_extents, &sample_jitter = m_sample_jitter, &world_bbox, &sampler, &min_value, &max_value](const openvdb::Coord& index, float& output) {

        auto jitter_stride = openvdb::Vec3i(1, N_FILT_SAMP, N_FILT_SAMP * N_FILT_SAMP);
        // Use prime number of different jitter arrangements to reduce the chance of banding -- hence the 5.
        const int permut = index.asVec3i().dot(texture_extents.asVec3i()) % 5;
        for (int i = 0; i < permut; ++i) {
            std::next_permutation(jitter_stride.asPointer(), jitter_stride.asPointer() + 3);
        }

        const auto domain_extents = texture_extents.asVec3d();
        const auto jitter_mag = (world_bbox.extents() / domain_extents) / N_FILT_SAMP;
        const auto sample_pos_ws = (index.asVec3d() + 0.5) / domain_extents * world_bbox.extents() + world_bbox.min();
        const float weight = 1.0f / (N_FILT_SAMP * N_FILT_SAMP * N_FILT_SAMP);
        const auto n_half = openvdb::Vec3d(0.5 * N_FILT_SAMP);
        float sample_value = 0;
        for (int k = 0; k < N_FILT_SAMP; ++k) {
            for (int j = 0; j < N_FILT_SAMP; ++j) {
                for (int i = 0; i < N_FILT_SAMP; ++i) {
                    const auto sample_idx = openvdb::Vec3i(i, j, k);
                    const auto jitter_idx = sample_idx.dot(jitter_stride);
                    const auto jitter = (sample_jitter[jitter_idx] + sample_idx - n_half) * jitter_mag;
                    sample_value += weight * sampler.wsSample(sample_pos_ws + jitter);
                }
            }
        }

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

void
VolumeSampler::initSampleJitter()
{
    std::default_random_engine gen;
    std::uniform_real_distribution<float> dist(0.0, 1.0);

    // Canonical arrangement.
    m_sample_jitter.resize(N_FILT_SAMP * N_FILT_SAMP * N_FILT_SAMP);
    auto it = m_sample_jitter.begin();
    for (int k = 0; k < N_FILT_SAMP; ++k) {
        for (int j = 0; j < N_FILT_SAMP; ++j) {
            for (int i = 0; i < N_FILT_SAMP; ++i) {
                const float x = (i + (j + (k + dist(gen)) / N_FILT_SAMP) / N_FILT_SAMP) / N_FILT_SAMP;
                const float y = (j + (k + (i + dist(gen)) / N_FILT_SAMP) / N_FILT_SAMP) / N_FILT_SAMP;
                const float z = (k + (i + (j + dist(gen)) / N_FILT_SAMP) / N_FILT_SAMP) / N_FILT_SAMP;
                *it++ = openvdb::Vec3d(x, y, z) - openvdb::Vec3d(0.5, 0.5, 0.5);
            }
        }
    }

    // Shuffle -- 1st dimension.
    for (int i = 0; i < N_FILT_SAMP; ++i) {
        const int i_new = static_cast<int>(i + dist(gen) * (N_FILT_SAMP - i));
        for (int k = 0; k < N_FILT_SAMP; ++k) {
            for (int j = 0; j < N_FILT_SAMP; ++j) {
                std::swap(m_sample_jitter[i     + j * N_FILT_SAMP + k * N_FILT_SAMP * N_FILT_SAMP],
                          m_sample_jitter[i_new + j * N_FILT_SAMP + k * N_FILT_SAMP * N_FILT_SAMP]);
            }
        }
    }

    // Shuffle -- 2nd dimension.
    for (int j = 0; j < N_FILT_SAMP; ++j) {
        const int j_new = static_cast<int>(j + dist(gen) * (N_FILT_SAMP - j));
        for (int k = 0; k < N_FILT_SAMP; ++k) {
            for (int i = 0; i < N_FILT_SAMP; ++i) {
                std::swap(m_sample_jitter[i + j     * N_FILT_SAMP + k * N_FILT_SAMP * N_FILT_SAMP],
                          m_sample_jitter[i + j_new * N_FILT_SAMP + k * N_FILT_SAMP * N_FILT_SAMP]);
            }
        }
    }

    // Shuffle -- 3rd dimension.
    for (int k = 0; k < N_FILT_SAMP; ++k) {
        const int k_new = static_cast<int>(k + dist(gen) * (N_FILT_SAMP - k));
        for (int j = 0; j < N_FILT_SAMP; ++j) {
            for (int i = 0; i < N_FILT_SAMP; ++i) {
                std::swap(m_sample_jitter[i + j * N_FILT_SAMP + k     * N_FILT_SAMP * N_FILT_SAMP],
                          m_sample_jitter[i + j * N_FILT_SAMP + k_new * N_FILT_SAMP * N_FILT_SAMP]);
            }
        }
    }
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
