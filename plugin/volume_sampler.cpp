#include "volume_sampler.h"

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

Texture
VolumeSampler::sampleMultiResGrid(const openvdb::tools::MultiResGrid<openvdb::FloatTree>& multires, const openvdb::Coord& texture_extents)
{
    const CoordBBox fine_bbox = multires.grid(0)->evalActiveVoxelBoundingBox();

    // Calculate LOD level.
    const auto coarse_voxel_size = fine_bbox.extents().asVec3d() / texture_extents.asVec3d();
    const auto max_component = std::max(std::max(coarse_voxel_size.x(), coarse_voxel_size.y()), coarse_voxel_size.z());
    const auto lod_level = boost::algorithm::clamp(std::log2(max_component), 0, double(multires.coarsestLevel()));

    float min_value = 0.0, max_value = 0.0;
    multires.grid(0)->evalMinMax(min_value, max_value);
    auto samples = samplingLoop(fine_bbox, texture_extents, [&multires, lod_level, min_value, max_value](const openvdb::Vec3d& sample_pos) {
            return unlerp(min_value, max_value, multires.sampleValue<1>(sample_pos, lod_level));
        });

    Texture res;
    res.texture = acquireVolumeTexture(texture_extents, samples.data());
    res.min_value = min_value;
    res.max_value = max_value;
    return res;
}

Texture
VolumeSampler::sampleGrid(const openvdb::FloatGrid& grid, const openvdb::Coord& texture_extents)
{
    const CoordBBox fine_bbox = grid.evalActiveVoxelBoundingBox();
    //auto sampler = openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler>(grid);
    auto sampler = openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::PointSampler>(grid);

    float min_value = 0.0, max_value = 0.0;
    grid.evalMinMax(min_value, max_value);
    auto samples = samplingLoop(fine_bbox, texture_extents,
        [&sampler, min_value, max_value, &sample_jitter=m_sample_jitter](const openvdb::Vec3d& sample_pos) {
        //const int sx = int(sample_pos.x());
        //const int sy = int(sample_pos.y());
        //const int sz = int(sample_pos.z());
        //[&sampler, min_value, max_value](const openvdb::Vec3d& sample_pos) {
        float res = 0;
        const float weight = 1.0f / (N_FILT_SAMP * N_FILT_SAMP * N_FILT_SAMP);
        auto jitter_it = sample_jitter.begin();
        for (int k = 0; k < N_FILT_SAMP; ++k) {
            //const int z = (k + sz) % N_FILT_SAMP;
            //const int base_z = z * N_FILT_SAMP * N_FILT_SAMP;
            for (int j = 0; j < N_FILT_SAMP; ++j) {
                //const int y = (j + sy) % N_FILT_SAMP;
                //const int base_y = y * N_FILT_SAMP + base_z;
                for (int i = 0; i < N_FILT_SAMP; ++i) {
                    //const int x = (i + sx) % N_FILT_SAMP;
                    //auto jitter = sample_jitter[x + base_y];
                    auto jitter = *jitter_it++;
                    res += weight * unlerp(min_value, max_value, sampler.isSample(sample_pos + jitter));
                }
            }
        }
        return res;
        //return unlerp(min_value, max_value, sampler.isSample(sample_pos));
    });

    Texture res;
    res.texture = acquireVolumeTexture(texture_extents, samples.data());
    res.min_value = min_value;
    res.max_value = max_value;
    return res;
}

std::vector<float>
VolumeSampler::samplingLoop(const openvdb::CoordBBox& input_domain, const openvdb::Coord& output_extents,
                            std::function<float(openvdb::Vec3d)> sampling_func)
{
    // BBox splitting treats both ends of the bbox as inclusive, so subtract 1 from the slice counts to get inclusive bbox.
    const auto coarse_bbox = openvdb::CoordBBox(openvdb::Coord(), output_extents - openvdb::Coord(1, 1, 1));

    // Create coarse -> fine index transform funcition.
    const auto coarse_voxel_size = input_domain.extents().asVec3d() / output_extents.asVec3d();
    const auto index_transform = [&input_domain, &coarse_voxel_size](const Coord& coarse_index) {
        return input_domain.min().asVec3d() + coarse_index.asVec3d() * coarse_voxel_size;
    };

    // Iterate through coarse voxels and pass in fine coords to sampling func.
    std::vector<float> output(output_extents.x() * output_extents.y() * output_extents.z());
    const auto stride = openvdb::Vec3i(1, output_extents.x(), output_extents.x() * output_extents.y());
    tbb::parallel_for(coarse_bbox, [&index_transform, &stride, &sampling_func, &output](const CoordBBox& bbox) {
        typedef int32_t index_t;
        for (index_t z = bbox.min().z(); z <= bbox.max().z(); ++z) {
            const index_t base_z = z * stride.z();
            for (index_t y = bbox.min().y(); y <= bbox.max().y(); ++y) {
                const index_t base_y = base_z + y * stride.y();
                for (index_t x = bbox.min().x(); x <= bbox.max().x(); ++x) {
                    const auto out_index = base_y + x;
                    const auto sample_pos = index_transform(openvdb::Coord(x, y, z));
                    output[out_index] = sampling_func(sample_pos);
                }
            }
        }
    });

    return output;
}

MHWRender::MTexture* VolumeSampler::acquireVolumeTexture(const openvdb::Coord& texture_extents, const float* pixel_data)
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
    return m_texture_manager->acquireTexture("", texture_desc, pixel_data, true);
}

void VolumeSampler::initSampleJitter()
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
                *it++ = N_FILT_SAMP * (openvdb::Vec3d(x, y, z) - openvdb::Vec3d(0.5, 0.5, 0.5));
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