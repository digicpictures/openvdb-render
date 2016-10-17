#pragma once

#include <type_traits>
#include <vector>

#include <openvdb/openvdb.h>
#include <openvdb/tools/MultiResGrid.h>

namespace MHWRender {
    class MTexture;
    class MTextureManager;
}

struct ValueRange
{
    float min, max;
    ValueRange() : min(0), max(0) {}
    ValueRange(float min_ , float max_) : min(min_), max(max_) {}
};

struct VolumeTexture
{
    MHWRender::MTexture* texture;
    ValueRange value_range;

    VolumeTexture() : texture(nullptr) {}
    VolumeTexture(const openvdb::Coord& extents, const float* data_norm, const ValueRange& value_range_, MHWRender::MTextureManager* texture_manager)
        : texture(acquireVolumeTexture(extents, data_norm, texture_manager)), value_range(value_range_) {}

private:
    MHWRender::MTexture* acquireVolumeTexture(const openvdb::Coord& texture_extents, const float* pixel_data, MHWRender::MTextureManager* texture_manager);
};

class VolumeSampler
{
public:
    VolumeSampler(MHWRender::MTextureManager* texture_manager) : m_texture_manager(texture_manager) { initSampleJitter(); }

    // Sample a single grid using simple box filter.
    VolumeTexture sampleGrid(const openvdb::FloatGrid& grid, const openvdb::Coord& texture_extents);
    // Sample a multi res grid. Filtering is done by the built-in sampling mechanism of MultiResGrid.
    VolumeTexture sampleMultiResGrid(const openvdb::FloatGrid& grid, const openvdb::Coord& texture_extents);

private:
    void
    samplingLoop(float* output, const openvdb::CoordBBox& domain, std::function<void(openvdb::Coord,float&)> sampling_func);

    MHWRender::MTextureManager* m_texture_manager;

    static constexpr int N_FILT_SAMP = 1;
    void initSampleJitter();
    std::vector<openvdb::Vec3d> m_sample_jitter;
    std::vector<float> m_buffer;
};
