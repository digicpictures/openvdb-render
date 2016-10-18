#pragma once

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
    VolumeTexture(const openvdb::Coord& extents, const float* data_norm, const ValueRange& value_range_)
        : texture(acquireVolumeTexture(extents, data_norm)), value_range(value_range_) {}
    VolumeTexture(VolumeTexture&&) = default;
    VolumeTexture& operator=(VolumeTexture&&) = default;

private:
    MHWRender::MTexture* acquireVolumeTexture(const openvdb::Coord& texture_extents, const float* pixel_data);
};

class VolumeSampler
{
public:
    // Sample a single grid using simple box filter.
    VolumeTexture sampleGrid(const openvdb::FloatGrid& grid, const openvdb::Coord& texture_extents);
    // Sample a multi res grid. Filtering is done by the built-in sampling mechanism of MultiResGrid.
    VolumeTexture sampleMultiResGrid(const openvdb::FloatGrid& grid, const openvdb::Coord& texture_extents);

private:
    void
    samplingLoop(float* output, const openvdb::CoordBBox& domain, std::function<void(openvdb::Coord,float&)> sampling_func);

    std::vector<float> m_buffer;
};
