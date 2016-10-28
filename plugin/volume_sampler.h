#pragma once

#include <vector>

#include <openvdb/openvdb.h>
#include <openvdb/tools/MultiResGrid.h>

#include "util.h"

class ProgressBar;
namespace MHWRender {
    class MTexture;
}

struct VolumeTexture
{
    MHWRender::MTexture* texture;
    FloatRange value_range;

    VolumeTexture() : texture(nullptr) {}
    VolumeTexture(const openvdb::Coord& extents, const float* data_norm, const FloatRange& value_range_)
        : texture(acquireVolumeTexture(extents, data_norm)), value_range(value_range_) {}
    VolumeTexture(VolumeTexture&&) = default;
    VolumeTexture& operator=(VolumeTexture&&) = default;

private:
    MHWRender::MTexture* acquireVolumeTexture(const openvdb::Coord& texture_extents, const float* pixel_data);
};


// Sample a 3D domain (e.g. openvdb grid) using a sampling function (e.g. box filter),
// and create a Maya texture using the samples.
class VolumeSampler
{
public:
    // Generic sampling function.
    template <typename SamplingFunc>
    VolumeTexture sampleVolume(const openvdb::Coord& extents, SamplingFunc sampling_func, ProgressBar *progress_bar = nullptr);

    // Sample a single grid using simple box filter.
    VolumeTexture sampleGridWithBoxFilter(const openvdb::FloatGrid& grid, const openvdb::Coord& texture_extents, ProgressBar *progress_bar = nullptr);
    // Create a MultiResGrid and sample it. Filtering is done by the built-in sampling mechanism of MultiResGrid.
    VolumeTexture sampleGridWithMipmapFilter(const openvdb::FloatGrid& grid, const openvdb::Coord& texture_extents, ProgressBar *progress_bar = nullptr);

private:
    std::vector<float> m_buffer;
};
