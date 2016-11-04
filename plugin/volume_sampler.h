#pragma once

#include <vector>

#include <openvdb/openvdb.h>
#include <openvdb/tools/MultiResGrid.h>

#include "vdb_subscene_utils.hpp"

class ProgressBar;
namespace MHWRender {
    class MTexture;
}

template <typename T>
struct ValueRange
{
    T min, max;
    ValueRange() noexcept : min(std::numeric_limits<T>::max()), max(std::numeric_limits<T>::min()) {}
    ValueRange(T min_ , T max_) noexcept : min(min_), max(max_) {}
    void update(T value) { min = std::min(min, value); max = std::max(max, value); }
    void update(const ValueRange<T>& value_range) { min = std::min(min, value_range.min); max = std::max(max, value_range.max); }
    ValueRange merge(const ValueRange<T>& value_range) { return { std::min(min, value_range.min), std::max(max, value_range.max) }; }
};
typedef ValueRange<float> FloatRange;

struct VolumeTexture
{
    TexturePtr texture;
    FloatRange value_range;

    VolumeTexture() : texture(nullptr) {}
    VolumeTexture(const openvdb::Coord& extents, const float* data_norm, const FloatRange& value_range_)
        : texture(acquireVolumeTexture(extents, data_norm)), value_range(value_range_) {}
    VolumeTexture(const VolumeTexture&) = delete;
    VolumeTexture& operator=(const VolumeTexture&) = delete;
    VolumeTexture(VolumeTexture&&) = default;
    VolumeTexture& operator=(VolumeTexture&&) = default;

    void clear() { texture.reset(); }
    operator bool() const { return texture.get() != nullptr; }

private:
    static MHWRender::MTexture* acquireVolumeTexture(const openvdb::Coord& texture_extents, const float* pixel_data);
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
    // Sample a MultiResGrid. Filtering is done by the built-in sampling mechanism of MultiResGrid.
    VolumeTexture sampleMultiResGrid(const openvdb::tools::MultiResGrid<openvdb::FloatTree>& multires, const openvdb::Coord& texture_extents, ProgressBar *progress_bar = nullptr);

private:
    std::vector<float> m_buffer;
};
