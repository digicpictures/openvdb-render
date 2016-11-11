#pragma once

#include <utility>

#include <maya/MBoundingBox.h>
#include <maya/MFloatVector.h>
#include <openvdb/io/File.h>
#include "gradient.hpp"

enum VDBDisplayMode {
    DISPLAY_AXIS_ALIGNED_BBOX = 0,
    DISPLAY_GRID_BBOX,
    DISPLAY_POINT_CLOUD,
    DISPLAY_NON_SHADED,
    DISPLAY_SHADED,
    DISPLAY_MESH,
    DISPLAY_SLICES
};

enum class EmissionMode {
    NONE = 0,
    DENSITY = 1,
    CHANNEL = 2,
    BLACKBODY = 3,
    DENSITY_AND_BLACKBODY = 4
};

enum class ChannelColorSource {
    COLOR = 0,
    RAMP = 1
};

template <typename T>
struct ChannelParams
{
    typename T value_type;
    T value;
    std::string name;
    Gradient gradient;
    ChannelColorSource color_source;
    template <typename... Args>
    ChannelParams(Args&&... args) : value(std::forward<Args>(args)...) {}
    ChannelParams(const ChannelParams<T>&) = default;
    ChannelParams<T>& operator=(const ChannelParams<T>&) = default;
    bool operator==(const ChannelParams<T>& rhs) const
    {
        return value == rhs.value && name == rhs.name && !gradient.is_different(rhs.gradient);
    }
    bool operator!=(const ChannelParams<T>& rhs) const { return !(*this == rhs); }
};

struct VDBVisualizerData {
    MBoundingBox bbox;

    std::string vdb_path;
    openvdb::io::File* vdb_file;

    ChannelParams<float>        density_channel;
    ChannelParams<MFloatVector> scattering_channel;
    ChannelParams<MFloatVector> attenuation_channel;
    ChannelParams<MFloatVector> emission_channel;
    ChannelParams<MFloatVector> transparency_channel;
    ChannelParams<float>        temperature_channel;
    float anisotropy;
    EmissionMode emission_mode;

    int update_trigger;
    VDBDisplayMode display_mode;

    // Point cloud params.
    float point_size;
    float point_jitter;
    int point_skip;

    // Sliced display params.
    int max_slice_count;
    float shadow_gain;
    int shadow_sample_count;

    VDBVisualizerData();
    ~VDBVisualizerData();

    void clear(const MBoundingBox& bb = MBoundingBox());
};
