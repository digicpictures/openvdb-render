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

struct ChannelParams
{
    float intensity;
    MFloatVector color;
    std::string name;
    Gradient gradient;
    ChannelColorSource color_source;
    ChannelParams() : intensity(0) {}
    ChannelParams(const ChannelParams&) = default;
    ChannelParams& operator=(const ChannelParams&) = default;
    bool operator==(const ChannelParams& rhs) const
    {
        return intensity == rhs.intensity && color == rhs.color && name == rhs.name && !gradient.is_different(rhs.gradient);
    }
    bool operator!=(const ChannelParams& rhs) const { return !(*this == rhs); }
};

struct VDBVisualizerData {
    MBoundingBox bbox;

    std::string vdb_path;
    openvdb::io::File* vdb_file;

    ChannelParams density_channel;
    ChannelParams scattering_channel;
    ChannelParams attenuation_channel;
    ChannelParams emission_channel;
    ChannelParams transparency_channel;
    ChannelParams temperature_channel;
    float anisotropy;
    float blackbody_intensity;
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
    bool per_slice_gamma;

    VDBVisualizerData();
    ~VDBVisualizerData();

    void clear(const MBoundingBox& bb = MBoundingBox());
};
