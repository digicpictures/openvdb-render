#pragma once

#include <maya/MBoundingBox.h>
#include <maya/MFloatVector.h>
#include <openvdb/io/File.h>
#include "gradient.hpp"

enum VDBDisplayMode{
    DISPLAY_AXIS_ALIGNED_BBOX = 0,
    DISPLAY_GRID_BBOX,
    DISPLAY_POINT_CLOUD,
    DISPLAY_NON_SHADED,
    DISPLAY_SHADED,
    DISPLAY_MESH,
    DISPLAY_SLICES
};

struct VDBVisualizerData{
    MBoundingBox bbox;

    // Sliced display mode
    float value_min;
    float value_max;
    std::string density_channel;
    float slice_size;
    MFloatVector light_color;
    MFloatVector light_direction;
    MFloatVector scattering;
    MFloatVector absorption;
    float shadow_gain;
    int shadow_sample_count;

    // ===
    MFloatVector scattering_color;
    MFloatVector attenuation_color;
    MFloatVector emission_color;

    std::string vdb_path;
    std::string attenuation_channel;
    std::string scattering_channel;
    std::string emission_channel;

    Gradient scattering_gradient;
    Gradient attenuation_gradient;
    Gradient emission_gradient;

    openvdb::io::File* vdb_file;

    float point_size;
    float point_jitter;

    int point_skip;
    int update_trigger;
    VDBDisplayMode display_mode;

    VDBVisualizerData();
    ~VDBVisualizerData();

    void clear(const MBoundingBox& bb = MBoundingBox());
};
