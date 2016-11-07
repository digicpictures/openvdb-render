#pragma once

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

struct SliceShaderParams {
    MFloatVector scattering;
    MFloatVector absorption;
    float shadow_gain;
    int shadow_sample_count;
    bool operator==(const SliceShaderParams& rhs) const {
        return
            scattering == rhs.scattering &&
            absorption == rhs.absorption &&
            shadow_gain == rhs.shadow_gain &&
            shadow_sample_count == rhs.shadow_sample_count; 
    }
    bool operator !=(const SliceShaderParams& rhs) const { return !(*this == rhs); }
};

struct VDBVisualizerData {
    MBoundingBox bbox;

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

    float anisotropy;

    openvdb::io::File* vdb_file;

    float point_size;
    float point_jitter;

    int point_skip;
    int update_trigger;
    VDBDisplayMode display_mode;

    // Sliced display params.
    std::string sliced_display_channel;
    int max_slice_count;
    SliceShaderParams sliced_display_shader_params;

    VDBVisualizerData();
    ~VDBVisualizerData();

    void clear(const MBoundingBox& bb = MBoundingBox());
};
