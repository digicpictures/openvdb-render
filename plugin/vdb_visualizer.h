#pragma once

#include <maya/MPxSurfaceShape.h>
#include <maya/MPxSurfaceShapeUI.h>
#include <maya/MBoundingBox.h>

#include <openvdb/openvdb.h>
// std regex in gcc 4.8.3 is broken
#include <boost/regex.hpp>

#include "vdb_sampler.h"

enum VDBDisplayMode{
    DISPLAY_AXIS_ALIGNED_BBOX = 0,
    DISPLAY_GRID_BBOX,
    DISPLAY_POINT_CLOUD,
    DISPLAY_NON_SHADED,
    DISPLAY_SHADED,
    DISPLAY_MESH
};

struct VDBVisualizerData{
    MBoundingBox bbox;

    std::string vdb_path;

    openvdb::io::File* vdb_file;

    float point_size;
    float point_jitter;

    int update_trigger;
    VDBDisplayMode display_mode;

    VDBVisualizerData();
    ~VDBVisualizerData();

    void clear(const MBoundingBox& _bbox = MBoundingBox());
};

class VDBVisualizerShapeUI : public MPxSurfaceShapeUI {
public:
    VDBVisualizerShapeUI();
    ~VDBVisualizerShapeUI();

    static void* creator();

    virtual bool select(
            MSelectInfo& selectInfo,
            MSelectionList& selectionList,
            MPointArray& worldSpaceSelectPts) const;

    virtual bool canDrawUV() const;
};

class VDBVisualizerShape : public MPxSurfaceShape {
public:
    VDBVisualizerShape();
    ~VDBVisualizerShape();

    static void* creator();

    virtual bool isBounded() const;
    virtual MBoundingBox boundingBox() const;

    MStatus compute(const MPlug& plug, MDataBlock& dataBlock);
    static MStatus initialize();

    static const MTypeId typeId;
    static const MString typeName;
    static const MString drawDbClassification;

    static MObject s_vdb_path;
    static MObject s_cache_time;
    static MObject s_cache_playback_start;
    static MObject s_cache_playback_end;
    static MObject s_cache_playback_offset;
    static MObject s_cache_before_mode;
    static MObject s_cache_after_mode;
    static MObject s_display_mode;
    static MObject s_update_trigger;
    // mainly for 3rd party renderers
    static MObject s_out_vdb_path;
    static MObject s_grid_names;
    static MObject s_bbox_min;
    static MObject s_bbox_max;
    static MObject s_channel_stats;
    static MObject s_voxel_size;

    // display parameters
    static MObject s_point_size;
    static MObject s_point_jitter;

    // shader parameters
    static MObject s_scattering_source;
    static MObject s_scattering;
    static MObject s_scattering_channel;
    static MObject s_scattering_color;
    static MObject s_scattering_intensity;
    static MObject s_anisotropy;
    static MObject s_attenuation_source;
    static MObject s_attenuation;
    static MObject s_attenuation_channel;
    static MObject s_attenuation_color;
    static MObject s_attenuation_intensity;
    static MObject s_attenuation_mode;
    static MObject s_emission_source;
    static MObject s_emission;
    static MObject s_emission_channel;
    static MObject s_emission_color;
    static MObject s_emission_intensity;
    static MObject s_position_offset;
    static MObject s_interpolation;
    static MObject s_compensate_scaling;
    static MObject s_sampling_quality;

    static VDBGradientParams s_scattering_gradient;
    static VDBGradientParams s_attenuation_gradient;
    static VDBGradientParams s_emission_gradient;

    static MObject s_additional_channel_export;

    static const boost::regex s_frame_expr;
    static const boost::regex s_hash_expr;

    VDBVisualizerData* get_update();
private:
    VDBVisualizerData m_vdb_data;
};
