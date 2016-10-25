#pragma once

#include <maya/MPxSurfaceShape.h>
#include <maya/MPxSurfaceShapeUI.h>

#include <openvdb/openvdb.h>
#include <maya/MNodeMessage.h>
#include <maya/MDGMessage.h>

#include "vdb_visualizer_data.h"
#include "vdb_sampler.h"
#include "vdb_shader.h"
#include "vdb_simple_shader.h"

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
    void postConstructor();

    static const MTypeId typeId;
    static const MString typeName;
    static const MString s_classification;

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
    static MObject s_matte;

    // sliced display parameters
    static MObject s_density_channel;
    static MObject s_slice_size;
    static MObject s_light_color;
    static MObject s_light_intensity;
    static MObject s_light_direction;
    static MObject s_scattering_color;
    static MObject s_scattering_intensity;
    static MObject s_absorption_color;
    static MObject s_absorption_intensity;
    static MObject s_shadow_gain;
    static MObject s_shadow_sample_count;

    // display parameters
    static MObject s_point_size;
    static MObject s_point_jitter;
    static MObject s_point_skip;

    static MObject s_override_shader;
    static MObject s_sampling_quality;
    static MObject s_additional_channel_export;

    // Velocity params
    static MObject s_velocity_grids;
    static MObject s_velocity_scale;
    static MObject s_velocity_fps;
    static MObject s_velocity_shutter_start;
    static MObject s_velocity_shutter_end;

    // disp params
    static MObject s_bounds_slack;

    // shader parameters
    static MObject s_shader_mode;
    static VDBShaderParams s_shader_params;
    static VDBSimpleShaderParams s_simple_shader_params;

    VDBVisualizerData* get_update();
private:
    static void attribute_changed(MNodeMessage::AttributeMessage msg, MPlug& plug, MPlug& other_plug, void* client_data);
    static void time_changed(MTime& time, void* client_data);

    VDBVisualizerData m_vdb_data;
    MCallbackId m_time_changed_id;
};
