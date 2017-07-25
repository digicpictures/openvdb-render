#pragma once

#include <maya/MPxSurfaceShape.h>
#include <maya/MPxSurfaceShapeUI.h>
#include <maya/MPxCommand.h>
#include <maya/MSyntax.h>

#include <openvdb/openvdb.h>
#include <maya/MNodeMessage.h>
#include <maya/MDGMessage.h>

#include "vdb_visualizer_data.h"
#include "vdb_sampler.h"
#include "vdb_shader.h"
#include "vdb_simple_shader.h"
#include "vdb_std_shader.h"

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

    virtual void getDrawRequests(
            const MDrawInfo & info,
            bool objectAndActiveOnly,
            MDrawRequestQueue & requests);
    virtual void draw(
            const MDrawRequest & request,
            M3dView & view) const;

private:
    std::vector<MFloatVector> m_bbox_verts;
    std::vector<MFloatVector> m_cross_verts;
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

    void updateMaxSliceCount();

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

    // sliced display parameters
    static MObject s_max_slice_count;
    static MObject s_apply_max_slice_count;
    static MObject s_shadow_gain;
    static MObject s_shadow_sample_count;
    static MObject s_per_slice_gamma;

    // Velocity params
    static MObject s_velocity_grids;
    static MObject s_velocity_scale;
    static MObject s_velocity_fps;
    static MObject s_velocity_shutter_start;
    static MObject s_velocity_shutter_end;

    // shader parameters
    static VDBVolumeStandardShaderParams s_volume_standard_params;

    VDBVisualizerData* get_update();

    void setChannel(const MObject& channel_name, const MObject& intensity, const MObject& color, const MObject& color_source, const VDBGradientParams *gradient, ChannelParams& output);

private:

    VDBVisualizerData m_vdb_data;
    MCallbackId m_time_changed_id;
};

class VDBVisualizerUpdateMaxSliceCountCmd : public MPxCommand {
public:
    VDBVisualizerUpdateMaxSliceCountCmd() {}
    ~VDBVisualizerUpdateMaxSliceCountCmd() {}

    static void* creator() { return new VDBVisualizerUpdateMaxSliceCountCmd(); }
    static MSyntax create_syntax();

    MStatus doIt(const MArgList& args);
};
