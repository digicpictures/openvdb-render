/*
 * The update trigger is used to trigger updates for viewport visualization
 * if it's changed that means something have changed, and the viewport data has to be reloaded.
 *
 * This is done, so we only update the node if we are interactively editing the node, also
 * so we can completely separate reading the data from the main node, all it does it's just loading
 * a vdb dataset and reading information about the contained channels, but not loading actual voxel data.
 *
 * TODO
 * * Check if the metadata reads are cached by the file or not, if not then cache it manually
 */

#include "vdb_visualizer.h"

#include "../util/node_ids.h"

#include <maya/MArgDatabase.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MRampAttribute.h>
#include <maya/MTime.h>
#include <maya/MPointArray.h>
#include <maya/MSelectionList.h>
#include <maya/MDagPath.h>
#include <maya/MFnStringData.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MViewport2Renderer.h>

#include <array>
#include <sstream>
#include <maya/MGlobal.h>
#include <maya/MNodeMessage.h>

#include "gradient.hpp"
#include "vdb_maya_utils.hpp"

const MTypeId VDBVisualizerShape::typeId(ID_VDB_VISUALIZER);
const MString VDBVisualizerShape::typeName("vdb_visualizer");
const MString VDBVisualizerShape::drawDbClassification("drawdb/subscene/volume/vdb_visualizer");

MObject VDBVisualizerShape::s_update_trigger;
MObject VDBVisualizerShape::s_vdb_path;
MObject VDBVisualizerShape::s_cache_time;
MObject VDBVisualizerShape::s_cache_playback_start;
MObject VDBVisualizerShape::s_cache_playback_end;
MObject VDBVisualizerShape::s_cache_playback_offset;
MObject VDBVisualizerShape::s_cache_before_mode;
MObject VDBVisualizerShape::s_cache_after_mode;
MObject VDBVisualizerShape::s_display_mode;
MObject VDBVisualizerShape::s_out_vdb_path;
MObject VDBVisualizerShape::s_grid_names;
MObject VDBVisualizerShape::s_bbox_min;
MObject VDBVisualizerShape::s_bbox_max;
MObject VDBVisualizerShape::s_channel_stats;
MObject VDBVisualizerShape::s_voxel_size;
MObject VDBVisualizerShape::s_matte;

// Sliced display params.
MObject VDBVisualizerShape::s_max_slice_count;
MObject VDBVisualizerShape::s_apply_max_slice_count;
MObject VDBVisualizerShape::s_shadow_gain;
MObject VDBVisualizerShape::s_shadow_sample_count;
MObject VDBVisualizerShape::s_per_slice_gamma;

MObject VDBVisualizerShape::s_point_size;
MObject VDBVisualizerShape::s_point_jitter;
MObject VDBVisualizerShape::s_point_skip;

MObject VDBVisualizerShape::s_override_shader;
MObject VDBVisualizerShape::s_sampling_quality;
MObject VDBVisualizerShape::s_additional_channel_export;

MObject VDBVisualizerShape::s_velocity_grids;
MObject VDBVisualizerShape::s_velocity_scale;
MObject VDBVisualizerShape::s_velocity_fps;
MObject VDBVisualizerShape::s_velocity_shutter_start;
MObject VDBVisualizerShape::s_velocity_shutter_end;
MObject VDBVisualizerShape::s_bounds_slack;

MObject VDBVisualizerShape::s_shader_mode;
VDBShaderParams VDBVisualizerShape::s_shader_params;
VDBSimpleShaderParams VDBVisualizerShape::s_simple_shader_params;
VDBVolumeStandardShaderParams VDBVisualizerShape::s_volume_standard_params;

namespace {
    enum {
        CACHE_OUT_OF_RANGE_MODE_NONE = 0,
        CACHE_OUT_OF_RANGE_MODE_HOLD,
        CACHE_OUT_OF_RANGE_MODE_REPEAT
    };

    enum {
        SHADER_MODE_VOLUME_COLLECTOR = 0,
        SHADER_MODE_SIMPLE = 1,
        SHADER_MODE_VOLUME_STANDARD = 2
    };

    // we have to do lots of line, rectangle intersection, so using the Cohen-Sutherland algorithm
    // https://en.wikipedia.org/wiki/Cohen%E2%80%93Sutherland_algorithm
    class SelectionRectangle {
    private:
        float xmin, ymin;
        float xmax, ymax;

        enum {
            OUTCODE_INSIDE = 0, // 0000
            OUTCODE_LEFT = 1,   // 0001
            OUTCODE_RIGHT = 2,  // 0010
            OUTCODE_BOTTOM = 4, // 0100
            OUTCODE_TOP = 8,    // 1000
        };

        int compute_out_code(float x, float y) const
        {
            int code = OUTCODE_INSIDE;          // initialised as being inside of clip window

            if (x < xmin)           // to the left of clip window
                code |= OUTCODE_LEFT;
            else if (x > xmax)      // to the right of clip window
                code |= OUTCODE_RIGHT;
            if (y < ymin)           // below the clip window
                code |= OUTCODE_BOTTOM;
            else if (y > ymax)      // above the clip window
                code |= OUTCODE_TOP;

            return code;
        }
    public:
        SelectionRectangle(MSelectInfo& selectInfo)
        {
            unsigned int orig_x = 0;
            unsigned int orig_y = 0;
            unsigned int size_x = 0;
            unsigned int size_y = 0;
            selectInfo.selectRect(orig_x, orig_y, size_x, size_y);

            xmin = static_cast<float>(orig_x);
            ymin = static_cast<float>(orig_y);
            xmax = static_cast<float>(orig_x + size_x);
            ymax = static_cast<float>(orig_y + size_y);
        }

        float getXMin() const { return xmin; }
        float getXMax() const { return xmax; }
        float getYMin() const { return ymin; }
        float getYMax() const { return ymax; }

        bool clip_line(float x0, float y0, float x1, float y1) const
        {
            // compute outcodes for P0, P1, and whatever point lies outside the clip rectangle
            int outcode0 = compute_out_code(x0, y0);
            int outcode1 = compute_out_code(x1, y1);

            while (true)
            {
                if (!(outcode0 | outcode1)) // Bitwise OR is 0. Trivially accept and get out of loop
                    return true;
                else if (outcode0 & outcode1) // Bitwise AND is not 0. Trivially reject and get out of loop
                    return false;
                else
                {
                    float x = 0.0f;
                    float y = 0.0f;
                    const int outcode_out = outcode0 ? outcode0 : outcode1;
                    if (outcode_out & OUTCODE_TOP)
                    {
                        x = x0 + (x1 - x0) * (ymax - y0) / (y1 - y0);
                        y = ymax;
                    }
                    else if (outcode_out & OUTCODE_BOTTOM)
                    {
                        x = x0 + (x1 - x0) * (ymin - y0) / (y1 - y0);
                        y = ymin;
                    }
                    else if (outcode_out & OUTCODE_RIGHT)
                    {
                        y = y0 + (y1 - y0) * (xmax - x0) / (x1 - x0);
                        x = xmax;
                    }
                    else if (outcode_out & OUTCODE_LEFT)
                    {
                        y = y0 + (y1 - y0) * (xmin - x0) / (x1 - x0);
                        x = xmin;
                    }

                    if (outcode_out == outcode0)
                    {
                        x0 = x;
                        y0 = y;
                        outcode0 = compute_out_code(x0, y0);
                    }
                    else
                    {
                        x1 = x;
                        y1 = y;
                        outcode1 = compute_out_code(x1, y1);
                    }
                }
            }
        }
    };

    constexpr int DEFAULT_MAX_SLICE_COUNT = 128;
}

VDBVisualizerData::VDBVisualizerData()
    : bbox(MPoint(-1.0, -1.0, -1.0), MPoint(1.0, 1.0, 1.0)),
    vdb_file(nullptr),
    point_size(std::numeric_limits<float>::infinity()),
    point_jitter(std::numeric_limits<float>::infinity()),
    point_skip(-1),
    update_trigger(-1),
    max_slice_count(-1),
    display_mode(DISPLAY_AXIS_ALIGNED_BBOX)
{
}

VDBVisualizerData::~VDBVisualizerData()
{
    clear();
}

void VDBVisualizerData::clear(const MBoundingBox& bb)
{
    if (vdb_file != nullptr)
    {
        if (vdb_file->isOpen())
            vdb_file->close();
        delete vdb_file;
        vdb_file = nullptr;
    }
    bbox = bb;
}

VDBVisualizerShape::VDBVisualizerShape()
{

}

VDBVisualizerShape::~VDBVisualizerShape()
{
    if (MGlobal::mayaState() == MGlobal::kInteractive)
        MDGMessage::removeCallback(m_time_changed_id);
}

void* VDBVisualizerShape::creator()
{
    return new VDBVisualizerShape();
}

MStatus VDBVisualizerShape::compute(const MPlug& plug, MDataBlock& dataBlock)
{
    MStatus status = MS::kSuccess;

    if (plug == s_out_vdb_path)
    {
        MayaPathSpec vdb_path_spec = dataBlock.inputValue(s_vdb_path).asString().asChar();
        std::string vdb_path = vdb_path_spec.getPath();

        if (vdb_path_spec.hasFrameField())
        {
            const double cache_time = dataBlock.inputValue(s_cache_time).asTime().as(MTime::uiUnit());
            const double cache_playback_offset = dataBlock.inputValue(s_cache_playback_offset).asTime().as(MTime::uiUnit());
            int cache_frame = static_cast<int>(cache_time - cache_playback_offset);
            const int cache_playback_start = dataBlock.inputValue(s_cache_playback_start).asInt();
            const int cache_playback_end = std::max(cache_playback_start, dataBlock.inputValue(s_cache_playback_end).asInt());
            bool frame_in_range = true;
            if (cache_frame < cache_playback_start)
            {
                const short cache_before_mode = dataBlock.inputValue(s_cache_before_mode).asShort();
                if (cache_before_mode == CACHE_OUT_OF_RANGE_MODE_NONE)
                    frame_in_range = false;
                else if (cache_before_mode == CACHE_OUT_OF_RANGE_MODE_HOLD)
                    cache_frame = cache_playback_start;
                else if (cache_before_mode == CACHE_OUT_OF_RANGE_MODE_REPEAT)
                {
                    const int cache_playback_range = cache_playback_end - cache_playback_start;
                    cache_frame = cache_playback_end - (cache_playback_start - cache_frame - 1) % (cache_playback_range + 1);
                }
            }
            else if (cache_frame > cache_playback_end)
            {
                const short cache_after_mode = dataBlock.inputValue(s_cache_after_mode).asShort();
                if (cache_after_mode == CACHE_OUT_OF_RANGE_MODE_NONE)
                    frame_in_range = false;
                else if (cache_after_mode == CACHE_OUT_OF_RANGE_MODE_HOLD)
                    cache_frame = cache_playback_end;
                else if (cache_after_mode == CACHE_OUT_OF_RANGE_MODE_REPEAT)
                {
                    const int cache_playback_range = cache_playback_end - cache_playback_start;
                    cache_frame = cache_playback_start + (cache_frame - cache_playback_end - 1) % (cache_playback_range + 1);
                }
            }
            cache_frame = std::max(0, cache_frame);

            if (frame_in_range) {
                vdb_path = vdb_path_spec.getPath(cache_frame);
            } else {
                vdb_path = "";
            }
        }

        if (vdb_path != m_vdb_data.vdb_path)
        {
            m_vdb_data.max_slice_count = dataBlock.inputValue(s_max_slice_count).asInt();

            m_vdb_data.clear();
            m_vdb_data.vdb_path = vdb_path;
            if (m_vdb_data.vdb_file != nullptr)
                delete m_vdb_data.vdb_file;
            try{
                m_vdb_data.vdb_file = new openvdb::io::File(vdb_path.c_str());
                m_vdb_data.vdb_file->open(false);
                if (m_vdb_data.vdb_file->isOpen())
                {
                    openvdb::GridPtrVecPtr grids = m_vdb_data.vdb_file->readAllGridMetadata();
                    for (openvdb::GridPtrVec::const_iterator it = grids->begin(); it != grids->end(); ++it)
                    {
                        if (openvdb::GridBase::ConstPtr grid = *it)
                            read_transformed_bounding_box(grid, m_vdb_data.bbox);
                    }
                }
                else
                    m_vdb_data.clear(MBoundingBox(MPoint(-1.0, -1.0, -1.0), MPoint(1.0, 1.0, 1.0)));
            }
            catch(...)
            {
                m_vdb_data.clear(MBoundingBox(MPoint(-1.0, -1.0, -1.0), MPoint(1.0, 1.0, 1.0)));
            }
        }
        MDataHandle out_vdb_path_handle = dataBlock.outputValue(s_out_vdb_path);
        out_vdb_path_handle.setString(vdb_path.c_str());
    }
    else
    {
        dataBlock.inputValue(s_out_vdb_path).asString(); // trigger cache reload
        if (plug == s_grid_names)
        {
            MDataHandle grid_names_handle = dataBlock.outputValue(s_grid_names);
            if (m_vdb_data.vdb_file != nullptr && m_vdb_data.vdb_file->isOpen())
            {
                std::stringstream grid_names;
                openvdb::GridPtrVecPtr grids = m_vdb_data.vdb_file->readAllGridMetadata();
                for (openvdb::GridPtrVec::const_iterator it = grids->begin(); it != grids->end(); ++it)
                {
                    if (openvdb::GridBase::ConstPtr grid = *it)
                        grid_names << grid->getName() << " ";
                }
                std::string grid_names_string = grid_names.str();
                grid_names_handle.setString(grid_names_string.length() ? grid_names_string.substr(0, grid_names_string.length() - 1).c_str() : "");
            }
            else
                grid_names_handle.setString("");
        }
        else if (plug == s_update_trigger)
        {
            MDataHandle update_trigger_handle = dataBlock.outputValue(s_update_trigger);
            update_trigger_handle.setInt(m_vdb_data.update_trigger + 1);
        }
        else if (plug == s_bbox_min)
        {
            // TODO : why the MDataBlock is buggy in this case?
            const MPoint mn = m_vdb_data.bbox.min();
            plug.child(0).setDouble(mn.x);
            plug.child(1).setDouble(mn.y);
            plug.child(2).setDouble(mn.z);
        }
        else if (plug == s_bbox_max)
        {
            // TODO : why the MDataBlock is buggy in this case?
            const MPoint mx = m_vdb_data.bbox.max();
            plug.child(0).setDouble(mx.x);
            plug.child(1).setDouble(mx.y);
            plug.child(2).setDouble(mx.z);
        }
        else if (plug == s_channel_stats)
        {
            std::stringstream ss;
            if (m_vdb_data.vdb_file != nullptr && m_vdb_data.vdb_file->isOpen())
            {
                ss << "Bounding box : " << "[ [";
                ss << m_vdb_data.bbox.min().x << ", " << m_vdb_data.bbox.min().y << ", " << m_vdb_data.bbox.min().z;
                ss << " ] [ ";
                ss << m_vdb_data.bbox.max().x << ", " << m_vdb_data.bbox.max().y << ", " << m_vdb_data.bbox.max().z;
                ss << " ] ]" << std::endl;
                ss << "Channels : " << std::endl;
                openvdb::GridPtrVecPtr grids = m_vdb_data.vdb_file->readAllGridMetadata();
                for (openvdb::GridPtrVec::const_iterator it = grids->begin(); it != grids->end(); ++it)
                {
                    if (openvdb::GridBase::ConstPtr grid = *it)
                        ss << " - " << grid->getName() << " (" << grid->valueType() << ")" << std::endl;
                }
            }
            dataBlock.outputValue(s_channel_stats).setString(ss.str().c_str());
        }
        else if (plug == s_voxel_size)
        {
            float voxel_size = std::numeric_limits<float>::max();
            if (m_vdb_data.vdb_file != nullptr && m_vdb_data.vdb_file->isOpen())
            {
                openvdb::GridPtrVecPtr grids = m_vdb_data.vdb_file->readAllGridMetadata();
                for (openvdb::GridPtrVec::const_iterator it = grids->begin(); it != grids->end(); ++it)
                {
                    if (openvdb::GridBase::ConstPtr grid = *it)
                    {
                        openvdb::Vec3d vs = grid->voxelSize();
                        if (vs.x() > 0.0)
                            voxel_size = std::min(static_cast<float>(vs.x()), voxel_size);
                        if (vs.y() > 0.0)
                            voxel_size = std::min(static_cast<float>(vs.y()), voxel_size);
                        if (vs.z() > 0.0)
                            voxel_size = std::min(static_cast<float>(vs.z()), voxel_size);
                    }
                }
            }
            else
                voxel_size = 1.0f;
            dataBlock.outputValue(s_voxel_size).setFloat(voxel_size);
        }
        else
            return MStatus::kUnknownParameter;
    }

    dataBlock.setClean(plug);

    return status;
}

MStatus VDBVisualizerShape::initialize()
{
    MFnNumericAttribute nAttr;
    MFnTypedAttribute tAttr;
    MFnEnumAttribute eAttr;
    MFnUnitAttribute uAttr;
    MFnStringData sData;

    MStatus status = MS::kSuccess;

    // input params

    s_vdb_path = tAttr.create("vdbPath", "vdb_path", MFnData::kString);

    s_cache_time = uAttr.create("cacheTime", "cache_time", MFnUnitAttribute::kTime, 0.0);

    s_cache_playback_start = nAttr.create("cachePlaybackStart", "cache_playback_start", MFnNumericData::kInt);
    nAttr.setDefault(1);
    nAttr.setMin(0);

    s_cache_playback_end = nAttr.create("cachePlaybackEnd", "cache_playback_end", MFnNumericData::kInt);
    nAttr.setDefault(100);
    nAttr.setMin(0);

    s_cache_playback_offset = uAttr.create("cachePlaybackOffset", "cache_playback_offset", MFnUnitAttribute::kTime, 0.0);

    s_cache_before_mode = eAttr.create("cacheBeforeMode", "cache_before_mode");
    eAttr.addField("None", CACHE_OUT_OF_RANGE_MODE_NONE);
    eAttr.addField("Hold", CACHE_OUT_OF_RANGE_MODE_HOLD);
    eAttr.addField("Repeat", CACHE_OUT_OF_RANGE_MODE_REPEAT);
    eAttr.setDefault(CACHE_OUT_OF_RANGE_MODE_HOLD);

    s_cache_after_mode = eAttr.create("cacheAfterMode", "cache_after_mode");
    eAttr.addField("None", CACHE_OUT_OF_RANGE_MODE_NONE);
    eAttr.addField("Hold", CACHE_OUT_OF_RANGE_MODE_HOLD);
    eAttr.addField("Repeat", CACHE_OUT_OF_RANGE_MODE_REPEAT);
    eAttr.setDefault(CACHE_OUT_OF_RANGE_MODE_HOLD);

    s_display_mode = eAttr.create("displayMode", "display_mode");
    eAttr.addField("Axis Aligned Bounding Box", DISPLAY_AXIS_ALIGNED_BBOX);
    eAttr.addField("Per Grid Bounding Box", DISPLAY_GRID_BBOX);
    eAttr.addField("Point Cloud", DISPLAY_POINT_CLOUD);
    eAttr.addField("Volumetric Non Shaded", DISPLAY_NON_SHADED);
    eAttr.addField("Volumetric Shaded", DISPLAY_SHADED);
    eAttr.addField("Mesh", DISPLAY_MESH);
    eAttr.addField("Slices", DISPLAY_SLICES);
    eAttr.setDefault(DISPLAY_SLICES);

    addAttribute(s_display_mode);

    // output params

    s_update_trigger = nAttr.create("updateTrigger", "update_trigger", MFnNumericData::kInt);
    nAttr.setDefault(0);
    nAttr.setStorable(false);
    nAttr.setReadable(false);
    nAttr.setWritable(false);

    s_grid_names = tAttr.create("gridNames", "grid_names", MFnData::kString);
    tAttr.setStorable(false);
    tAttr.setWritable(false);
    tAttr.setReadable(true);

    s_out_vdb_path = tAttr.create("outVdbPath", "out_vdb_path", MFnData::kString);
    tAttr.setStorable(false);
    tAttr.setWritable(false);
    tAttr.setReadable(true);

    s_bbox_min = nAttr.createPoint("bboxMin", "bbox_min");
    nAttr.setStorable(false);
    nAttr.setWritable(false);
    nAttr.setReadable(true);

    s_bbox_max = nAttr.createPoint("bboxMax", "bbox_max");
    nAttr.setStorable(false);
    nAttr.setWritable(false);
    nAttr.setReadable(true);

    s_channel_stats = tAttr.create("channelStats", "channel_stats", MFnData::kString);
    tAttr.setStorable(false);
    tAttr.setWritable(false);
    tAttr.setReadable(true);

    s_voxel_size = nAttr.create("voxelSize", "voxel_size", MFnNumericData::kFloat);
    nAttr.setStorable(false);
    nAttr.setWritable(false);
    nAttr.setReadable(true);

    MObject input_params[] = {
        s_vdb_path, s_cache_time, s_cache_playback_start, s_cache_playback_end,
        s_cache_playback_offset, s_cache_before_mode, s_cache_after_mode
    };

    MObject output_params[] = {
        s_update_trigger, s_grid_names, s_out_vdb_path, s_bbox_min, s_bbox_max, s_channel_stats, s_voxel_size
    };

    for (const auto& output_param : output_params)
        addAttribute(output_param);

    for (const auto& input_param : input_params)
    {
        addAttribute(input_param);

        for (const auto& output_param : output_params)
            attributeAffects(input_param, output_param);
    }

    attributeAffects(s_display_mode, s_update_trigger);

    s_matte = nAttr.create("matte", "matte", MFnNumericData::kBoolean);
    nAttr.setDefault(false);
    addAttribute(s_matte);

    s_point_size = nAttr.create("pointSize", "point_size", MFnNumericData::kFloat);
    nAttr.setMin(1.0f);
    nAttr.setSoftMax(10.0f);
    nAttr.setDefault(2.0f);
    nAttr.setChannelBox(true);

    s_point_jitter = nAttr.create("pointJitter", "point_jitter", MFnNumericData::kFloat);
    nAttr.setMin(0.0f);
    nAttr.setSoftMax(0.5f);
    nAttr.setDefault(0.15f);

    s_point_skip = nAttr.create("pointSkip", "point_skip", MFnNumericData::kInt);
    nAttr.setMin(1);
    nAttr.setSoftMax(20);
    nAttr.setDefault(10);
    nAttr.setChannelBox(true);

    // Sliced display attributes.

    s_max_slice_count = nAttr.create("maxSliceCount", "max_slice_count", MFnNumericData::kInt);
    nAttr.setDefault(DEFAULT_MAX_SLICE_COUNT);
    nAttr.setMin(4);
    nAttr.setSoftMax(256);
    addAttribute(s_max_slice_count);

    s_apply_max_slice_count = nAttr.create("applyMaxSliceCount", "apply_max_slice_count", MFnNumericData::kBoolean);
    nAttr.setDefault(false);
    addAttribute(s_apply_max_slice_count);
    attributeAffects(s_apply_max_slice_count, s_update_trigger);

    s_shadow_gain = nAttr.create("shadowGain", "shadow_gain", MFnNumericData::kFloat);
    nAttr.setDefault(0.2);
    nAttr.setMin(0.0);
    nAttr.setSoftMax(4.0);
    addAttribute(s_shadow_gain);
    attributeAffects(s_shadow_gain, s_update_trigger);

    s_shadow_sample_count = nAttr.create("shadowSampleCount", "shadow_sample_count", MFnNumericData::kInt);
    nAttr.setDefault(8);
    nAttr.setMin(0);
    nAttr.setSoftMax(16);
    addAttribute(s_shadow_sample_count);
    attributeAffects(s_shadow_sample_count, s_update_trigger);

    s_per_slice_gamma = nAttr.create("perSliceGamma", "per_slice_gamma", MFnNumericData::kBoolean);
    nAttr.setDefault(false);
    addAttribute(s_per_slice_gamma);
    attributeAffects(s_per_slice_gamma, s_update_trigger);

    //

    s_override_shader = nAttr.create("overrideShader", "override_shader", MFnNumericData::kBoolean);
    nAttr.setDefault(false);
    nAttr.setChannelBox(true);

    s_sampling_quality = nAttr.create("samplingQuality", "sampling_quality", MFnNumericData::kFloat);
    nAttr.setDefault(100.0f);
    nAttr.setMin(1.0f);
    nAttr.setSoftMax(300.0f);
    nAttr.setChannelBox(true);
    addAttribute(s_sampling_quality);

    s_additional_channel_export = tAttr.create("additionalChannelExport", "additional_channel_export", MFnData::kString);
    addAttribute(s_additional_channel_export);

    s_velocity_grids = tAttr.create("velocityGrids", "velocity_grids", MFnData::kString);
    addAttribute(s_velocity_grids);

    s_velocity_scale = nAttr.create("velocityScale", "velocity_scale", MFnNumericData::kFloat);
    nAttr.setSoftMin(0.0f);
    nAttr.setSoftMax(2.0f);
    nAttr.setDefault(1.0f);
    nAttr.setChannelBox(true);
    addAttribute(s_velocity_scale);

    s_velocity_fps = nAttr.create("velocityFps", "velocity_fps", MFnNumericData::kFloat);
    nAttr.setMin(0.0f);
    nAttr.setSoftMax(30.0f);
    nAttr.setDefault(24.0f);
    addAttribute(s_velocity_fps);

    s_velocity_shutter_start = nAttr.create("velocityShutterStart", "velocity_shutter_start", MFnNumericData::kFloat);
    nAttr.setDefault(-0.25f);
    nAttr.setSoftMin(-1.0f);
    nAttr.setSoftMax(0.0f);
    addAttribute(s_velocity_shutter_start);

    s_velocity_shutter_end = nAttr.create("velocityShutterEnd", "velocity_shutter_end", MFnNumericData::kFloat);
    nAttr.setDefault(0.25f);
    nAttr.setSoftMin(0.0f);
    nAttr.setSoftMax(1.0f);
    addAttribute(s_velocity_shutter_end);

    s_bounds_slack = nAttr.create("boundsSlack", "bounds_slack", MFnNumericData::kFloat);
    nAttr.setDefault(0.0f);
    nAttr.setMin(0.0f);
    nAttr.setSoftMax(1.0f);
    addAttribute(s_bounds_slack);

    s_shader_mode = eAttr.create("shaderMode", "shader_mode");
    eAttr.addField("Arnold Volume Collector", SHADER_MODE_VOLUME_COLLECTOR);
    eAttr.addField("Simple Shader", SHADER_MODE_SIMPLE);
    eAttr.addField("Volume Standard Shader", SHADER_MODE_VOLUME_STANDARD);
    eAttr.setDefault(2);
    addAttribute(s_shader_mode);

    s_shader_params.create_params();
    s_simple_shader_params.create_params(false);
    s_volume_standard_params.create_params(&s_shader_params);

    MObject display_params[] = {
            s_point_size, s_point_jitter, s_point_skip, s_override_shader, s_shader_mode
    };

    for (const auto& shader_param : display_params)
    {
        addAttribute(shader_param);
        attributeAffects(shader_param, s_update_trigger);
    }

    s_shader_params.affect_output(s_update_trigger);
    s_simple_shader_params.affect_output(s_update_trigger);
    s_volume_standard_params.affect_output(s_update_trigger);

    return status;
}

namespace {

    MFloatVector plugAsFloatVector(const MPlug plug)
    {
        return { plug.child(0).asFloat(), plug.child(1).asFloat(), plug.child(2).asFloat() };
    }

} // unnamed namespace

void VDBVisualizerShape::setChannel(const MObject& channel_name, const MObject& intensity, const MObject& color, const MObject& color_source, const VDBGradientParams *gradient, ChannelParams& output)
{
    MObject tmo = thisMObject();
    output.name = MPlug(tmo, channel_name).asString().asChar();
    if (!intensity.isNull())
        output.intensity = MPlug(tmo, intensity).asFloat();
    if (!color.isNull())
        output.color = LinearFromSRGB(plugAsFloatVector(MPlug(tmo, color)));
    if (!color_source.isNull())
        output.color_source = ChannelColorSource(MPlug(tmo, color_source).asInt());
    if (gradient)
        output.gradient.update(*gradient, tmo);
}

VDBVisualizerData* VDBVisualizerShape::get_update()
{
    const int update_trigger = MPlug(thisMObject(), s_update_trigger).asInt();

    if (update_trigger != m_vdb_data.update_trigger)
    {
        MObject tmo = thisMObject();
        auto display_mode = static_cast<VDBDisplayMode>(MPlug(tmo, s_display_mode).asShort());
        if (display_mode != m_vdb_data.display_mode)
        {
            if (display_mode == DISPLAY_SLICES)
            {
                // Always use RGB ramp in sliced display mode.
                m_vdb_data.scattering_channel.gradient.set_channel_mode_override(4);
                m_vdb_data.emission_channel.gradient.set_channel_mode_override(4);
            }
            else
            {
                m_vdb_data.scattering_channel.gradient.clear_channel_mode_override();
                m_vdb_data.emission_channel.gradient.clear_channel_mode_override();
            }
            m_vdb_data.display_mode = display_mode;
        }

        m_vdb_data.point_size = MPlug(tmo, s_point_size).asFloat();
        m_vdb_data.point_jitter = MPlug(tmo, s_point_jitter).asFloat();
        m_vdb_data.point_skip = MPlug(tmo, s_point_skip).asInt();
        m_vdb_data.update_trigger = update_trigger;

        if (m_vdb_data.display_mode < DISPLAY_POINT_CLOUD)
        {
            return &m_vdb_data;
        }
        else if (m_vdb_data.display_mode == DISPLAY_SLICES)
        {
            m_vdb_data.shadow_gain =         MPlug(tmo, s_shadow_gain).asFloat();
            m_vdb_data.shadow_sample_count = MPlug(tmo, s_shadow_sample_count).asInt();
            m_vdb_data.per_slice_gamma =     MPlug(tmo, s_per_slice_gamma).asFloat();
        }

        if (m_vdb_data.display_mode >= DISPLAY_POINT_CLOUD)
        {
            const short shader_mode = MPlug(tmo, s_shader_mode).asShort();
            if (shader_mode == SHADER_MODE_VOLUME_COLLECTOR)
            {
                const short scattering_mode = MPlug(tmo, s_shader_params.scattering_source).asShort();
                MPlug scattering_color_plug(tmo, s_shader_params.scattering_color);
                const float scattering_intensity = MPlug(tmo, s_shader_params.scattering_intensity).asFloat();
                m_vdb_data.scattering_channel.color.x = scattering_color_plug.child(0).asFloat() * scattering_intensity;
                m_vdb_data.scattering_channel.color.y = scattering_color_plug.child(1).asFloat() * scattering_intensity;
                m_vdb_data.scattering_channel.color.z = scattering_color_plug.child(2).asFloat() * scattering_intensity;

                if (scattering_mode == 1)
                    m_vdb_data.scattering_channel.name = MPlug(tmo, s_shader_params.scattering_channel).asString().asChar();
                else
                {
                    MPlug scattering_plug(tmo, s_shader_params.scattering);
                    if (!scattering_plug.isConnected()) // TODO: handle this
                    {
                        m_vdb_data.scattering_channel.color.x *= scattering_plug.child(0).asFloat();
                        m_vdb_data.scattering_channel.color.y *= scattering_plug.child(1).asFloat();
                        m_vdb_data.scattering_channel.color.z *= scattering_plug.child(2).asFloat();
                    }

                    m_vdb_data.scattering_channel.name = "";
                }

                const short attenuation_mode = MPlug(tmo, s_shader_params.attenuation_source).asShort();
                MPlug attenuation_color_plug(tmo, s_shader_params.attenuation_color);
                const float attenuation_intensity = MPlug(tmo, s_shader_params.attenuation_intensity).asFloat();
                m_vdb_data.attenuation_channel.color.x = attenuation_color_plug.child(0).asFloat() * attenuation_intensity;
                m_vdb_data.attenuation_channel.color.y = attenuation_color_plug.child(1).asFloat() * attenuation_intensity;
                m_vdb_data.attenuation_channel.color.z = attenuation_color_plug.child(2).asFloat() * attenuation_intensity;

                if (attenuation_mode == 1)
                    m_vdb_data.attenuation_channel.name = MPlug(tmo, s_shader_params.attenuation_channel).asString().asChar();
                else if (attenuation_mode == 0)
                {
                    MPlug attenuation_plug(tmo, s_shader_params.attenuation);
                    if (!attenuation_plug.isConnected()) // TODO: handle this
                    {
                        m_vdb_data.attenuation_channel.color.x *= attenuation_plug.child(0).asFloat();
                        m_vdb_data.attenuation_channel.color.y *= attenuation_plug.child(1).asFloat();
                        m_vdb_data.attenuation_channel.color.z *= attenuation_plug.child(2).asFloat();
                    }

                    m_vdb_data.attenuation_channel.name = "";
                }
                else
                    m_vdb_data.attenuation_channel.name = m_vdb_data.scattering_channel.name;

                const short emission_mode = MPlug(tmo, s_shader_params.emission_source).asShort();
                MPlug emission_color_plug(tmo, s_shader_params.emission_color);
                const float emission_intensity = MPlug(tmo, s_shader_params.emission_intensity).asFloat();
                m_vdb_data.emission_channel.color.x = emission_color_plug.child(0).asFloat() * emission_intensity;
                m_vdb_data.emission_channel.color.y = emission_color_plug.child(1).asFloat() * emission_intensity;
                m_vdb_data.emission_channel.color.z = emission_color_plug.child(2).asFloat() * emission_intensity;

                if (emission_mode == 1)
                    m_vdb_data.emission_channel.name = MPlug(tmo, s_shader_params.emission_channel).asString().asChar();
                else
                {
                    MPlug emission_plug(tmo, s_shader_params.emission);
                    if (!emission_plug.isConnected()) // TODO: handle this
                    {
                        m_vdb_data.emission_channel.color.x *= emission_plug.child(0).asFloat();
                        m_vdb_data.emission_channel.color.y *= emission_plug.child(1).asFloat();
                        m_vdb_data.emission_channel.color.z *= emission_plug.child(2).asFloat();
                    }

                    m_vdb_data.emission_channel.name = "";
                }

                m_vdb_data.scattering_channel.gradient.update(s_shader_params.scattering_gradient, tmo);
                m_vdb_data.attenuation_channel.gradient.update(s_shader_params.attenuation_gradient, tmo);
                m_vdb_data.emission_channel.gradient.update(s_shader_params.emission_gradient, tmo);

                m_vdb_data.anisotropy = MPlug(tmo, s_shader_params.anisotropy).asFloat();
            }
            else if (shader_mode == SHADER_MODE_SIMPLE)
            {
                m_vdb_data.scattering_channel.name = MPlug(tmo, s_simple_shader_params.smoke_channel).asString().asChar();
                const float smoke_intensity = MPlug(tmo, s_simple_shader_params.smoke_intensity).asFloat();
                MPlug smoke_plug(tmo, s_simple_shader_params.smoke);
                if (!smoke_plug.isConnected())
                {
                    m_vdb_data.scattering_channel.color.x = smoke_plug.child(0).asFloat() * smoke_intensity;
                    m_vdb_data.scattering_channel.color.y = smoke_plug.child(1).asFloat() * smoke_intensity;
                    m_vdb_data.scattering_channel.color.z = smoke_plug.child(2).asFloat() * smoke_intensity;
                }
                else
                    m_vdb_data.scattering_channel.color.x = m_vdb_data.scattering_channel.color.y = m_vdb_data.scattering_channel.color.z = smoke_intensity;

                m_vdb_data.attenuation_channel.name = MPlug(tmo, s_simple_shader_params.opacity_channel).asString().asChar();
                const float opacity_intensity = MPlug(tmo, s_simple_shader_params.opacity_intensity).asFloat();
                MPlug opacity_plug(tmo, s_simple_shader_params.opacity);
                if (!opacity_plug.isConnected())
                {
                    m_vdb_data.attenuation_channel.color.x = opacity_plug.child(0).asFloat() * opacity_intensity;
                    m_vdb_data.attenuation_channel.color.y = opacity_plug.child(1).asFloat() * opacity_intensity;
                    m_vdb_data.attenuation_channel.color.z = opacity_plug.child(2).asFloat() * opacity_intensity;
                }
                else
                    m_vdb_data.attenuation_channel.color.x = m_vdb_data.attenuation_channel.color.y = m_vdb_data.attenuation_channel.color.z = opacity_intensity;

                m_vdb_data.emission_channel.name = MPlug(tmo, s_simple_shader_params.fire_channel).asString().asChar();
                const float fire_intensity = MPlug(tmo, s_simple_shader_params.fire_intensity).asFloat();
                MPlug fire_plug(tmo, s_simple_shader_params.fire);
                if (!fire_plug.isConnected())
                {
                    m_vdb_data.emission_channel.color.x = fire_plug.child(0).asFloat() * fire_intensity;
                    m_vdb_data.emission_channel.color.y = fire_plug.child(1).asFloat() * fire_intensity;
                    m_vdb_data.emission_channel.color.z = fire_plug.child(2).asFloat() * fire_intensity;
                }
                else
                    m_vdb_data.emission_channel.color.x = m_vdb_data.emission_channel.color.y = m_vdb_data.emission_channel.color.z = fire_intensity;

                m_vdb_data.scattering_channel.gradient.clear();
                m_vdb_data.attenuation_channel.gradient.clear();
                m_vdb_data.emission_channel.gradient.clear();

                m_vdb_data.anisotropy = MPlug(tmo, s_shader_params.anisotropy).asFloat();
            }
            else if (shader_mode == SHADER_MODE_VOLUME_STANDARD)
            {
                // Density.
                setChannel(s_volume_standard_params.density_channel,
                           s_volume_standard_params.density,
                           {}, {}, nullptr,
                           m_vdb_data.density_channel);

                // Scattering.
                setChannel(s_volume_standard_params.scattering_channel,
                           s_volume_standard_params.scattering_intensity,
                           s_volume_standard_params.scattering_color,
                           s_volume_standard_params.scattering_source,
                           &s_volume_standard_params.scattering_gradient,
                           m_vdb_data.scattering_channel);
                m_vdb_data.anisotropy = MPlug(tmo, s_volume_standard_params.anisotropy).asFloat();

                // Transparency.
                setChannel(s_volume_standard_params.transparency_channel,
                           {}, s_volume_standard_params.transparency,
                           {}, nullptr,
                           m_vdb_data.transparency_channel);

                // Emission.
                setChannel(s_volume_standard_params.emission_channel,
                           s_volume_standard_params.emission_intensity,
                           s_volume_standard_params.emission_color,
                           s_volume_standard_params.emission_source,
                           &s_volume_standard_params.emission_gradient,
                           m_vdb_data.emission_channel);
                m_vdb_data.emission_mode = EmissionMode(MPlug(tmo, s_volume_standard_params.emission_mode).asShort());

                // Temperature.
                setChannel(s_volume_standard_params.temperature_channel,
                           s_volume_standard_params.temperature,
                           {}, {}, nullptr,
                           m_vdb_data.temperature_channel);
                m_vdb_data.blackbody_intensity = MPlug(tmo, s_volume_standard_params.blackbody_intensity).asFloat();
            }
        }
    }
    return &m_vdb_data; // this data will be checked in the subscene override
}

bool VDBVisualizerShape::isBounded() const
{
    return true;
}

MBoundingBox VDBVisualizerShape::boundingBox() const
{
    MPlug(thisMObject(), s_out_vdb_path).asString();
    return m_vdb_data.bbox;
}

void VDBVisualizerShape::postConstructor()
{
    setRenderable(true);
    MObject tmo = thisMObject();
    s_shader_params.scattering_gradient.post_constructor(tmo);
    s_shader_params.attenuation_gradient.post_constructor(tmo);
    s_shader_params.emission_gradient.post_constructor(tmo);
}

VDBVisualizerShapeUI::VDBVisualizerShapeUI()
{

}

VDBVisualizerShapeUI::~VDBVisualizerShapeUI()
{
}

void* VDBVisualizerShapeUI::creator()
{
    return new VDBVisualizerShapeUI();
}

bool VDBVisualizerShapeUI::select(MSelectInfo& selectInfo, MSelectionList& selectionList, MPointArray& worldSpaceSelectPts) const
{
    if (selectInfo.isRay())
        return false;

    const MBoundingBox bbox = surfaceShape()->boundingBox();

    const MPoint min = bbox.min();
    const MPoint max = bbox.max();

    M3dView view = selectInfo.view();
    const MDagPath object_path = selectInfo.selectPath();
    const MMatrix object_matrix = object_path.inclusiveMatrix();

    SelectionRectangle rect(selectInfo);

    auto convert_world_to_screen = [&](MPoint point) -> openvdb::Vec2f {
        point *= object_matrix;
        short x_pos = 0; short y_pos = 0;
        view.worldToView(point, x_pos, y_pos);
        return { static_cast<float>(x_pos), static_cast<float>(y_pos) };
    };

    const std::array<MPoint, 8> world_points = {
            min,                          // 000
            MPoint(min.x, max.y, min.z),  // 010
            MPoint(min.x, max.y, max.z),  // 011
            MPoint(min.x, min.y, max.z),  // 001
            MPoint(max.x, min.y, min.z),  // 100
            MPoint(max.x, max.y, min.z),  // 110
            max,                          // 111
            MPoint(max.x, min.y, max.z)   // 101
    };

    typedef openvdb::Vec2f Point2D;
    Point2D points[8];

    for (int i = 0; i < 8; ++i)
        points[i] = convert_world_to_screen(world_points[i]);

    static const std::array<std::pair<int, int>, 12> line_array = {
            std::make_pair(0, 1), std::make_pair(1, 2), std::make_pair(2, 3), std::make_pair(3, 0),
            std::make_pair(4, 5), std::make_pair(5, 6), std::make_pair(6, 7), std::make_pair(7, 4),
            std::make_pair(0, 4), std::make_pair(1, 5), std::make_pair(2, 6), std::make_pair(3, 7)
    };

    for (auto line : line_array)
    {
        const auto& p0 = points[line.first];
        const auto& p1 = points[line.second];
        if (rect.clip_line(p0.x(), p0.y(), p1.x(), p1.y()))
        {
            MSelectionList item;
            item.add(object_path);
            selectInfo.addSelection(item, (world_points[line.first] + world_points[line.second]) * 0.5, selectionList, worldSpaceSelectPts, MSelectionMask::kSelectMeshes, false);
            return true;
        }
    }

    // No lines intersect the selection rectangle.
    // Check if the selection rectangle is within the projected bbox (enough to check a single point).
    auto vec2f_from_pair = [](const std::pair<float, float>& in) -> Point2D {
        return { in.first, in.second };
    };

    typedef std::vector<Point2D> Polygon2D;
    auto is_point_inside_polygon = [](const Point2D& test_point, const Polygon2D& polygon)
    {
        auto cross2d = [](const auto& v1, const auto& v2)
        {
            return v1.x() * v2.y() - v1.y() * v2.x();
        };

        const auto& last_point = polygon.back();
        const bool leftOfLastEdge = cross2d(polygon[0] - last_point, test_point - last_point) >= 0;
        for (int point_idx = 0; point_idx < polygon.size() - 1; ++point_idx)
        {
            const auto& v1 = polygon[point_idx];
            const auto& v2 = polygon[point_idx + 1];
            const bool leftOfEdge = cross2d(v2 - v1, test_point - v1) > 0;
            if (leftOfEdge != leftOfLastEdge)
                return false;
        }

        return true;
    };

    static const int index_array_for_faces[] = { 0, 1, 2, 3, 7, 6, 5, 4,
                                                 0, 3, 7, 4, 1, 5, 6, 2,
                                                 0, 4, 5, 1, 3, 2, 6, 7 };

    const auto test_point = openvdb::Vec2f(rect.getXMin(), rect.getYMin());
    for (int face_idx = 0; face_idx < 6; ++face_idx)
    {
        const auto index_base = index_array_for_faces + 4 * face_idx;
        auto polygon = Polygon2D(4);
        std::transform(index_base, index_base + 4, polygon.begin(), [&points](int point_index) {
            return points[point_index];
        });
        if (is_point_inside_polygon(test_point, polygon)) {
            MFloatVector face_midpoint;
            for (int i = 0; i < 4; ++i)
                face_midpoint += 0.25f * MFloatVector(world_points[index_base[i]]);

            MSelectionList item;
            item.add(object_path);
            selectInfo.addSelection(item, MPoint(face_midpoint), selectionList, worldSpaceSelectPts, MSelectionMask::kSelectMeshes, false);

            return true;
        }
    }

    return false;
}

bool VDBVisualizerShapeUI::canDrawUV() const
{
    return false;
}

namespace {
    const char *resample_cmd_node_flag_short = "n";
    const char *resample_cmd_node_flag_long = "node";
} // unnamed namespace

MSyntax VDBVisualizerUpdateMaxSliceCountCmd::create_syntax()
{
    MSyntax syntax;
    syntax.addFlag(resample_cmd_node_flag_short, resample_cmd_node_flag_long, MSyntax::kSelectionItem);
    return syntax;
}

MStatus VDBVisualizerUpdateMaxSliceCountCmd::doIt(const MArgList& args)
{
    MArgDatabase arg_data(syntax(), args);

    if (!arg_data.isFlagSet(resample_cmd_node_flag_short))
    {
        MGlobal::displayError("[openvdb] No node was passed to the command, use the -node(n) flag.");
        return MS::kFailure;
    }

    MStatus status;
    MSelectionList slist;
    status = arg_data.getFlagArgument(resample_cmd_node_flag_short, 0, slist);
    if (!status)
        return status;
    MObject node;
    status = slist.getDependNode(0, node);
    if (!status)
        return status;
    MFnDependencyNode dnode(node, &status);
    if (!status)
        return status;

    if (dnode.typeName() != VDBVisualizerShape::typeName)
    {
        MGlobal::displayError("[openvdb] Wrong node was passed to the command : " + dnode.name());
        return MS::kFailure;
    }
    auto vdb_visualizer_shape = dynamic_cast<VDBVisualizerShape*>(dnode.userNode());
    if (!vdb_visualizer_shape)
        return MS::kFailure;

    vdb_visualizer_shape->updateMaxSliceCount();

    return MS::kSuccess;
}

void VDBVisualizerShape::updateMaxSliceCount()
{
    m_vdb_data.max_slice_count = MPlug(thisMObject(), s_max_slice_count).asInt();
    // Trigger an update.
    MPlug(thisMObject(), s_apply_max_slice_count).setBool(true);
}