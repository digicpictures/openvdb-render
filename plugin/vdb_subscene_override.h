#pragma once

#include <array>
#include <memory>

#include <maya/MPxSubSceneOverride.h>
#include <maya/MPxCommand.h>
#include <maya/MSyntax.h>

#include "vdb_visualizer.h"
#include "vdb_subscene_utils.hpp"

struct VDBSubSceneOverrideData;
class SlicedDisplay;

class VDBSubSceneOverride : public MHWRender::MPxSubSceneOverride {
public:
    static MPxSubSceneOverride* creator(const MObject& obj);

    VDBSubSceneOverride(const MObject& obj);
    virtual ~VDBSubSceneOverride();

    virtual MHWRender::DrawAPI supportedDrawAPIs() const;
    virtual void update(MHWRender::MSubSceneContainer& container, const MHWRender::MFrameContext& frameContext);
    virtual bool requiresUpdate(const MHWRender::MSubSceneContainer& container, const MHWRender::MFrameContext& frameContext) const;

    static MString registrantId;
private:

    MObject m_object;
    VDBVisualizerShape* p_vdb_visualizer;
    std::unique_ptr<VDBSubSceneOverrideData> p_data;

    std::unique_ptr<MHWRender::MVertexBuffer> p_bbox_position;
    std::unique_ptr<MHWRender::MIndexBuffer> p_bbox_indices;

    std::unique_ptr<MHWRender::MVertexBuffer> p_position_buffer;
    std::unique_ptr<MHWRender::MVertexBuffer> p_color_buffer;
    std::unique_ptr<MHWRender::MIndexBuffer> p_index_buffer;

    struct shader_instance_deleter {
        void operator()(MHWRender::MShaderInstance* p);
    };

    std::unique_ptr<MHWRender::MShaderInstance, shader_instance_deleter> p_point_cloud_shader;
    // max is not constexpr
    static constexpr size_t sampler_mem_size = sizeof(FloatToRGBSampler) > sizeof(Vec3SToRGBSampler)
                                               ? sizeof(FloatToRGBSampler) : sizeof(Vec3SToRGBSampler);

    typedef std::array<char, sampler_mem_size> sampler_mem_area;
    sampler_mem_area m_scattering_sampler;
    sampler_mem_area m_emission_sampler;
    sampler_mem_area m_attenuation_sampler;

    std::unique_ptr<SlicedDisplay> m_sliced_display;
};

class VDBVolumeCacheMemonyLimitCmd : public MPxCommand
{
public:
    VDBVolumeCacheMemonyLimitCmd() {}
    ~VDBVolumeCacheMemonyLimitCmd() {}

    static void* creator() { return new VDBVolumeCacheMemonyLimitCmd(); }
    static MSyntax create_syntax();

    MStatus doIt(const MArgList& args);
};
