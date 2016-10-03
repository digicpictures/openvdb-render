#pragma once

#include <memory>
#include <boost/container/flat_map.hpp>
#include <maya/MPxSubSceneOverride.h>
#include "vdb_visualizer_data.h"

class VDBSubSceneOverride : public MHWRender::MPxSubSceneOverride
{
public:
    static MString registrantId;
    static MPxSubSceneOverride* creator(const MObject& obj)
    {
        return new VDBSubSceneOverride(obj);
    }

    virtual ~VDBSubSceneOverride() override;

    virtual MHWRender::DrawAPI supportedDrawAPIs() const override
    {
        return MHWRender::kOpenGL;
    }

    virtual bool requiresUpdate(
        const MHWRender::MSubSceneContainer& container,
        const MHWRender::MFrameContext& frameContext) const override;

    virtual void update(
        MHWRender::MSubSceneContainer& container,
        const MHWRender::MFrameContext& frameContext) override;

private:
    VDBSubSceneOverride(const MObject& obj)
        : MPxSubSceneOverride(obj),
        m_object(obj),
        m_num_instances(0),
        m_volume_render_item(nullptr),
        m_volume_shader(nullptr),
        m_volume_position_buffer(nullptr),
        m_volume_index_buffer(nullptr) {}

    void initRenderItem();
    void updateShaderParams(const VDBVisualizerData* data);
    void updateGeometry(const VDBVisualizerData* data);

    MObject m_object;

    struct InstanceInfo
    {
        MMatrix m_transform;
        bool m_is_selected;

        InstanceInfo() : m_is_selected(false) {}
        InstanceInfo(const MMatrix& matrix, bool selected) : m_transform(matrix), m_is_selected(selected) {}
        bool operator==(const InstanceInfo& rhs) const { return m_is_selected == rhs.m_is_selected && m_transform.isEquivalent(rhs.m_transform); }
        bool operator!=(const InstanceInfo& rhs) const { return !(*this == rhs); }
    };
    //typedef std::map<unsigned int, InstanceInfo> InstanceInfoMap;
    typedef boost::container::flat_map<unsigned int, InstanceInfo> InstanceInfoMap;
    InstanceInfoMap m_instance_info_cache;
    int m_num_instances;

    struct ShaderInstanceDeleter {
        void operator()(MHWRender::MShaderInstance* ptr) const; 
    };
    typedef std::unique_ptr<MHWRender::MShaderInstance, ShaderInstanceDeleter> ShaderPtr;

    struct TextureDeleter {
        void operator()(MHWRender::MTexture* ptr) const; 
    };
    typedef std::unique_ptr<MHWRender::MTexture, TextureDeleter> TexturePtr;

    // Rendering resources.
    ShaderPtr m_volume_shader;
    TexturePtr m_volume_texture;
    MHWRender::MRenderItem* m_volume_render_item;
    std::unique_ptr<MHWRender::MVertexBuffer> m_volume_position_buffer;
    std::unique_ptr<MHWRender::MIndexBuffer> m_volume_index_buffer;
};
