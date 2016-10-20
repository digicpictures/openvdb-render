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
    VDBSubSceneOverride(const MObject& obj);

    bool initVolumeRenderItem();
    bool initBBoxRenderItem();
    void updateShaderParams(const SliceShaderParams& shader_params);
    void updateDensityVolume(const DensityGridSpec& grid_spec);
    void updateBBoxGeometry(const openvdb::BBoxd& bbox);

    MObject m_object;

    // Rendering resources.

    struct ShaderInstanceDeleter {
        void operator()(MHWRender::MShaderInstance* ptr) const; 
    };
    typedef std::unique_ptr<MHWRender::MShaderInstance, ShaderInstanceDeleter> ShaderPtr;

    struct TextureDeleter {
        void operator()(MHWRender::MTexture* ptr) const; 
    };
    typedef std::unique_ptr<MHWRender::MTexture, TextureDeleter> TexturePtr;

    ShaderPtr m_volume_shader;
    TexturePtr m_volume_texture;

    struct RenderItem {
        MHWRender::MRenderItem* render_item;
        MHWRender::MVertexBufferArray vertex_buffer_array;
        std::unique_ptr<MHWRender::MVertexBuffer> position_buffer;
        std::unique_ptr<MHWRender::MIndexBuffer> index_buffer;
        MPxSubSceneOverride* parent;

        RenderItem(MPxSubSceneOverride* parent_) : render_item(nullptr), parent(parent_) {}
        void updateGeometry(const MBoundingBox& bbox);
    };

    RenderItem m_volume_render_item;
    RenderItem m_bbox_render_item;
};
