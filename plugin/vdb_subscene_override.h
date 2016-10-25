#pragma once

#include <memory>
#include <boost/container/flat_map.hpp>
#include <maya/MPxSubSceneOverride.h>
#include "util.h"
#include "vdb_visualizer_data.h"

class VDBSubSceneOverride : public MHWRender::MPxSubSceneOverride
{
public:
    static MString s_registrant_id;
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

    bool initSliceRenderables(MHWRender::MSubSceneContainer& container);
    bool initBBoxRenderable(MHWRender::MSubSceneContainer& container);
    void updateShaderParams(const SliceShaderParams& shader_params);
    void updateDensityVolume(const DensityGridSpec& grid_spec);
    void updateBBoxGeometry(const openvdb::BBoxd& bbox);

    MObject m_object;

    // Rendering resources.
    ShaderPtr m_volume_shader;
    TexturePtr m_volume_texture;

    struct Renderable {
        MHWRender::MRenderItem* render_item;
        MHWRender::MVertexBufferArray vertex_buffer_array;
        std::unique_ptr<MHWRender::MVertexBuffer> position_buffer;
        std::unique_ptr<MHWRender::MIndexBuffer> index_buffer;

        Renderable() : render_item(nullptr) {}
    };

    void updateGeometry(const Renderable& renderable, const MBoundingBox& bbox);

    Renderable m_slices_renderable;
    Renderable m_bbox_renderable;
};
