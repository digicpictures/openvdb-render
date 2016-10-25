#include "vdb_subscene_override.h"

#include <cmath>
#include <chrono>
#include <sstream>

#include <maya/MFnDagNode.h>
#include <maya/MHWGeometryUtilities.h>
#include <maya/MGlobal.h>
#include <maya/MSelectionList.h>
#include <maya/MShaderManager.h>

#include "paths.h"
#include "util.h"
#include "vdb_visualizer.h"
#include "volume_sampler.h"

#define MAX_SLICE_COUNT 64

namespace {

    MFloatVector inline mayavecFromVec3f(const openvdb::Vec3f& vec)
    {
        return { vec.x(), vec.y(), vec.z() };
    }

    MBoundingBox inline mayabboxFromBBoxd(const openvdb::BBoxd& bbox)
    {
        return { mayavecFromVec3f(bbox.min()), mayavecFromVec3f(bbox.max()) };
    }

    const MHWRender::MShaderManager* getShaderManager()
    {
        auto renderer = MHWRender::MRenderer::theRenderer();
        if (!renderer) return nullptr;
        return renderer->getShaderManager();
    }

} // unnamed namespace

VDBSubSceneOverride::VDBSubSceneOverride(const MObject& obj)
    : MPxSubSceneOverride(obj),
    m_object(obj)
{
    const MHWRender::MShaderManager* shader_manager = getShaderManager();
    assert(shader_manager);

    // Load volume shader from effect file.
    auto effect_file_path = Paths::getVolumeEffectFile();
    m_volume_shader.reset(shader_manager->getEffectsFileShader(effect_file_path.c_str(), "Main", 0, 0, false));
    if (!m_volume_shader) {
        std::stringstream ss;
        ss << "Cannot load cgfx file: " << effect_file_path;
        LOG_ERROR(ss.str());
        return;
    }
    m_volume_shader->setIsTransparent(true);
    m_volume_shader->setParameter("max_slice_count", MAX_SLICE_COUNT);
}

VDBSubSceneOverride::~VDBSubSceneOverride()
{
}

bool VDBSubSceneOverride::requiresUpdate(const MHWRender::MSubSceneContainer& container, const MHWRender::MFrameContext& frameContext) const
{
    return true;
}

bool VDBSubSceneOverride::initSliceRenderables(MHWRender::MSubSceneContainer& container)
{
    if (m_slices_renderable.render_item != nullptr) {
        // Already initialized.
        return true;
    }

    if (!m_volume_shader) {
        // Couldn't load the shader.
        return false;
    }

    m_slices_renderable.render_item = MHWRender::MRenderItem::Create("vdb_volume_slices", MHWRender::MRenderItem::RenderItemType::MaterialSceneItem, MHWRender::MGeometry::kTriangles);
    m_slices_renderable.render_item->setDrawMode(MHWRender::MGeometry::kAll);
    m_slices_renderable.render_item->castsShadows(false);
    m_slices_renderable.render_item->receivesShadows(false);
    if (!m_slices_renderable.render_item->setShader(m_volume_shader.get())) {
        LOG_ERROR("Could not set shader for volume render item.");
        return false;
    }

    // Build geometry.
    const unsigned int slice_count = MAX_SLICE_COUNT;

    // - Vertices
    // Note: descriptor name (first ctor arg) MUST be "", or setGeometryForRenderItem will return kFailure.
    const MHWRender::MVertexBufferDescriptor pos_desc("", MHWRender::MGeometry::kPosition, MHWRender::MGeometry::kFloat, 3);
    m_slices_renderable.position_buffer.reset(new MHWRender::MVertexBuffer(pos_desc));
    const auto vertex_count = slice_count * 4;
    MFloatVector* positions = reinterpret_cast<MFloatVector*>(m_slices_renderable.position_buffer->acquire(vertex_count, true));
    for (unsigned int i = 0; i < slice_count; ++i) {
        const float z = i * 1.0f / (slice_count - 1.0f);
        positions[4*i+0] = MFloatVector(0.0, 0.0, z);
        positions[4*i+1] = MFloatVector(1.0, 0.0, z);
        positions[4*i+2] = MFloatVector(0.0, 1.0, z);
        positions[4*i+3] = MFloatVector(1.0, 1.0, z);
    }
    m_slices_renderable.position_buffer->commit(positions);

    // - Indices
    m_slices_renderable.index_buffer.reset(new MHWRender::MIndexBuffer(MHWRender::MGeometry::kUnsignedInt32));
    const auto index_count = slice_count * 6;
    unsigned int* indices = reinterpret_cast<unsigned int*>(m_slices_renderable.index_buffer->acquire(index_count, true));
    for (unsigned int i = 0; i < slice_count; ++i) {
        indices[6*i+0] = 4*i+0;
        indices[6*i+1] = 4*i+1;
        indices[6*i+2] = 4*i+3;
        indices[6*i+3] = 4*i+0;
        indices[6*i+4] = 4*i+3;
        indices[6*i+5] = 4*i+2;
    }
    m_slices_renderable.index_buffer->commit(indices);

    m_slices_renderable.vertex_buffer_array.clear();
    CHECK_MSTATUS(m_slices_renderable.vertex_buffer_array.addBuffer("pos_model", m_slices_renderable.position_buffer.get()));

    if (!container.add(m_slices_renderable.render_item)) {
        LOG_ERROR("Could not add m_slices_renderable render item.");
        return false;
    }

    return true;
}

bool VDBSubSceneOverride::initBBoxRenderable(MHWRender::MSubSceneContainer& container)
{
    if (m_bbox_renderable.render_item) {
        // Already initialized.
        return true;
    }

    const MHWRender::MShaderManager* shader_manager = getShaderManager();
    if (!shader_manager) {
        return false;
    }

    static auto shader = shader_manager->getStockShader(MHWRender::MShaderManager::k3dSolidShader);
    if (!shader) {
        LOG_ERROR("Couldn't get stock shader: k3dSolidShader.");
        return false;
    }

    m_bbox_renderable.render_item = MHWRender::MRenderItem::Create("bounding_box", MHWRender::MRenderItem::RenderItemType::DecorationItem, MHWRender::MGeometry::kLines);
    if (!m_bbox_renderable.render_item) {
        LOG_ERROR("Failed to create bbox render item.");
        return false;
    }

    if (!m_bbox_renderable.render_item->setShader(shader)) {
        LOG_ERROR("Failed to set shader for bbox render item.");
        return false;
    }

    m_bbox_renderable.render_item->setDrawMode(MHWRender::MGeometry::kAll);
    m_bbox_renderable.render_item->depthPriority(MHWRender::MRenderItem::sActiveWireDepthPriority);

    // Note: descriptor name (first ctor arg) MUST be "", or setGeometryForRenderItem will return kFailure.
    const MHWRender::MVertexBufferDescriptor pos_desc("", MHWRender::MGeometry::kPosition, MHWRender::MGeometry::kFloat, 3);
    m_bbox_renderable.position_buffer.reset(new MHWRender::MVertexBuffer(pos_desc));
    m_bbox_renderable.vertex_buffer_array.clear();
    CHECK_MSTATUS(m_bbox_renderable.vertex_buffer_array.addBuffer("pos_model", m_bbox_renderable.position_buffer.get()));

    m_bbox_renderable.index_buffer.reset(new MHWRender::MIndexBuffer(MHWRender::MGeometry::kUnsignedInt32));
    constexpr auto index_count = 2 * 12;
    static const uint32_t BOX_WIREFRAME_INDICES[] = { 0, 1, 1, 3, 3, 2, 2, 0, 4, 5, 5, 7, 7, 6, 6, 4, 0, 4, 1, 5, 3, 7, 2, 6 };
    CHECK_MSTATUS(m_bbox_renderable.index_buffer->update(BOX_WIREFRAME_INDICES, 0, index_count, true));

    if (!container.add(m_bbox_renderable.render_item)) {
        LOG_ERROR("Could not add bbox render item.");
        return false;
    }

    return true;
}

void VDBSubSceneOverride::updateBBoxGeometry(const openvdb::BBoxd& bbox)
{
    constexpr auto vertex_count = 8;
    MFloatVector* positions = static_cast<MFloatVector*>(m_bbox_renderable.position_buffer->acquire(vertex_count, true));
    for (unsigned int i = 0; i < vertex_count; ++i) {
        positions[i] = mayavecFromVec3f(openvdb::Vec3d(i & 1, (i & 2) >> 1, (i & 4) >> 2) * bbox.extents() + bbox.min());
    }
    m_bbox_renderable.position_buffer->commit(positions);

    updateGeometry(m_bbox_renderable, mayabboxFromBBoxd(bbox));
}

void VDBSubSceneOverride::updateShaderParams(const SliceShaderParams& params)
{
    // Shading parameters.
    CHECK_MSTATUS(m_volume_shader->setParameter("light_dir", params.light_direction));
    CHECK_MSTATUS(m_volume_shader->setParameter("light_color", params.light_color));
    CHECK_MSTATUS(m_volume_shader->setParameter("scattering", params.scattering));
    CHECK_MSTATUS(m_volume_shader->setParameter("absorption", params.absorption));
    CHECK_MSTATUS(m_volume_shader->setParameter("shadow_gain", params.shadow_gain));
    CHECK_MSTATUS(m_volume_shader->setParameter("shadow_sample_count", params.shadow_sample_count));
    CHECK_MSTATUS(m_volume_shader->setParameter("slice_size_model", params.slice_size));
}

namespace {

    openvdb::FloatGrid::ConstPtr loadDensityGrid(const DensityGridSpec& grid_spec)
    {
        if (!grid_spec.vdb_file || !grid_spec.vdb_file->isOpen()) {
            return nullptr;
        }

        openvdb::GridBase::ConstPtr grid_base_ptr;
        try {
            grid_base_ptr = grid_spec.vdb_file->readGrid(grid_spec.grid_name);
        }
        catch (const openvdb::Exception& e) {
            std::stringstream ss;
            ss << "Error reading grid " << grid_spec.grid_name << ": " << e.what();
            LOG_ERROR(ss.str());
            return nullptr;
        }

        auto grid_ptr = openvdb::gridConstPtrCast<openvdb::FloatGrid>(grid_base_ptr);
        if (!grid_ptr) {
            LOG_ERROR("Grid is not a FloatGrid.");
            return nullptr;
        }

        return grid_ptr;
    }

}

void VDBSubSceneOverride::updateDensityVolume(const DensityGridSpec& grid_spec)
{
    auto grid_ptr = loadDensityGrid(grid_spec);
    if (!grid_ptr) {
        m_volume_texture.reset();
        MHWRender::MTextureAssignment volume_texture_resource;
        volume_texture_resource.texture = nullptr;
        CHECK_MSTATUS(m_volume_shader->setParameter("volume_texture", volume_texture_resource));
        return;
    }

    const auto bbox_is = getIndexSpaceBoundingBox(grid_ptr.get());
    const auto bbox_ws = grid_ptr->transform().indexToWorld(bbox_is);
    CHECK_MSTATUS(m_volume_shader->setParameter("volume_size_model", mayavecFromVec3f(bbox_ws.extents())));
    CHECK_MSTATUS(m_volume_shader->setParameter("volume_origin_model", mayavecFromVec3f(bbox_ws.min())));

    // Set bbox for render items.
    updateGeometry(m_slices_renderable, mayabboxFromBBoxd(bbox_ws));
    updateBBoxGeometry(bbox_ws);

    // Sample the multi resolution grid at regular intervals.
    VolumeSampler volume_sampler;
    const auto texture_extents = openvdb::Coord(MAX_SLICE_COUNT, MAX_SLICE_COUNT, MAX_SLICE_COUNT);
    auto volume = volume_sampler.sampleGridWithMipmapFilter(*grid_ptr, texture_extents);
    m_volume_texture.reset(volume.texture);

    // Set shader params.
    CHECK_MSTATUS(m_volume_shader->setParameter("min_voxel_value", volume.value_range.min));
    CHECK_MSTATUS(m_volume_shader->setParameter("max_voxel_value", volume.value_range.max));

    MHWRender::MTextureAssignment volume_texture_resource;
    volume_texture_resource.texture = m_volume_texture.get();
    m_volume_shader->setParameter("volume_texture", volume_texture_resource);

    MHWRender::MSamplerStateDesc volume_sampler_state_desc;
    volume_sampler_state_desc.filter = MHWRender::MSamplerState::kMinMagMipLinear;
    volume_sampler_state_desc.addressU = MHWRender::MSamplerState::kTexClamp;
    volume_sampler_state_desc.addressV = MHWRender::MSamplerState::kTexClamp;
    volume_sampler_state_desc.addressW = MHWRender::MSamplerState::kTexClamp;
    const MHWRender::MSamplerState *volume_sampler_resource = MHWRender::MStateManager::acquireSamplerState(volume_sampler_state_desc);
    m_volume_shader->setParameter("volume_sampler", *volume_sampler_resource);
}

namespace {

    bool isPathSelected(MDagPath path)
    {
        MSelectionList selectedList;
        MGlobal::getActiveSelectionList(selectedList);
        do {
            if (selectedList.hasItem(path)) {
                return true;
            }
        } while (path.pop());
        return false;
    }
}

void VDBSubSceneOverride::update(MHWRender::MSubSceneContainer& container, const MHWRender::MFrameContext& frameContext)
{
    if (m_slices_renderable.render_item == nullptr) {
        if (!initSliceRenderables(container)) {
            return;
        }
    }

    if (m_bbox_renderable.render_item == nullptr) {
        if (!initBBoxRenderable(container)) {
            return;
        }
    }

    MStatus status;
    MFnDagNode node(m_object, &status);
    if (!status) return;

    MDagPath path;
    node.getPath(path);

    // Handle selection.
    bool volume_is_selected = isPathSelected(path);
    m_bbox_renderable.render_item->enable(volume_is_selected);

    // Set wireframe color.
    const auto color = MHWRender::MGeometryUtilities::wireframeColor(path);
    const float color_as_array[] = { color.r, color.g, color.b, color.a };
    CHECK_MSTATUS(m_bbox_renderable.render_item->getShader()->setParameter("solidColor", color_as_array));

    // Set world matrix.
    m_slices_renderable.render_item->setMatrix(&path.inclusiveMatrix());
    m_bbox_renderable.render_item->setMatrix(&path.inclusiveMatrix());

    // Get VDB shape object from associated node.
    typedef VDBVisualizerShape ShapeType;
    ShapeType* shape_node = dynamic_cast<ShapeType*>(node.userNode());
    if (!shape_node) return;

    auto data = shape_node->get_update();
    if (!data) {
        return;
    }
    if (data->density_grid_data.isDirty()) {
        updateDensityVolume(data->density_grid_data.clearAndGet());
    }
    if (data->slice_shader_params.isDirty()) {
        updateShaderParams(data->slice_shader_params.clearAndGet());
    }
}

MString VDBSubSceneOverride::s_registrant_id = "VDBVisualizerSubSceneOverride";

void VDBSubSceneOverride::updateGeometry(const Renderable& renderable, const MBoundingBox& bbox)
{
    // Note: render item has to be added to the MSubSceneContainer before calling setGeometryForRenderItem.
    CHECK_MSTATUS(setGeometryForRenderItem(*renderable.render_item, renderable.vertex_buffer_array, *renderable.index_buffer, &bbox));
}
