#include "vdb_subscene_override.h"

#include <cmath>
#include <chrono>
#include <sstream>

#include <boost/filesystem.hpp>

#include <maya/MFnDagNode.h>
#include <maya/MGlobal.h>
#include <maya/MSelectionList.h>
#include <maya/MShaderManager.h>

#include "paths.h"
#include "util.h"
#include "vdb_visualizer.h"
#include "volume_sampler.h"

// Disable "decorated name length exceeded, name was truncated" warning.
#pragma warning(disable: 4503)

#define MAX_SLICE_COUNT 64

using namespace MHWRender;

namespace {

    MFloatVector inline mayavecFromVec3f(const openvdb::Vec3f& vec)
    {
        return { vec.x(), vec.y(), vec.z() };
    }

    MBoundingBox inline mayabboxFromBBoxd(const openvdb::BBoxd& bbox)
    {
        return { mayavecFromVec3f(bbox.min()), mayavecFromVec3f(bbox.max()) };
    }

    const MShaderManager* getShaderManager()
    {
        auto renderer = MHWRender::MRenderer::theRenderer();
        if (!renderer) return nullptr;
        return renderer->getShaderManager();
    }

} // unnamed namespace

VDBSubSceneOverride::VDBSubSceneOverride(const MObject& obj)
    : MPxSubSceneOverride(obj),
    m_object(obj),
    m_volume_shader(nullptr),
    m_volume_render_item(this),
    m_bbox_render_item(this)
{
    const MShaderManager* shader_manager = getShaderManager();
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

bool VDBSubSceneOverride::initVolumeRenderItem()
{
    if (m_volume_render_item.render_item) {
        // Already initialized.
        return true;
    }

    if (!m_volume_shader) {
        // Couldn't load the shader.
        return false;
    }

    // Create render item.
    m_volume_render_item.render_item = MRenderItem::Create("vdb_volume", MRenderItem::RenderItemType::MaterialSceneItem, MHWRender::MGeometry::kTriangles);
    m_volume_render_item.render_item->setDrawMode(MGeometry::kAll);
    m_volume_render_item.render_item->castsShadows(false);
    m_volume_render_item.render_item->receivesShadows(false);
    if (!m_volume_render_item.render_item->setShader(m_volume_shader.get())) {
        LOG_ERROR("Could not set shader for volume render item.");
        return false;
    }

    // Build geometry.
    const unsigned int slice_count = MAX_SLICE_COUNT;

    // - Vertices
    // Note: descriptor name (first ctor arg) MUST be "", or setGeometryForRenderItem will return kFailure.
    const MVertexBufferDescriptor pos_desc("", MGeometry::kPosition, MGeometry::kFloat, 3);
    m_volume_render_item.position_buffer.reset(new MVertexBuffer(pos_desc));
    const auto vertex_count = slice_count * 4;
    MFloatVector* positions = reinterpret_cast<MFloatVector*>(m_volume_render_item.position_buffer->acquire(static_cast<unsigned int>(vertex_count), false));
    for (unsigned int i = 0; i < slice_count; ++i) {
        const float z = i * 1.0f / (slice_count - 1.0f);
        positions[4*i+0] = MFloatVector(0.0, 0.0, z);
        positions[4*i+1] = MFloatVector(1.0, 0.0, z);
        positions[4*i+2] = MFloatVector(0.0, 1.0, z);
        positions[4*i+3] = MFloatVector(1.0, 1.0, z);
    }
    m_volume_render_item.position_buffer->commit(positions);

    // - Indices
    m_volume_render_item.index_buffer.reset(new MIndexBuffer(MGeometry::kUnsignedInt32));
    const auto index_count = slice_count * 6;
    unsigned int* indices = reinterpret_cast<unsigned int*>(m_volume_render_item.index_buffer->acquire(index_count, false));
    for (unsigned int i = 0; i < slice_count; ++i) {
        indices[6*i+0] = 4*i+0;
        indices[6*i+1] = 4*i+1;
        indices[6*i+2] = 4*i+3;
        indices[6*i+3] = 4*i+0;
        indices[6*i+4] = 4*i+3;
        indices[6*i+5] = 4*i+2;
    }
    m_volume_render_item.index_buffer->commit(indices);

    m_volume_render_item.vertex_buffer_array.clear();
    CHECK_MSTATUS(m_volume_render_item.vertex_buffer_array.addBuffer("pos_model", m_volume_render_item.position_buffer.get()));
    return true;
}

bool VDBSubSceneOverride::initBBoxRenderItem()
{
    if (m_bbox_render_item.render_item) {
        // Already initialized.
        return true;
    }

    const MShaderManager* shader_manager = getShaderManager();
    if (!shader_manager) {
        return false;
    }

    static auto shader = shader_manager->getStockShader(MShaderManager::k3dSolidShader);
    if (!shader) {
        LOG_ERROR("Couldn't get stock shader: k3dSolidShader.");
        return false;
    }
    static const float color[] = {67.0f / 255.0f, 1.0f, 163.0f / 255.0f, 1.0f};
    CHECK_MSTATUS(shader->setParameter("solidColor", color));

    m_bbox_render_item.render_item = MRenderItem::Create("bounding_box", MRenderItem::RenderItemType::DecorationItem, MHWRender::MGeometry::kLines);
    if (!m_bbox_render_item.render_item) {
        LOG_ERROR("Failed to create bbox render item.");
        return false;
    }
    m_bbox_render_item.render_item->setDrawMode(MGeometry::kAll);
    if (!m_bbox_render_item.render_item->setShader(shader)) {
        LOG_ERROR("Failed to set shader for bbox render item.");
        return false;
    }

    // Build geometry.
    // Note: descriptor name (first ctor arg) MUST be "", or setGeometryForRenderItem will return kFailure.
    const MVertexBufferDescriptor pos_desc("", MGeometry::kPosition, MGeometry::kFloat, 3);
    m_bbox_render_item.position_buffer.reset(new MVertexBuffer(pos_desc));

    m_bbox_render_item.index_buffer.reset(new MIndexBuffer(MGeometry::kUnsignedInt32));
    constexpr auto index_count = 2 * 12;
    static const uint32_t BOX_INDICES[] = { 0, 1, 1, 3, 3, 2, 2, 0, 4, 5, 5, 7, 7, 6, 6, 4, 0, 4, 1, 5, 3, 7, 2, 6 };
    CHECK_MSTATUS(m_bbox_render_item.index_buffer->update(BOX_INDICES, 0, index_count, true));

    return true;
}

void VDBSubSceneOverride::updateBBoxGeometry(const openvdb::BBoxd& bbox)
{
    constexpr auto vertex_count = 8;
    static MFloatVector positions[vertex_count];
    for (unsigned int i = 0; i < vertex_count; ++i) {
        positions[i] = mayavecFromVec3f(openvdb::Vec3d(i & 1, (i & 2) >> 1, (i & 4) >> 2) * bbox.extents() + bbox.min());
    }
    CHECK_MSTATUS(m_bbox_render_item.position_buffer->update(positions, 0, vertex_count, true));

    m_bbox_render_item.vertex_buffer_array.clear();
    CHECK_MSTATUS(m_bbox_render_item.vertex_buffer_array.addBuffer("pos_model", m_bbox_render_item.position_buffer.get()));

    m_bbox_render_item.updateGeometry(mayabboxFromBBoxd(bbox));
}

void VDBSubSceneOverride::updateShaderParams(const SliceShaderParams& params)
{
    // Shading parameters.
    m_volume_shader->setParameter("light_dir", params.light_direction);
    m_volume_shader->setParameter("light_color", params.light_color);
    m_volume_shader->setParameter("scattering", params.scattering);
    m_volume_shader->setParameter("absorption", params.absorption);
    m_volume_shader->setParameter("shadow_gain", params.shadow_gain);
    m_volume_shader->setParameter("shadow_sample_count", params.shadow_sample_count);
    m_volume_shader->setParameter("slice_size_model", params.slice_size);
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
        m_volume_shader->setParameter("volume_texture", volume_texture_resource);
        return;
    }

    const auto bbox_is = getIndexSpaceBoundingBox(grid_ptr.get());
    const auto bbox_ws = grid_ptr->transform().indexToWorld(bbox_is);
    CHECK_MSTATUS(m_volume_shader->setParameter("volume_size_model", mayavecFromVec3f(bbox_ws.extents())));
    CHECK_MSTATUS(m_volume_shader->setParameter("volume_origin_model", mayavecFromVec3f(bbox_ws.min())));

    // Set bbox for render items.
    m_volume_render_item.updateGeometry(mayabboxFromBBoxd(bbox_ws));
    updateBBoxGeometry(bbox_ws);

    // Sample the multi resolution grid at regular intervals.
    VolumeSampler volume_sampler;
    const auto texture_extents = openvdb::Coord(MAX_SLICE_COUNT, MAX_SLICE_COUNT, MAX_SLICE_COUNT);
    auto volume = volume_sampler.sampleGridWithMipmapFilter(*grid_ptr, texture_extents);
    m_volume_texture.reset(volume.texture);

    // Set shader params.
    m_volume_shader->setParameter("min_voxel_value", volume.value_range.min);
    m_volume_shader->setParameter("max_voxel_value", volume.value_range.max);

    MHWRender::MTextureAssignment volume_texture_resource;
    volume_texture_resource.texture = m_volume_texture.get();
    m_volume_shader->setParameter("volume_texture", volume_texture_resource);

    MHWRender::MSamplerStateDesc volume_sampler_state_desc;
    volume_sampler_state_desc.filter = MHWRender::MSamplerState::kMinMagMipLinear;
    volume_sampler_state_desc.addressU = MHWRender::MSamplerState::kTexClamp;
    volume_sampler_state_desc.addressV = MHWRender::MSamplerState::kTexClamp;
    volume_sampler_state_desc.addressW = MHWRender::MSamplerState::kTexClamp;
    const MHWRender::MSamplerState *volume_sampler_resource = MStateManager::acquireSamplerState(volume_sampler_state_desc);
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
    if (!m_volume_render_item.render_item) {
        if (!initVolumeRenderItem()) {
            return;
        }
        if (!container.add(m_volume_render_item.render_item)) {
            LOG_ERROR("Could not add volume render item.");
            return;
        }
    }

    if (!m_bbox_render_item.render_item) {
        if (!initBBoxRenderItem()) {
            return;
        }

        if (!container.add(m_bbox_render_item.render_item)) {
            LOG_ERROR("Could not add bbox render item.");
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
    m_bbox_render_item.render_item->enable(volume_is_selected);

    // Set world matrix.
    m_volume_render_item.render_item->setMatrix(&path.inclusiveMatrix());
    m_bbox_render_item.render_item->setMatrix(&path.inclusiveMatrix());

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

MString VDBSubSceneOverride::registrantId = "VDBVisualizerSubSceneOverride";

void VDBSubSceneOverride::ShaderInstanceDeleter::operator()(MShaderInstance* ptr) const
{
    if (ptr) {
        MHWRender::MRenderer::theRenderer()->getShaderManager()->releaseShader(ptr);
    }
}

void VDBSubSceneOverride::TextureDeleter::operator()(MTexture* ptr) const
{
    if (ptr) {
        MHWRender::MRenderer::theRenderer()->getTextureManager()->releaseTexture(ptr);
    }
}

void VDBSubSceneOverride::RenderItem::updateGeometry(const MBoundingBox& bbox)
{
    // Note: render item has to be added to the MSubSceneContainer before calling setGeometryForRenderItem.
    CHECK_MSTATUS(parent->setGeometryForRenderItem(*render_item, vertex_buffer_array, *index_buffer, &bbox));
}
