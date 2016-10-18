#include "vdb_subscene_override.h"

#include <cmath>
#include <chrono>

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
        return{ vec.x(), vec.y(), vec.z() };
    }

    MBoundingBox inline mayabboxFromBBoxd(const openvdb::BBoxd& bbox)
    {
        return { mayavecFromVec3f(bbox.min()), mayavecFromVec3f(bbox.max()) };
    }

} // unnamed namespace

struct VDBSubSceneOverrideError
{
    VDBSubSceneOverrideError(const char *what) : m_what(what) {}
    const char *what() const { return m_what; }
    const char *m_what;
};

VDBSubSceneOverride::~VDBSubSceneOverride()
{
}

bool VDBSubSceneOverride::requiresUpdate(const MHWRender::MSubSceneContainer& container, const MHWRender::MFrameContext& frameContext) const
{
    return true;
}

bool VDBSubSceneOverride::initRenderItem()
{
    assert(!m_volume_render_item);

    // Obtain pointers to rendering managers.
    MRenderer* renderer = MRenderer::theRenderer();
    if (!renderer) {
        return false;
    }
    const MShaderManager* shader_manager = renderer->getShaderManager();
    if (!shader_manager) {
        std::cerr << "Unable to obtain pointer to shader manager." << std::endl;
        return false;
    }

    if (!m_volume_shader) {
        // Create shader instance.
        auto effect_file_path = Paths::getVolumeEffectFile();
        m_volume_shader.reset(shader_manager->getEffectsFileShader(effect_file_path.c_str(), "Main", 0, 0, false));
        if (!m_volume_shader) {
            static bool first_time = true;
            if (first_time) {
                std::cerr << "Cannot load cgfx file: " << effect_file_path << std::endl;
                first_time = false;
            }
            return false;
        }
        m_volume_shader->setIsTransparent(true);
        m_volume_shader->setParameter("max_slice_count", MAX_SLICE_COUNT);
    }

    // Create render item.
    m_volume_render_item = MRenderItem::Create("vdb_volume", MRenderItem::RenderItemType::MaterialSceneItem, MHWRender::MGeometry::kTriangles);
    m_volume_render_item->setDrawMode(MGeometry::kAll);
    m_volume_render_item->castsShadows(false);
    m_volume_render_item->receivesShadows(false);
    if (!m_volume_render_item->setShader(m_volume_shader.get())) {
        std::cerr << "VDBSubSceneOverride::initRenderItem: Could not set volume shader." << std::endl;
        return false;
    }

    return true;
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

void VDBSubSceneOverride::updateDensityVolume(const DensityGridData& grid_data)
{
    // Clean-up routine.
    const auto clean_texture = [&m_volume_texture = m_volume_texture, &m_volume_shader = m_volume_shader]() {
        m_volume_texture.reset();
        MHWRender::MTextureAssignment volume_texture_resource;
        volume_texture_resource.texture = nullptr;
        m_volume_shader->setParameter("volume_texture", volume_texture_resource);
    };

    // Bail unless vdb_file is specified and loaded.
    if (!grid_data.vdb_file || !grid_data.vdb_file->isOpen()) {
        clean_texture();
        return;
    }

    // Load the density grid.
    openvdb::FloatGrid::ConstPtr grid_ptr;
    try {
        grid_ptr = openvdb::gridConstPtrCast<openvdb::FloatGrid>(grid_data.vdb_file->readGrid(grid_data.grid_name));
    } catch (const openvdb::Exception& e) {
        std::cerr << "error reading grid " << grid_data.grid_name << ": " << e.what() << std::endl;
        clean_texture();
        return;
    }

    // Bail unless we have a valid grid pointer.
    if (!grid_ptr) {
        clean_texture();
        return;
    }

    // Set model space volume bounds (world space in the vdb grid's perspective).
    const auto bbox_is = getIndexSpaceBoundingBox(grid_ptr.get());
    const auto grid_bbox_ws = grid_ptr->transform().indexToWorld(bbox_is);
    CHECK_MSTATUS(m_volume_shader->setParameter("volume_size_model", mayavecFromVec3f(grid_bbox_ws.extents())));
    CHECK_MSTATUS(m_volume_shader->setParameter("volume_origin_model", mayavecFromVec3f(grid_bbox_ws.min())));

    // Set render item bbox.
    const auto bbox_ws_maya = mayabboxFromBBoxd(grid_ptr->transform().indexToWorld(bbox_is));
    CHECK_MSTATUS(setGeometryForRenderItem(*m_volume_render_item, m_volume_vertex_buffers, *m_volume_index_buffer, &bbox_ws_maya));

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

void VDBSubSceneOverride::updateGeometry(unsigned int slice_count)
{
    // - Vertices
    // Note: descriptor name (first ctor arg) MUST be "", or setGeometryForRenderItem will return kFailure.
    const MVertexBufferDescriptor pos_desc("", MGeometry::kPosition, MGeometry::kFloat, 3);
    m_volume_position_buffer.reset(new MVertexBuffer(pos_desc));
    const auto vertex_count = slice_count * 4;
    MFloatVector* positions = reinterpret_cast<MFloatVector*>(m_volume_position_buffer->acquire(static_cast<unsigned int>(vertex_count), false));
    for (unsigned int i = 0; i < slice_count; ++i) {
        const float z = i * 1.0f / (slice_count - 1.0f);
        positions[4*i+0] = MFloatVector(0.0, 0.0, z);
        positions[4*i+1] = MFloatVector(1.0, 0.0, z);
        positions[4*i+2] = MFloatVector(0.0, 1.0, z);
        positions[4*i+3] = MFloatVector(1.0, 1.0, z);
    }
    m_volume_position_buffer->commit(positions);

    // - Indices
    m_volume_index_buffer.reset(new MIndexBuffer(MGeometry::kUnsignedInt32));
    const auto index_count = slice_count * 6;
    unsigned int* indices = reinterpret_cast<unsigned int*>(m_volume_index_buffer->acquire(index_count, false));
    for (unsigned int i = 0; i < slice_count; ++i) {
        indices[6*i+0] = 4*i+0;
        indices[6*i+1] = 4*i+1;
        indices[6*i+2] = 4*i+3;
        indices[6*i+3] = 4*i+0;
        indices[6*i+4] = 4*i+3;
        indices[6*i+5] = 4*i+2;
    }
    m_volume_index_buffer->commit(indices);

    m_volume_vertex_buffers.clear();
    CHECK_MSTATUS(m_volume_vertex_buffers.addBuffer("pos_model", m_volume_position_buffer.get()));
}

void VDBSubSceneOverride::update(MHWRender::MSubSceneContainer& container, const MHWRender::MFrameContext& frameContext)
{
    if (!m_volume_render_item) {
        if (!initRenderItem()) {
            return;
        }
        updateGeometry(MAX_SLICE_COUNT);

        if (!container.add(m_volume_render_item)) {
            std::cerr << "VDBSubSceneOverride::update: could not add render item." << std::endl;
        }
    }

    MStatus status;
    MFnDagNode node(m_object, &status);
    if (!status) return;

    MDagPath path;
    node.getPath(path);
    m_volume_render_item->setMatrix(&path.inclusiveMatrix());

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
