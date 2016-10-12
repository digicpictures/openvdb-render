#include "vdb_subscene_override.h"

#include <cmath>
#include <chrono>

#include <boost/filesystem.hpp>

#include <maya/MFnDagNode.h>
#include <maya/MGlobal.h>
#include <maya/MSelectionList.h>
#include <maya/MShaderManager.h>

#include <openvdb/tools/MultiResGrid.h>

#include "vdb_visualizer.h"
#include "volume_sampler.h"

// Disable "decorated name length exceeded, name was truncated" warning.
#pragma warning(disable: 4503)

using namespace MHWRender;

#define NUM_SLICES 80

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
        // TODO: compile in effect source.
        using boost::filesystem::path;
        path effect_file = path(__FILE__).parent_path() / "volume.cgfx";
        m_volume_shader.reset(shader_manager->getEffectsFileShader(effect_file.c_str(), "Main", 0, 0, false));
        if (!m_volume_shader) {
            static bool first_time = true;
            if (first_time) {
                std::cerr << "Cannot load cgfx file." << std::endl;
                first_time = false;
            }
            return false;
        }
        m_volume_shader->setIsTransparent(true);
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

void VDBSubSceneOverride::updateShaderParams(const VDBVisualizerData* data)
{
    assert(data);

    const MFloatVector bbox_size = data->bbox.max() - data->bbox.min();
    const MFloatVector bbox_origin = data->bbox.min();
    CHECK_MSTATUS(m_volume_shader->setParameter("volume_size_model", bbox_size));
    CHECK_MSTATUS(m_volume_shader->setParameter("volume_origin_model", bbox_origin));

    auto texture_manager = MRenderer::theRenderer()->getTextureManager();

    // Sample grid
    if (data->vdb_file)
    {
        if (!data->vdb_file->isOpen()) {
            data->vdb_file->open(false);
        }
        auto grid_ptr = openvdb::gridPtrCast<openvdb::FloatGrid>(data->vdb_file->readGrid("density")); // TODO: use cached attribute
        const auto grid_extent = grid_ptr->evalActiveVoxelBoundingBox().extents();
        const size_t levels = openvdb::math::Ceil(std::log2(grid_extent[grid_extent.maxIndex()]));

        auto start_time = std::chrono::steady_clock::now();
        openvdb::tools::MultiResGrid<openvdb::FloatTree> multi_res_grid(levels, grid_ptr);
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time).count();
        std::cout << "MultiResGrid generation time: " << elapsed << "ms." << std::endl;

        // TODO: use multi res grid
        start_time = std::chrono::steady_clock::now();
        m_volume_texture.reset(volumeTextureFromGrid(multi_res_grid, openvdb::Coord(NUM_SLICES, NUM_SLICES, NUM_SLICES), texture_manager));
        elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time).count();
        std::cout << "Texture generation time: " << elapsed << "ms." << std::endl;
    }

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

void VDBSubSceneOverride::updateGeometry(const VDBVisualizerData* data)
{
    assert(data);

    // Create geometry.
    const int num_slices = NUM_SLICES; // TODO: get from node attribute
    const int num_vertices = num_slices * 4;
    const int num_indices = num_slices * 6;

    // - Vertices
    // Note: descriptor name (first ctor arg) MUST be "", or setGeometryForRenderItem will return kFailure.
    const MVertexBufferDescriptor pos_desc("", MGeometry::kPosition, MGeometry::kFloat, 3);
    m_volume_position_buffer.reset(new MVertexBuffer(pos_desc));
    MFloatVector* positions = reinterpret_cast<MFloatVector*>(m_volume_position_buffer->acquire(static_cast<unsigned int>(num_vertices), false));
    for (int i = 0; i < num_slices; ++i)
    {
        const float z = i * 1.0f / (num_slices - 1.0f);
        positions[4*i+0] = MFloatVector(0.0, 0.0, z);
        positions[4*i+1] = MFloatVector(1.0, 0.0, z);
        positions[4*i+2] = MFloatVector(0.0, 1.0, z);
        positions[4*i+3] = MFloatVector(1.0, 1.0, z);
    }
    m_volume_position_buffer->commit(positions);

    // - Indices
    m_volume_index_buffer.reset(new MIndexBuffer(MGeometry::kUnsignedInt32));
    unsigned int* indices = reinterpret_cast<unsigned int*>(m_volume_index_buffer->acquire(num_indices, false));
    for (unsigned int i = 0; i < num_slices; ++i) {
        indices[6*i+0] = 4*i+0;
        indices[6*i+1] = 4*i+1;
        indices[6*i+2] = 4*i+3;
        indices[6*i+3] = 4*i+0;
        indices[6*i+4] = 4*i+3;
        indices[6*i+5] = 4*i+2;
    }
    m_volume_index_buffer->commit(indices);

    MVertexBufferArray volume_vertex_buffers;
    CHECK_MSTATUS(volume_vertex_buffers.addBuffer("pos_model", m_volume_position_buffer.get()));
    CHECK_MSTATUS(setGeometryForRenderItem(*m_volume_render_item, volume_vertex_buffers, *m_volume_index_buffer, &data->bbox));
}

void VDBSubSceneOverride::update(MHWRender::MSubSceneContainer& container, const MHWRender::MFrameContext& frameContext)
{
    if (!m_volume_render_item) {
        if (!initRenderItem()) {
            return;
        }

        if (!container.add(m_volume_render_item)) {
            std::cerr << "VDBSubSceneOverride::update: could not add render item." << std::endl;
        }
    }

    MStatus status;
    MFnDagNode node(m_object, &status);
    if (!status) return;

    // Get VDB shape object from associated node.
    typedef VDBVisualizerShape ShapeType;
    ShapeType* shape_node = dynamic_cast<ShapeType*>(node.userNode());
    if (!shape_node) return;

    auto data = shape_node->get_update();
    if (data) {
        updateShaderParams(data);
        updateGeometry(data);
    }

    // TODO: handle selection.
    MSelectionList selected_list;
    MGlobal::getActiveSelectionList(selected_list);

    // Handle instancing.
    MDagPathArray instances;
    if (!node.getAllPaths(instances) || instances.length() == 0) {
        return;
    }

    auto is_any_ancestor_selected = [&selected_list](MDagPath path) {
        MStatus status;
        do {
            if (selected_list.hasItem(path)) {
                return true;
            }
            status = path.pop();
        } while (status);
        return false;
    };

    MMatrixArray instance_transform_array(instances.length());
    int num_instances = 0;
    bool any_instance_changed = false;
    for (auto i = 0u; i < instances.length(); ++i) {
        auto& instance = instances[i];
        if (!instance.isValid() || !instance.isVisible()) {
            continue;
        }

        // Add or update instance info.
        InstanceInfo instance_info(instance.inclusiveMatrix(), is_any_ancestor_selected(instance));
        const int instance_number = instance.instanceNumber();
        InstanceInfoMap::iterator iter = m_instance_info_cache.find(instance_number);
        if (iter == m_instance_info_cache.end()) {
            m_instance_info_cache[instance_number] = instance_info;
            any_instance_changed = true;
        } else if (iter->second != instance_info) {
            iter->second = instance_info;
            any_instance_changed = true;
        }

        instance_transform_array[num_instances++] = instance_info.m_transform;
    }
    instance_transform_array.setLength(num_instances);

    if (m_num_instances == 0 && num_instances == 1) {
        // First instance.
        addInstanceTransform(*m_volume_render_item, instance_transform_array[0]);
        m_num_instances = num_instances;
    } else if (num_instances > m_num_instances) {
        // Multiple instances are not supported (yet?).
        std::cerr << "VDBVisualizer node does not support instancing." << std::endl;
        m_num_instances = num_instances;
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
