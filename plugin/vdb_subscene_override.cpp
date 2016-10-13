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

#define NUM_SLICES 64

using namespace MHWRender;

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
        m_volume_shader->setParameter("max_slice_count", NUM_SLICES);
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

    // Shading parameters.
    m_volume_shader->setParameter("light_dir",   data->light_direction);
    auto v = data->light_direction;
    m_volume_shader->setParameter("light_color", data->light_color);
    m_volume_shader->setParameter("scattering",  data->scattering);
    m_volume_shader->setParameter("absorption",  data->absorption);
    m_volume_shader->setParameter("shadow_gain", data->shadow_gain);
    m_volume_shader->setParameter("shadow_sample_count", data->shadow_sample_count);
    m_volume_shader->setParameter("slice_size_model", data->slice_size);

    // Bail if no vdb file is specified.
    if (!data->vdb_file)
    {
        m_volume_texture.reset();
        MHWRender::MTextureAssignment volume_texture_resource;
        volume_texture_resource.texture = nullptr;
        m_volume_shader->setParameter("volume_texture", volume_texture_resource);
        return;
    }

    // Check if we need to recreate the volume texture.

    struct VDBChannel {
        std::string m_file_name;
        std::string m_density_channel_name;

        VDBChannel(const std::string& file_name = "", const std::string& density_channel_name = "")
            : m_file_name(file_name), m_density_channel_name(density_channel_name) {}
        bool operator==(const VDBChannel& rhs) const { return m_file_name == rhs.m_file_name && m_density_channel_name == rhs.m_density_channel_name; }
    };

    VDBChannel channel(data->vdb_path, data->density_channel);
    static VDBChannel cache;
    if (cache == channel) {
        return;
    }
    cache = channel;

    // Load channel data.
    openvdb::FloatGrid::ConstPtr grid_ptr;
    try {
        if (!data->vdb_file->isOpen()) {
            data->vdb_file->open(false);
        }
        // TODO: uncomment
        grid_ptr = openvdb::gridPtrCast<openvdb::FloatGrid>(data->vdb_file->readGrid(data->density_channel));
        //grid_ptr = openvdb::gridPtrCast<openvdb::FloatGrid>(data->vdb_file->readGrid("density"));
    } catch (const openvdb::Exception& e) {
        std::cerr << "error reading file " << data->vdb_path << " : " << e.what() << std::endl;
        return;
    }

    // TODO
    auto texture_extents = openvdb::Coord(NUM_SLICES, NUM_SLICES, NUM_SLICES);

    /*
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
    */

    // Sample the density channel grid and create a volume texture from the samples.

    auto start_time = std::chrono::steady_clock::now();

    auto volume_sampler = VolumeSampler(texture_manager);
    auto volume = volume_sampler.sampleGrid(*grid_ptr, texture_extents);
    m_volume_texture.reset(volume.texture);

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time).count();
    std::cout << "Sampling time: " << elapsed << "ms." << std::endl;

    m_volume_shader->setParameter("min_voxel_value", volume.min_value);
    m_volume_shader->setParameter("max_voxel_value", volume.max_value);

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

    MDagPath path;
    node.getPath(path);
    m_volume_render_item->setMatrix(&path.inclusiveMatrix());
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
