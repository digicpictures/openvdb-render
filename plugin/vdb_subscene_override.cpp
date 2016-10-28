#include "vdb_subscene_override.h"

#include <new>
#include <random>

#include <maya/MHwGeometryUtilities.h>

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>

#include "progress_bar.h"
#include "volume_sampler.h"
#include "vdb_maya_utils.hpp"


namespace {
    // We have to options to code shaders, either cgfx, which is deprecated since 2012
    // or ogsfx, which is a glslfx like thing and severely underdocumented.
    // I decided to go with ogsfx, that can be reused easier later on in other
    // packages like katana. -- Pal
    // Best example for ogsfx https://knowledge.autodesk.com/search-result/caas/CloudHelp/cloudhelp/2016/ENU/Maya-SDK/files/GUID-94505429-12F9-4F04-A4D9-B80880AD0BA1-htm.html

    // Fun part comes, when maya is not giving you error messages when the shader is invalid.
    // Awesome, right?

    const char* point_cloud_technique = R"ogsfx(
uniform mat4 wvp : WorldViewProjection;

attribute vs_input
{
    vec3 in_position : POSITION;
};

attribute vs_to_ps
{
    vec4 point_color;
};

attribute ps_output
{
    vec4 out_color : COLOR0;
}

GLSLShader VS
{
    void main()
    {
        gl_Position = wvp * vec4(in_position, 1.0);
        vsOut.point_color = vec(0.0, 1.0, 0.0, 1.0);
    }
}

GLSLShader PS
{
    void main()
    {
        out_color = vec4(1.0, 0.0, 0.0, 1.0);
    }
}

technique Main
{
    pass p0
    {
        VertexShader(in vs_input, out vs_to_ps vsOut) = VS;
        PixelShader(in vs_to_ps psIn, out ps_output) = PS;
    }
}
    )ogsfx";

    const MHWRender::MShaderManager* get_shader_manager()
    {
        auto renderer = MHWRender::MRenderer::theRenderer();
        if (renderer == nullptr)
            return nullptr;

        auto shader_manager = renderer->getShaderManager();
        if (shader_manager == nullptr)
            return nullptr;

        return shader_manager;
    }

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

struct GridSpec {
    openvdb::io::File* vdb_file;
    std::string grid_name;
    GridSpec(openvdb::io::File* vdb_file_, const std::string& grid_name_) : vdb_file(vdb_file_), grid_name(grid_name_) {}
};

namespace MHWRender {

    // === Sliced display mode =================================================

    struct Renderable
    {
        MHWRender::MRenderItem* render_item;
        MHWRender::MVertexBufferArray vertex_buffer_array;
        std::unique_ptr<MHWRender::MVertexBuffer> position_buffer;
        std::unique_ptr<MHWRender::MIndexBuffer> index_buffer;

        Renderable() : render_item(nullptr) {}
        void update(MHWRender::MPxSubSceneOverride& subscene_override, const MBoundingBox& bbox);
    };

    class SlicedDisplay
    {
    public:
        SlicedDisplay(MHWRender::MPxSubSceneOverride& parent);
        bool update(MHWRender::MSubSceneContainer& container, const VDBSubSceneOverrideData& data);
        void enable(bool enable);

    private:
        bool initSliceRenderables(MHWRender::MSubSceneContainer& container);
        bool initBBoxRenderable(MHWRender::MSubSceneContainer& container);
        void updateShaderParams(const SliceShaderParams& shader_params);
        void updateDensityVolume(const GridSpec& grid_spec);
        void updateBBoxGeometry(const openvdb::BBoxd& bbox);

        MHWRender::MPxSubSceneOverride& m_parent;
        ShaderPtr m_volume_shader;
        TexturePtr m_volume_texture;
        Renderable m_slices_renderable;
        Renderable m_bbox_renderable;
        bool m_enabled;

        static const std::string s_effect_code;
    };

    // === SubSceneOverride implementation =====================================

    void VDBSubSceneOverride::shader_instance_deleter::operator()(MShaderInstance* p)
    {
        auto shmgr = get_shader_manager();
        if (shmgr != nullptr)
            shmgr->releaseShader(p);
    }

    struct VDBSubSceneOverrideData {
        MBoundingBox bbox;

        MMatrix world_matrix;
        bool is_selected;
        MColor wireframe_color;

        MFloatVector scattering_color;
        MFloatVector attenuation_color;
        MFloatVector emission_color;

        std::string attenuation_channel;
        std::string scattering_channel;
        std::string emission_channel;

        Gradient scattering_gradient;
        Gradient attenuation_gradient;
        Gradient emission_gradient;

        std::string vdb_path;
        openvdb::io::File* vdb_file;
        openvdb::GridBase::ConstPtr scattering_grid;
        openvdb::GridBase::ConstPtr attenuation_grid;
        openvdb::GridBase::ConstPtr emission_grid;

        float point_size;
        float point_jitter;

        int point_skip;
        int update_trigger;
        VDBDisplayMode display_mode;

        std::string sliced_display_channel;
        SliceShaderParams sliced_display_shader_params;

        enum class ChangeSet : unsigned int {
            NO_CHANGES = 0,
            GENERIC_ATTRIBUTE = 1,
            VDB_FILE = 2,
            SLICED_DISPLAY_GRID = 4
        };
        ChangeSet change_set;

        VDBSubSceneOverrideData();
        ~VDBSubSceneOverrideData();
        void clear();
        bool update(const VDBVisualizerData* data, const MObject& obj);
    };

    inline VDBSubSceneOverrideData::ChangeSet& operator|=(VDBSubSceneOverrideData::ChangeSet& lhs, VDBSubSceneOverrideData::ChangeSet rhs)
    {
        return lhs = VDBSubSceneOverrideData::ChangeSet(unsigned(lhs) | unsigned(rhs));
    }
    inline VDBSubSceneOverrideData::ChangeSet operator|(VDBSubSceneOverrideData::ChangeSet lhs, VDBSubSceneOverrideData::ChangeSet rhs)
    {
        lhs |= rhs;
        return lhs;
    }
    inline VDBSubSceneOverrideData::ChangeSet& operator&=(VDBSubSceneOverrideData::ChangeSet& lhs, VDBSubSceneOverrideData::ChangeSet rhs)
    {
        return lhs = VDBSubSceneOverrideData::ChangeSet(unsigned(lhs) & unsigned(rhs));
    }
    inline VDBSubSceneOverrideData::ChangeSet operator&(VDBSubSceneOverrideData::ChangeSet lhs, VDBSubSceneOverrideData::ChangeSet rhs)
    {
        lhs &= rhs;
        return lhs;
    }

    VDBSubSceneOverrideData::VDBSubSceneOverrideData() : point_size(std::numeric_limits<float>::infinity()),
        point_jitter(std::numeric_limits<float>::infinity()),
        point_skip(-1),
        update_trigger(-1),
        display_mode(DISPLAY_AXIS_ALIGNED_BBOX),
        change_set(ChangeSet::NO_CHANGES)
    {
        for (unsigned int x = 0; x < 4; ++x)
        {
            for (unsigned int y = 0; y < 4; ++y)
                world_matrix(x, y) = std::numeric_limits<float>::infinity();
        }
    }

    VDBSubSceneOverrideData::~VDBSubSceneOverrideData()
    {
        clear();
    }

    void VDBSubSceneOverrideData::clear()
    {
        scattering_grid = 0;
        attenuation_grid = 0;
        emission_grid = 0;
        vdb_file = nullptr;
    }

    bool VDBSubSceneOverrideData::update(const VDBVisualizerData* data, const MObject& obj)
    {
        // TODO: we can limit some of the comparisons to the display mode
        // ie, we don't need to compare certain things if we are using the bounding
        // box mode
        change_set = ChangeSet::NO_CHANGES;

        MDagPath dg = MDagPath::getAPathTo(obj);
        const MMatrix inc_world_matrix = dg.inclusiveMatrix();
        change_set |= setup_parameter(world_matrix, inc_world_matrix, ChangeSet::GENERIC_ATTRIBUTE);
        change_set |= setup_parameter(is_selected, isPathSelected(dg), ChangeSet::GENERIC_ATTRIBUTE);
        change_set |= setup_parameter(wireframe_color, MHWRender::MGeometryUtilities::wireframeColor(dg), ChangeSet::GENERIC_ATTRIBUTE);

        if (data == nullptr || update_trigger == data->update_trigger)
            return change_set != ChangeSet::NO_CHANGES;

        update_trigger = data->update_trigger;

        change_set |= setup_parameter(display_mode, data->display_mode, ChangeSet::GENERIC_ATTRIBUTE);
        change_set |= setup_parameter(bbox, data->bbox, ChangeSet::GENERIC_ATTRIBUTE);
        change_set |= setup_parameter(scattering_color, data->scattering_color, ChangeSet::GENERIC_ATTRIBUTE);
        change_set |= setup_parameter(attenuation_color, data->attenuation_color, ChangeSet::GENERIC_ATTRIBUTE);
        change_set |= setup_parameter(emission_color, data->emission_color, ChangeSet::GENERIC_ATTRIBUTE);
        change_set |= setup_parameter(attenuation_channel, data->attenuation_channel, ChangeSet::GENERIC_ATTRIBUTE);
        change_set |= setup_parameter(scattering_channel, data->scattering_channel, ChangeSet::GENERIC_ATTRIBUTE);
        change_set |= setup_parameter(emission_channel, data->emission_channel, ChangeSet::GENERIC_ATTRIBUTE);
        change_set |= setup_parameter(scattering_gradient, data->scattering_gradient, ChangeSet::GENERIC_ATTRIBUTE);
        change_set |= setup_parameter(attenuation_gradient, data->attenuation_gradient, ChangeSet::GENERIC_ATTRIBUTE);
        change_set |= setup_parameter(emission_gradient, data->emission_gradient, ChangeSet::GENERIC_ATTRIBUTE);
        change_set |= setup_parameter(point_skip, data->point_skip, ChangeSet::GENERIC_ATTRIBUTE);

        change_set |= setup_parameter(vdb_path, data->vdb_path, ChangeSet::VDB_FILE);
        vdb_file = data->vdb_file;
        if (!vdb_file)
            clear();

        if (display_mode == VDBDisplayMode::DISPLAY_SLICES) {
            change_set |= setup_parameter(sliced_display_shader_params, data->sliced_display_shader_params, ChangeSet::GENERIC_ATTRIBUTE);
            change_set |= setup_parameter(sliced_display_channel, data->sliced_display_channel, ChangeSet::SLICED_DISPLAY_GRID);
        }

        point_size = data->point_size;
        point_jitter = data->point_jitter; // We can jitter in the vertex shader. Hopefully

        return change_set != ChangeSet::NO_CHANGES;
    }

    // === VDBSubSceneOverride implementation ==================================

    MString VDBSubSceneOverride::registrantId("VDBVisualizerSubSceneOverride");

    MPxSubSceneOverride* VDBSubSceneOverride::creator(const MObject& obj)
    {
        return new VDBSubSceneOverride(obj);
    }

    VDBSubSceneOverride::VDBSubSceneOverride(const MObject& obj) : MPxSubSceneOverride(obj),
        p_data(new VDBSubSceneOverrideData),
        m_sliced_display(new SlicedDisplay(*this))
    {
        m_object = obj;
        MFnDependencyNode dnode(obj);
        p_vdb_visualizer = dynamic_cast<VDBVisualizerShape*>(dnode.userNode());
        auto shmgr = get_shader_manager();
        if (shmgr != nullptr)
        {
            p_point_cloud_shader.reset(shmgr->getEffectsBufferShader(
                point_cloud_technique, static_cast<unsigned int>(strlen(point_cloud_technique)), "Main", 0, 0, false));

            if (p_point_cloud_shader != nullptr)
                p_point_cloud_shader->setIsTransparent(true);
            else
                std::cerr << "Point cloud shader is zero!" << std::endl;
        }
    }

    VDBSubSceneOverride::~VDBSubSceneOverride()
    {
    }

    MHWRender::DrawAPI VDBSubSceneOverride::supportedDrawAPIs() const
    {
#if MAYA_API_VERSION >= 201600
        return kOpenGLCoreProfile | kOpenGL;
#else
        return kOpenGL;
#endif
    }

    void VDBSubSceneOverride::update(MSubSceneContainer& container, const MFrameContext& /*frameContext*/)
    {
        VDBSubSceneOverrideData* data = p_data.get();

        MHWRender::MRenderer* renderer = MHWRender::MRenderer::theRenderer();
        if (renderer == nullptr)
            return;

        const MHWRender::MShaderManager* shader_manager = renderer->getShaderManager();
        if (shader_manager == nullptr)
            return;

        MRenderItem* bounding_box = container.find("bounding_box");
        if (bounding_box == nullptr)
        {
            bounding_box = MHWRender::MRenderItem::Create("bounding_box",
                MRenderItem::NonMaterialSceneItem,
                MGeometry::kLines);
            bounding_box->enable(false);
            bounding_box->setDrawMode(MGeometry::kAll);
            bounding_box->depthPriority(MRenderItem::sDormantWireDepthPriority);

            MHWRender::MShaderInstance* shader = shader_manager->getStockShader(
                MHWRender::MShaderManager::k3dSolidShader, nullptr, nullptr);
            if (shader)
            {
                // Set the color on the shader instance using the parameter interface
                static const float color[] = { 0.0f, 1.0f, 0.0f, 1.0f };
                shader->setParameter("solidColor", color);

                // Assign the shader to the custom render item
                bounding_box->setShader(shader);
            }

            container.add(bounding_box);
        }

        MHWRender::MRenderItem* point_cloud = container.find("point_cloud");
        if (point_cloud == nullptr)
        {
            point_cloud = MHWRender::MRenderItem::Create("point_cloud",
                MHWRender::MGeometry::kPoints,
                MHWRender::MGeometry::kAll,
                false);
            point_cloud->enable(false);
            point_cloud->setDrawMode(MGeometry::kAll);
            point_cloud->depthPriority(MRenderItem::sDormantPointDepthPriority);

            //if (p_point_cloud_shader == nullptr)
            {
                MHWRender::MShaderInstance* shader = shader_manager->getStockShader(
                    MHWRender::MShaderManager::k3dCPVFatPointShader, nullptr, nullptr);
                if (shader)
                    point_cloud->setShader(shader);
            }
            /*else
                point_cloud->setShader(p_point_cloud_shader.get());*/

            container.add(point_cloud);
        }

        bounding_box->setMatrix(&data->world_matrix);
        point_cloud->setMatrix(&data->world_matrix);

        if (data->change_set == VDBSubSceneOverrideData::ChangeSet::NO_CHANGES)
            return;

        const bool file_exists = data->vdb_file != nullptr;

        const static MVertexBufferDescriptor position_buffer_desc("", MGeometry::kPosition, MGeometry::kFloat, 3);
        const static MVertexBufferDescriptor color_buffer_desc("", MGeometry::kColor, MGeometry::kFloat, 4);

        if (!file_exists || data->display_mode <= DISPLAY_GRID_BBOX)
        {
            point_cloud->enable(false);
            bounding_box->enable(true);
            m_sliced_display->enable(false);

            MVertexBufferArray vertex_buffers;
            p_bbox_position.reset(new MVertexBuffer(position_buffer_desc));
            p_bbox_indices.reset(new MIndexBuffer(MGeometry::kUnsignedInt32));

            if ((data->display_mode == DISPLAY_AXIS_ALIGNED_BBOX) || !file_exists)
            {
                MFloatVector* bbox_vertices = reinterpret_cast<MFloatVector*>(p_bbox_position->acquire(8, true));
                MFloatVector min = data->bbox.min();
                MFloatVector max = data->bbox.max();
                bbox_vertices[0] = MFloatVector(min.x, min.y, min.z);
                bbox_vertices[1] = MFloatVector(min.x, max.y, min.z);
                bbox_vertices[2] = MFloatVector(min.x, max.y, max.z);
                bbox_vertices[3] = MFloatVector(min.x, min.y, max.z);
                bbox_vertices[4] = MFloatVector(max.x, min.y, min.z);
                bbox_vertices[5] = MFloatVector(max.x, max.y, min.z);
                bbox_vertices[6] = MFloatVector(max.x, max.y, max.z);
                bbox_vertices[7] = MFloatVector(max.x, min.y, max.z);
                p_bbox_position->commit(bbox_vertices);
                set_bbox_indices(1, p_bbox_indices.get());
            }
            else if (data->display_mode == DISPLAY_GRID_BBOX)
            {
                try
                {
                    if (!data->vdb_file->isOpen())
                        data->vdb_file->open(false);
                    openvdb::GridPtrVecPtr grids = data->vdb_file->readAllGridMetadata();
                    if (grids->size() == 0)
                        return;
                    std::vector<MFloatVector> vertices;
                    vertices.reserve(grids->size() * 8);

                    for (openvdb::GridPtrVec::const_iterator it = grids->begin(); it != grids->end(); ++it)
                    {
                        if (openvdb::GridBase::ConstPtr grid = *it)
                        {
                            std::array<MFloatVector, 8> _vertices;
                            if (read_grid_transformed_bbox_wire(grid, _vertices))
                            {
                                for (int v = 0; v < 8; ++v)
                                    vertices.push_back(_vertices[v]);
                            }
                        }
                    }

                    const unsigned int vertex_count = static_cast<unsigned int>(vertices.size());

                    if (vertex_count > 0)
                    {
                        p_bbox_position.reset(new MVertexBuffer(position_buffer_desc));
                        p_bbox_indices.reset(new MIndexBuffer(MGeometry::kUnsignedInt32));
                        MFloatVector* bbox_vertices = reinterpret_cast<MFloatVector*>(p_bbox_position->acquire(vertex_count, true));
                        for (unsigned int i = 0; i < vertex_count; ++i)
                            bbox_vertices[i] = vertices[i];
                        p_bbox_position->commit(bbox_vertices);
                        set_bbox_indices(vertex_count / 8, p_bbox_indices.get());
                    }
                }
                catch (...)
                {
                }
            }

            vertex_buffers.addBuffer("", p_bbox_position.get());
            setGeometryForRenderItem(*bounding_box, vertex_buffers, *p_bbox_indices.get(), &data->bbox);
        }
        else
        {
            bounding_box->enable(false);
            if (data->display_mode == DISPLAY_POINT_CLOUD)
            {
                try {
                    if (!data->vdb_file->isOpen())
                        data->vdb_file->open(false);
                    if (data->attenuation_grid == nullptr || data->attenuation_grid->getName() != data->attenuation_channel)
                        data->attenuation_grid = data->vdb_file->readGrid(data->attenuation_channel);
                }
                catch (...) {
                    data->attenuation_grid = nullptr;
                    data->scattering_grid = nullptr;
                    data->emission_grid = nullptr;
                    return;
                }

                point_cloud->enable(true);
                m_sliced_display->enable(false);

                const openvdb::Vec3d voxel_size = data->attenuation_grid->voxelSize();

                FloatVoxelIterator* iter = nullptr;

                if (data->attenuation_grid->valueType() == "float")
                    iter = new FloatToFloatVoxelIterator(openvdb::gridConstPtrCast<openvdb::FloatGrid>(data->attenuation_grid));
                else if (data->attenuation_grid->valueType() == "vec3s")
                    iter = new Vec3SToFloatVoxelIterator(openvdb::gridConstPtrCast<openvdb::Vec3SGrid>(data->attenuation_grid));
                else
                    iter = new FloatVoxelIterator();

                // setting up vertex buffers
                std::vector<MFloatVector> vertices;
                vertices.reserve(iter->get_active_voxels());
                const openvdb::math::Transform attenuation_transform = data->attenuation_grid->transform();

                std::mt19937 mt_generator;
                std::uniform_real_distribution<float> uniform_0_1_dist(0.0f, 1.0f);
                const float point_skip_ratio = 1.0f / static_cast<float>(std::max(data->point_skip, 1));
                for (; iter->is_valid(); iter->get_next())
                {
                    // this gives a better distribution than skipping based on index
                    if (uniform_0_1_dist(mt_generator) > point_skip_ratio)
                        continue;
                    openvdb::Vec3d vdb_pos = attenuation_transform.indexToWorld(iter->get_coord());
                    vertices.push_back(MFloatVector(static_cast<float>(vdb_pos.x()), static_cast<float>(vdb_pos.y()),
                        static_cast<float>(vdb_pos.z())));
                }

                vertices.shrink_to_fit();
                const unsigned int vertex_count = static_cast<unsigned int>(vertices.size());

                if (vertex_count == 0)
                    return;

                delete iter;

                tbb::task_scheduler_init task_init;

                p_index_buffer.reset(new MIndexBuffer(MGeometry::kUnsignedInt32));
                unsigned int* indices = reinterpret_cast<unsigned int*>(p_index_buffer->acquire(vertex_count, true));
                for (unsigned int i = 0; i < vertex_count; ++i)
                    indices[i] = i;
                p_index_buffer->commit(indices);


                p_position_buffer.reset(new MVertexBuffer(position_buffer_desc));
                MFloatVector* pc_vertices = reinterpret_cast<MFloatVector*>(p_position_buffer->acquire(static_cast<unsigned int>(vertices.size()), true));
                const float point_jitter = data->point_jitter;
                if (point_jitter > 0.001f)
                {
                    std::uniform_real_distribution<float> distributionX(-point_jitter * static_cast<float>(voxel_size.x()), point_jitter * static_cast<float>(voxel_size.x()));
                    std::uniform_real_distribution<float> distributionY(-point_jitter * static_cast<float>(voxel_size.y()), point_jitter * static_cast<float>(voxel_size.y()));
                    std::uniform_real_distribution<float> distributionZ(-point_jitter * static_cast<float>(voxel_size.z()), point_jitter * static_cast<float>(voxel_size.z()));

                    tbb::parallel_for(tbb::blocked_range<unsigned int>(0, vertex_count), [&](const tbb::blocked_range<unsigned int>& r) {
                        std::minstd_rand generatorX(42 + r.begin()); // LCG
                        std::minstd_rand generatorY(137 + r.begin());
                        std::minstd_rand generatorZ(1337 + r.begin());
                        for (unsigned int i = r.begin(); i != r.end(); ++i)
                        {
                            MFloatVector pos = vertices[i];
                            pos.x += distributionX(generatorX);
                            pos.y += distributionY(generatorY);
                            pos.z += distributionZ(generatorZ);
                            pc_vertices[i] = pos;
                        }
                    });
                }
                else
                {
                    tbb::parallel_for(tbb::blocked_range<unsigned int>(0, vertex_count), [&](const tbb::blocked_range<unsigned int>& r) {
                        memcpy(&pc_vertices[r.begin()], &vertices[r.begin()], (r.end() - r.begin()) * sizeof(MFloatVector));
                    });
                }
                p_position_buffer->commit(pc_vertices);

                // setting up color buffers

                try {
                    if (data->scattering_grid == nullptr || data->scattering_grid->getName() != data->scattering_channel)
                        data->scattering_grid = data->vdb_file->readGrid(data->scattering_channel);
                }
                catch (...) {
                    data->scattering_grid = nullptr;
                }

                RGBSampler* scattering_sampler = nullptr;

                if (data->scattering_grid == nullptr)
                    scattering_sampler = new (m_scattering_sampler.data()) RGBSampler();
                else
                {
                    if (data->scattering_grid->valueType() == "float")
                        scattering_sampler = new (m_scattering_sampler.data()) FloatToRGBSampler(openvdb::gridConstPtrCast<openvdb::FloatGrid>(data->scattering_grid));
                    else if (data->scattering_grid->valueType() == "vec3s")
                        scattering_sampler = new (m_scattering_sampler.data()) Vec3SToRGBSampler(openvdb::gridConstPtrCast<openvdb::Vec3SGrid>(data->scattering_grid));
                    else
                        scattering_sampler = new (m_scattering_sampler.data()) RGBSampler();
                }

                try {
                    if (data->emission_grid == nullptr || data->emission_grid->getName() != data->emission_channel)
                        data->emission_grid = data->vdb_file->readGrid(data->emission_channel);
                }
                catch (...) {
                    data->emission_grid = nullptr;
                }

                RGBSampler* emission_sampler = nullptr;

                if (data->emission_grid == nullptr)
                    emission_sampler = new (m_emission_sampler.data()) RGBSampler(MFloatVector(0.0f, 0.0f, 0.0f));
                else
                {
                    if (data->emission_grid->valueType() == "float")
                        emission_sampler = new (m_emission_sampler.data()) FloatToRGBSampler(openvdb::gridConstPtrCast<openvdb::FloatGrid>(data->emission_grid));
                    else if (data->emission_grid->valueType() == "vec3s")
                        emission_sampler = new (m_emission_sampler.data()) Vec3SToRGBSampler(openvdb::gridConstPtrCast<openvdb::Vec3SGrid>(data->emission_grid));
                    else
                        emission_sampler = new (m_emission_sampler.data()) RGBSampler(MFloatVector(0.0f, 0.0f, 0.0f));
                }

                RGBSampler* attenuation_sampler = 0;

                if (data->attenuation_grid->valueType() == "float")
                    attenuation_sampler = new (m_attenuation_sampler.data()) FloatToRGBSampler(openvdb::gridConstPtrCast<openvdb::FloatGrid>(data->attenuation_grid));
                else if (data->attenuation_grid->valueType() == "vec3s")
                    attenuation_sampler = new (m_attenuation_sampler.data()) Vec3SToRGBSampler(openvdb::gridConstPtrCast<openvdb::Vec3SGrid>(data->attenuation_grid));
                else
                    attenuation_sampler = new (m_attenuation_sampler.data()) RGBSampler(MFloatVector(1.0f, 1.0f, 1.0f));

                p_color_buffer.reset(new MVertexBuffer(color_buffer_desc));
                MColor* colors = reinterpret_cast<MColor*>(p_color_buffer->acquire(vertex_count, true));
                tbb::parallel_for(tbb::blocked_range<unsigned int>(0, vertex_count), [&](const tbb::blocked_range<unsigned int>& r) {
                    for (unsigned int i = r.begin(); i < r.end(); ++i)
                    {
                        const MFloatVector& v = vertices[i];
                        const openvdb::Vec3d pos(v.x, v.y, v.z);
                        MColor& color = colors[i];
                        const MFloatVector scattering_color = data->scattering_gradient.evaluate(scattering_sampler->get_rgb(pos));
                        const MFloatVector emission_color = data->emission_gradient.evaluate(emission_sampler->get_rgb(pos));
                        const MFloatVector attenuation_color = data->attenuation_gradient.evaluate(attenuation_sampler->get_rgb(pos));
                        color.r = scattering_color.x * data->scattering_color.x + emission_color.x * data->emission_color.x;
                        color.g = scattering_color.y * data->scattering_color.y + emission_color.y * data->emission_color.y;
                        color.b = scattering_color.z * data->scattering_color.z + emission_color.z * data->emission_color.z;
                        color.a = (attenuation_color.x * data->attenuation_color.x +
                            attenuation_color.y * data->attenuation_color.y +
                            attenuation_color.z * data->attenuation_color.z) / 3.0f;
                    }
                });
                p_color_buffer->commit(colors);

                scattering_sampler->~RGBSampler();
                emission_sampler->~RGBSampler();
                attenuation_sampler->~RGBSampler();

                MVertexBufferArray vertex_buffers;
                vertex_buffers.addBuffer("", p_position_buffer.get());
                vertex_buffers.addBuffer("", p_color_buffer.get());

                /*const MVertexBufferDescriptorList& reqs = point_cloud->requiredVertexBuffers();
                for (int i = 0; i < reqs.length(); ++i)
                {
                    MVertexBufferDescriptor desc;
                    reqs.getDescriptor(i, desc);
                    std::cerr << "Requirement : (" << desc.name() << ") " << desc.semantic() << std::endl;
                }*/
                setGeometryForRenderItem(*point_cloud, vertex_buffers, *p_index_buffer.get(), &data->bbox);
            }
            else if (data->display_mode == DISPLAY_SLICES)
            {
                point_cloud->enable(false);

                m_sliced_display->enable(true);
                m_sliced_display->update(container, *data);
            }
        }

        data->change_set = VDBSubSceneOverrideData::ChangeSet::NO_CHANGES;
    }

    bool VDBSubSceneOverride::requiresUpdate(const MSubSceneContainer& /*container*/, const MFrameContext& /*frameContext*/) const
    {
        return p_data->update(p_vdb_visualizer->get_update(), m_object);
    }

    // === Sliced display mode implementation ===================================

    constexpr unsigned int MAX_SLICE_COUNT = 64;

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

    const std::string SlicedDisplay::s_effect_code = R"cgfx(
float3 view_dir_world : ViewDirection < string Space = "World"; >;
float4x4 world_mat : World < string UIWidget = "None"; >;
float4x4 world_view_proj_mat : WorldViewProjection < string UIWidget = "None"; >;

float3 volume_size_model < string UIWidget = "None"; >;   // Size of volume in model space.
float3 volume_origin_model < string UIWidget = "None"; >; // Origin of volume in model space.
texture volume_texture < string UIWidget = "None"; string TextureType = "3D"; >;
sampler3D volume_sampler = sampler_state {
    Texture = <volume_texture>;
};

uniform int max_slice_count; // How many slices does the vertex buffer has vertices for.
uniform float slice_size_model;
uniform float min_voxel_value = 0;
uniform float max_voxel_value = 1;

uniform float3 light_dir = float3(0.3, 0.3, 0);
uniform float3 light_color = float3(0.7, 0.7, 0.7);
uniform float3 scattering = float3(1.5, 1.5, 1.5);
uniform float3 absorption = float3(0.1, 0.1, 0.1);
uniform float shadow_gain = 0.2;
uniform int shadow_sample_count = 4;

float SampleVolume(float3 tex_coords, float lod)
{
    return lerp(min_voxel_value, max_voxel_value, tex3Dlod(volume_sampler, float4(tex_coords, lod)).r);
}

// ======== VERTEX SHADER ========

struct VERT_INPUT
{
    float3 dummy_pos : Position;
    int id : VERTEXID;
};

struct VERT_OUTPUT
{
    float4 pos_clip : Position;
    float3 pos_model : Texcoord0;
    float3 slice_vector_world : Texcoord1;
};

VERT_OUTPUT VolumeVertexShader(VERT_INPUT input)
{
    VERT_OUTPUT output;

    int slice_idx = input.id >> 2;
    // Use dummy attrib to prevent the compiler optimizing it away.
    slice_idx += int(0.000001 * input.dummy_pos.x);
    float2 pos_slice = float2(input.id & 1, (input.id >> 1) & 1);
    float3 pos_model;

    // Define a macro for position calculation.
    // Slice axis can be x, y or z, other_axes can be yz, zx or xy.
#define ALIGN_SLICES(slice_axis, other_axes) do { \
        /* Reorder slices if needed for correct transparency. */ \
        if (slice_axis > 0) slice_idx = (max_slice_count - 1) - slice_idx; \
        /* Constrain slice distance to get at most max_slice_count slices. */ \
        float slice_distance = max(slice_size_model, volume_size_model.slice_axis / (max_slice_count - 1)); \
        /* Calculate model position so that a slice will always intersect the origin. */ \
        int slice_idx_ofs = -int(-(volume_origin_model.slice_axis / slice_distance)); \
        pos_model.slice_axis = (slice_idx + slice_idx_ofs) * slice_distance; \
        /* Kill vertex if out of bounding box. */ \
        if (pos_model.slice_axis < volume_origin_model.slice_axis || \
            pos_model.slice_axis > volume_origin_model.slice_axis + volume_size_model.slice_axis) { \
            output.pos_clip = float4(0, 0, 0, 0); \
            return output; \
        } \
        /* Model position along non-slice axes. */ \
        pos_model.other_axes = pos_slice * volume_size_model.other_axes + volume_origin_model.other_axes; \
        /* The slice vector points toward the slice axis and has slice distance magnitude. */ \
        output.slice_vector_world = local_##slice_axis##_world * slice_distance; \
    } while (false) // Stupid but classic macro trick.

    // Slice along axis of steepest angle w.r.t. camera.
    float3 local_x_world = normalize(world_mat._m00_m10_m20);
    float3 local_y_world = normalize(world_mat._m01_m11_m21);
    float3 local_z_world = normalize(world_mat._m02_m12_m22);
    float x = dot(local_x_world, view_dir_world),
          y = dot(local_y_world, view_dir_world),
          z = dot(local_z_world, view_dir_world);
    float ax = abs(x), ay = abs(y), az = abs(z);
    if (ax > ay && ax > az) {
        ALIGN_SLICES(x, yz);
    } else if (ay > ax && ay > az) {
        ALIGN_SLICES(y, zx);
    } else {
        ALIGN_SLICES(z, xy);
    }

    // Texture coordinates are based on model space position, so that slicing axis doesn't affect volume orientation.
    output.pos_model = pos_model;
    output.pos_clip = mul(world_view_proj_mat, float4(output.pos_model, 1));

    return output;
}

// ======== FRAGMENT SHADER ========

#define OOR(x) ((x) < 0.0f || (x) > 1.0f)
#define MAXV(v) max(max(v.x, v.y), v.z)

float3 shadow_raymarch(float3 pos, float3 dir, float3 extinction)
{
    float shadow_step_model = MAXV(volume_size_model) / float(shadow_sample_count);
    int3 volume_size_voxels = tex3Dsize(volume_sampler, 0);
    float3 step_voxels = volume_size_voxels / float(shadow_sample_count);
    float step_voxels_max = MAXV(step_voxels);
    float lod = log(step_voxels_max) / log(2.f);

    float3 step_texcoords = dir / float(shadow_sample_count);

    float3 transmittance = float3(1, 1, 1);
    float3 sample_texcoords = (pos - volume_origin_model) / volume_size_model;
    for (int i = 0; i < shadow_sample_count; ++i) {
        if (OOR(sample_texcoords.x) || OOR(sample_texcoords.y) || OOR(sample_texcoords.z)) break;
        float density = SampleVolume(sample_texcoords, lod);
        transmittance *= exp(extinction * density * shadow_step_model / (1.0 + float(i) * shadow_step_model * shadow_gain));
        sample_texcoords += step_texcoords;
    }
    return transmittance;
}

typedef VERT_OUTPUT FRAG_INPUT;

struct FRAG_OUTPUT
{
    float4 color : Color0;
};

FRAG_OUTPUT VolumeFragmentShader(FRAG_INPUT input)
{
    FRAG_OUTPUT output;

    float3 extinction = -scattering - absorption;
    float3 albedo = light_color * scattering / (scattering + absorption);

    float3 light_dir_norm = normalize(light_dir);
    float3 shadow = shadow_raymarch(input.pos_model, -light_dir_norm, extinction);
    float3 lumi = albedo * shadow;

    float3 tex_coord = (input.pos_model - volume_origin_model) / volume_size_model;
    float density = SampleVolume(tex_coord, 0);
    float slice_thickness = dot(input.slice_vector_world, input.slice_vector_world) / abs(dot(input.slice_vector_world, view_dir_world));
    float3 transmittance = exp(extinction * density * slice_thickness);

    output.color = float4(lumi, 1 - dot(transmittance, float3(1, 1, 1)/3));

    return output;
}

technique Main < int isTransparent = 1; >
{
    pass P0
    {
        BlendEnable = true;
        BlendFunc = int2(SrcAlpha, OneMinusSrcAlpha);
        CullFaceEnable = false;
        VertexProgram = compile gp5vp VolumeVertexShader();
        FragmentProgram = compile gp5fp VolumeFragmentShader();
    }
}
)cgfx";

    } // unnamed namespace

    SlicedDisplay::SlicedDisplay(MHWRender::MPxSubSceneOverride& parent) : m_parent(parent), m_enabled(false)
    {
        const MHWRender::MShaderManager* shader_manager = getShaderManager();
        assert(shader_manager);

        // Load volume shader from effect file.
        m_volume_shader.reset(shader_manager->getEffectsBufferShader(s_effect_code.c_str(), s_effect_code.size(), "Main", 0, 0, false));
        if (!m_volume_shader) {
            LOG_ERROR("Cannot compile cgfx.");
            return;
        }
        m_volume_shader->setIsTransparent(true);
        m_volume_shader->setParameter("max_slice_count", int(MAX_SLICE_COUNT));
    }

    void SlicedDisplay::enable(bool enable)
    {
        m_enabled = enable;
        if (m_slices_renderable.render_item)
            m_slices_renderable.render_item->enable(enable);
        if (m_bbox_renderable.render_item)
            m_bbox_renderable.render_item->enable(enable);
    }

    bool SlicedDisplay::initSliceRenderables(MHWRender::MSubSceneContainer& container)
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
            positions[4 * i + 0] = MFloatVector(0.0, 0.0, z);
            positions[4 * i + 1] = MFloatVector(1.0, 0.0, z);
            positions[4 * i + 2] = MFloatVector(0.0, 1.0, z);
            positions[4 * i + 3] = MFloatVector(1.0, 1.0, z);
        }
        m_slices_renderable.position_buffer->commit(positions);

        // - Indices
        m_slices_renderable.index_buffer.reset(new MHWRender::MIndexBuffer(MHWRender::MGeometry::kUnsignedInt32));
        const auto index_count = slice_count * 6;
        unsigned int* indices = reinterpret_cast<unsigned int*>(m_slices_renderable.index_buffer->acquire(index_count, true));
        for (unsigned int i = 0; i < slice_count; ++i) {
            indices[6 * i + 0] = 4 * i + 0;
            indices[6 * i + 1] = 4 * i + 1;
            indices[6 * i + 2] = 4 * i + 3;
            indices[6 * i + 3] = 4 * i + 0;
            indices[6 * i + 4] = 4 * i + 3;
            indices[6 * i + 5] = 4 * i + 2;
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

    bool SlicedDisplay::initBBoxRenderable(MHWRender::MSubSceneContainer& container)
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

        auto render_item = MHWRender::MRenderItem::Create(
                "sliced_display_bounding_box",
                MHWRender::MRenderItem::RenderItemType::DecorationItem,
                MHWRender::MGeometry::kLines);
        if (!render_item) {
            LOG_ERROR("Failed to create bbox render item.");
            return false;
        }

        if (!render_item->setShader(shader)) {
            LOG_ERROR("Failed to set shader for bbox render item.");
            return false;
        }

        render_item->setDrawMode(MHWRender::MGeometry::kAll);
        render_item->depthPriority(MHWRender::MRenderItem::sActiveWireDepthPriority);

        // Note: descriptor name (first ctor arg) MUST be "", or setGeometryForRenderItem will return kFailure.
        const MHWRender::MVertexBufferDescriptor pos_desc("", MHWRender::MGeometry::kPosition, MHWRender::MGeometry::kFloat, 3);
        m_bbox_renderable.position_buffer.reset(new MHWRender::MVertexBuffer(pos_desc));
        m_bbox_renderable.vertex_buffer_array.clear();
        CHECK_MSTATUS(m_bbox_renderable.vertex_buffer_array.addBuffer("pos_model", m_bbox_renderable.position_buffer.get()));

        m_bbox_renderable.index_buffer.reset(new MHWRender::MIndexBuffer(MHWRender::MGeometry::kUnsignedInt32));
        constexpr auto index_count = 2 * 12;
        static const uint32_t BOX_WIREFRAME_INDICES[] = { 0, 1, 1, 3, 3, 2, 2, 0, 4, 5, 5, 7, 7, 6, 6, 4, 0, 4, 1, 5, 3, 7, 2, 6 };
        CHECK_MSTATUS(m_bbox_renderable.index_buffer->update(BOX_WIREFRAME_INDICES, 0, index_count, true));

        if (!container.add(render_item)) {
            LOG_ERROR("Could not add bbox render item.");
            m_bbox_renderable.position_buffer.reset();
            m_bbox_renderable.index_buffer.reset();
            m_bbox_renderable.vertex_buffer_array.clear();
            return false;
        }
        m_bbox_renderable.render_item = render_item;

        return true;
    }

    void SlicedDisplay::updateBBoxGeometry(const openvdb::BBoxd& bbox)
    {
        constexpr auto vertex_count = 8;
        MFloatVector* positions = static_cast<MFloatVector*>(m_bbox_renderable.position_buffer->acquire(vertex_count, true));
        for (unsigned int i = 0; i < vertex_count; ++i) {
            positions[i] = mayavecFromVec3f(openvdb::Vec3d(i & 1, (i & 2) >> 1, (i & 4) >> 2) * bbox.extents() + bbox.min());
        }
        m_bbox_renderable.position_buffer->commit(positions);

        m_bbox_renderable.update(m_parent, mayabboxFromBBoxd(bbox));
    }

    void SlicedDisplay::updateShaderParams(const SliceShaderParams& params)
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

        openvdb::FloatGrid::ConstPtr loadDensityGrid(const GridSpec& grid_spec)
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

    } // unnamed namespace

    void SlicedDisplay::updateDensityVolume(const GridSpec& grid_spec)
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
        m_slices_renderable.update(m_parent, mayabboxFromBBoxd(bbox_ws));
        updateBBoxGeometry(bbox_ws);

        // Sample the multi resolution grid at regular intervals.
        ProgressBar pb("vdb_visualizer: sampling density grid");
        VolumeSampler volume_sampler;
        const auto texture_extents = openvdb::Coord(MAX_SLICE_COUNT, MAX_SLICE_COUNT, MAX_SLICE_COUNT);
        auto volume = volume_sampler.sampleGridWithMipmapFilter(*grid_ptr, texture_extents, &pb);
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

    } // unnamed namespace

    bool SlicedDisplay::update(MHWRender::MSubSceneContainer& container, const VDBSubSceneOverrideData& data)
    {
        if (m_slices_renderable.render_item == nullptr) {
            if (!initSliceRenderables(container)) {
                return false;
            }
        }

        if (m_bbox_renderable.render_item == nullptr) {
            if (!initBBoxRenderable(container)) {
                return false;
            }
        }

        // Handle selection.
        m_bbox_renderable.render_item->enable(data.is_selected && m_enabled);
        m_slices_renderable.render_item->enable(m_enabled);

        // Set wireframe color.
        const auto& color = data.wireframe_color;
        const float color_as_array[] = { color.r, color.g, color.b, color.a };
        CHECK_MSTATUS(m_bbox_renderable.render_item->getShader()->setParameter("solidColor", color_as_array));

        // Set world matrix.
        m_slices_renderable.render_item->setMatrix(&data.world_matrix);
        m_bbox_renderable.render_item->setMatrix(&data.world_matrix);

        updateShaderParams(data.sliced_display_shader_params);

        typedef VDBSubSceneOverrideData::ChangeSet ChangeSet;
        const auto rebuild_volume_change_mask = ChangeSet::VDB_FILE | ChangeSet::SLICED_DISPLAY_GRID;
        if ((data.change_set & rebuild_volume_change_mask) != ChangeSet::NO_CHANGES)
            updateDensityVolume(GridSpec(data.vdb_file, data.sliced_display_channel));

        return true;
    }

    void Renderable::update(MHWRender::MPxSubSceneOverride& subscene_override, const MBoundingBox& bbox)
    {
        // Note: render item has to be added to the MSubSceneContainer before calling setGeometryForRenderItem.
        CHECK_MSTATUS(subscene_override.setGeometryForRenderItem(*render_item, vertex_buffer_array, *index_buffer, &bbox));
    }
}
