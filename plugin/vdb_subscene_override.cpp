#include "vdb_subscene_override.h"

#include <new>
#include <random>
#include <unordered_map>

#include <maya/MHwGeometryUtilities.h>
#include <maya/MDrawContext.h>

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

#define LOG_ERROR(msg) log_error(msg, __FILE__, __LINE__)
    inline void log_error(const std::string& msg, const char *file_name, int line_no)
    {
#if _DEBUG
        std::cerr << "openvdb_render error: " << file_name << ": line " << line_no << ": " << msg << std::endl;
#else
        std::cerr << "openvdb_render error: " << msg << std::endl;
#endif
    }

} // unnamed namespace

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
        operator bool() const { return render_item != nullptr; }
    };

    struct VolumeChannel
    {
        typedef std::shared_ptr<VolumeChannel> Ptr;
        typedef std::shared_ptr<const VolumeChannel> ConstPtr;

        openvdb::FloatGrid::ConstPtr grid;
        openvdb::tools::MultiResGrid<openvdb::FloatTree>::ConstPtr multires;
        VolumeTexture volume_texture;

        VolumeChannel() {}
        VolumeChannel(openvdb::io::File* vdb_file, const std::string& channel_name);
        VolumeChannel(VolumeChannel&&) = default;
        VolumeChannel& operator=(VolumeChannel&&) = default;
        bool isValid() const { return grid.get() != nullptr; }
        void sample(VolumeSampler& volume_sampler, int slice_count);
    };

    class SlicedDisplay
    {
    public:
        SlicedDisplay(MHWRender::MPxSubSceneOverride& parent);
        bool update(MHWRender::MSubSceneContainer& container, const VDBSubSceneOverrideData& data);
        void enable(bool enable);

    private:
        bool initRenderItems(MHWRender::MSubSceneContainer& container, const VDBSubSceneOverrideData& data);

        void updateBBox(const MBoundingBox& bbox);
        void updateSliceGeo(const VDBSubSceneOverrideData& data);

        MHWRender::MPxSubSceneOverride& m_parent;
        ShaderPtr m_volume_shader;

        std::unordered_map<std::string, std::weak_ptr<VolumeChannel>> m_channel_cache;
        VolumeChannel::Ptr m_scattering_channel;
        VolumeChannel::Ptr m_absorption_channel;
        VolumeChannel::Ptr m_emission_channel;
        VolumeSampler m_volume_sampler;

        Renderable m_slices_renderable;
        Renderable m_bbox_renderable;
        const MHWRender::MSamplerState *m_volume_sampler_state;
        bool m_enabled;
        bool m_selected;

        static void preDrawCallback(MHWRender::MDrawContext& context, const MHWRender::MRenderItemList& renderItemList, MHWRender::MShaderInstance* shaderInstance);
        static const std::string s_effect_code;
    };

    // === VDBSubSceneOverrideData implementation =====================================

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
        int max_slice_count;
        SliceShaderParams sliced_display_shader_params;

        enum class ChangeSet : unsigned int {
            NO_CHANGES = 0,
            GENERIC_ATTRIBUTE = 1<<0,
            VDB_FILE =          1<<1,
            MAX_SLICE_COUNT =   1<<2,
            SCATTERING_GRID =   1<<3,
            ATTENUATION_GRID =  1<<4,
            EMISSION_GRID =     1<<5
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
    inline bool hasChange(VDBSubSceneOverrideData::ChangeSet change_set, VDBSubSceneOverrideData::ChangeSet mask)
    {
        return (change_set & mask) != VDBSubSceneOverrideData::ChangeSet::NO_CHANGES;
    }

    VDBSubSceneOverrideData::VDBSubSceneOverrideData() : point_size(std::numeric_limits<float>::infinity()),
        is_selected(false),
        point_jitter(std::numeric_limits<float>::infinity()),
        point_skip(-1),
        update_trigger(-1),
        max_slice_count(-1),
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
        change_set |= setup_parameter(attenuation_channel, data->attenuation_channel, ChangeSet::ATTENUATION_GRID);
        change_set |= setup_parameter(scattering_channel, data->scattering_channel, ChangeSet::SCATTERING_GRID);
        change_set |= setup_parameter(emission_channel, data->emission_channel, ChangeSet::EMISSION_GRID);
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
            change_set |= setup_parameter(max_slice_count, data->max_slice_count, ChangeSet::MAX_SLICE_COUNT);
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

    namespace {

        openvdb::FloatGrid::ConstPtr loadFloatGrid(openvdb::io::File* vdb_file, const std::string& grid_name)
        {
            if (!vdb_file || !vdb_file->isOpen()) {
                return nullptr;
            }

            openvdb::GridBase::ConstPtr grid_base_ptr;
            try {
                grid_base_ptr = vdb_file->readGrid(grid_name);
            }
            catch (const openvdb::Exception& e) {
                std::stringstream ss;
                ss << "Error reading grid " << grid_name << ": " << e.what();
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

    VolumeChannel::VolumeChannel(openvdb::io::File* vdb_file, const std::string& channel_name)
    {
        // Load grid.
        grid = loadFloatGrid(vdb_file, channel_name);
        if (!grid)
            return;

        // Create multires grid.
        const auto grid_extents = getIndexSpaceBoundingBox(grid.get()).extents().asVec3d();
        const auto num_levels = size_t(openvdb::math::Ceil(std::log2(maxComponentValue(grid_extents))));
        multires.reset(new openvdb::tools::MultiResGrid<openvdb::FloatTree>(num_levels, *grid.get()));
    }

    void VolumeChannel::sample(VolumeSampler& volume_sampler, int slice_count)
    {
        ProgressBar pb("vdb_visualizer: sampling grid");
        const auto extents = openvdb::Coord(slice_count, slice_count, slice_count);
        volume_texture = volume_sampler.sampleMultiResGrid(*multires, extents, &pb);
    }


    namespace {

        MFloatVector inline mayavecFromVec3f(const openvdb::Vec3f& vec)
        {
            return { vec.x(), vec.y(), vec.z() };
        }

        MFloatVector inline mayavecFromFloatRange(const FloatRange& float_range)
        {
            return { float_range.min, float_range.max };
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
float4x4 world_inverse_mat : WorldInverse < string UIWidget = "None"; >;
float4x4 world_view_proj_mat : WorldViewProjection < string UIWidget = "None"; >;

// The following sizes and positions are in model space.

float3 volume_size < string UIWidget = "None"; >;
float3 volume_origin < string UIWidget = "None"; >;

float3    scattering_color = float3(1, 1, 1);
float     scattering_anisotropy = 0;
bool      use_scattering_texture = true;
float2    scattering_value_range = float2(0, 1);
float3    scattering_volume_size;
float3    scattering_volume_origin;
texture   scattering_texture < string TextureType = "3D"; >;
sampler3D scattering_sampler = sampler_state {
    Texture = <scattering_texture>;
};

float3    absorption_color = float3(1, 1, 1);
bool      use_absorption_texture = true;
float2    absorption_value_range = float2(0, 1);
float3    absorption_volume_size;
float3    absorption_volume_origin;
texture   absorption_texture < string TextureType = "3D"; >;
sampler3D absorption_sampler = sampler_state {
    Texture = <absorption_texture>;
};

float3    emission_color = float3(1, 1, 1);
bool      use_emission_texture = true;
float2    emission_value_range = float2(0, 1);
float3    emission_volume_size;
float3    emission_volume_origin;
texture   emission_texture < string TextureType = "3D"; >;
sampler3D emission_sampler = sampler_state {
    Texture = <emission_texture>;
};

uniform int max_slice_count; // How many slices does the vertex buffer has vertices for.

// Lights.
#define MAX_POINT_LIGHTS 8
uniform int    point_light_count;
uniform float3 point_light_positions[MAX_POINT_LIGHTS];
uniform float3 point_light_colors[MAX_POINT_LIGHTS];
uniform float  point_light_intensities[MAX_POINT_LIGHTS];

#define MAX_DIRECTIONAL_LIGHTS 8
uniform int    directional_light_count;
uniform float3 directional_light_directions[MAX_DIRECTIONAL_LIGHTS];
uniform float3 directional_light_colors[MAX_DIRECTIONAL_LIGHTS];
uniform float  directional_light_intensities[MAX_DIRECTIONAL_LIGHTS];

uniform float shadow_gain = 0.2;
uniform int shadow_sample_count = 4;

float3 SampleScatteringTexture(float3 pos_model, float lod)
{
    float3 tex_coords = (pos_model - scattering_volume_origin) / scattering_volume_size;
    if (use_scattering_texture)
        return scattering_color * lerp(scattering_value_range.x, scattering_value_range.y, tex3Dlod(scattering_sampler, float4(tex_coords, lod)).r);
    else
        return scattering_color;
}

float3 SampleAbsorptionTexture(float3 pos_model, float lod)
{
    float3 tex_coords = (pos_model - absorption_volume_origin) / absorption_volume_size;
    if (use_absorption_texture)
        return absorption_color * lerp(absorption_value_range.x, absorption_value_range.y, tex3Dlod(absorption_sampler, float4(tex_coords, lod)).r);
    else
        return absorption_color;
}

float3 SampleEmissionTexture(float3 pos_model, float lod)
{
    float3 tex_coords = (pos_model - emission_volume_origin) / emission_volume_size;
    if (use_emission_texture)
        return emission_color * lerp(emission_value_range.x, emission_value_range.y, tex3Dlod(emission_sampler, float4(tex_coords, lod)).r);
    else
        return emission_color;
}

float3 DominantAxis(float3 dir, float3 v)
{
    float ax = abs(dir.x), ay = abs(dir.y), az = abs(dir.z);
    if (ax > ay && ax > az)
        return v.zxy;
    else if (ay > ax && ay > az)
        return v.yzx;
    else
        return v.xyz;
}

// ======== VERTEX SHADER ========

struct VERT_INPUT
{
    float3 pos : Position;
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
    float2 pos_slice = input.pos.xy;

    float3 view_dir_model = normalize(mul(world_inverse_mat, float4(view_dir_world, 0)).xyz);
    if (dot(DominantAxis(view_dir_model, float3(0, 0, 1)), view_dir_model) > 0)
        slice_idx = max_slice_count - 1 - slice_idx;
    float3 pos_dom = float3(pos_slice, float(slice_idx) / float(max_slice_count));
    output.pos_model = DominantAxis(view_dir_model, pos_dom) * volume_size + volume_origin;
    output.pos_clip = mul(world_view_proj_mat, float4(output.pos_model, 1));

    return output;
}

// ======== FRAGMENT SHADER ========

#define OOR(x) ((x) < 0.0f || (x) > 1.0f)
#define MAXV(v) max(max(v.x, v.y), v.z)

float CalcLOD(float distance_model, float3 size_model)
{
    float3 distance_voxels = (distance_model / size_model) * float(max_slice_count);
    return log(MAXV(distance_voxels)) / log(2.f);
}

float3 RayTransmittance(float3 from_model, float3 to_model)
{
    float3 step_model = (to_model - from_model) / float(shadow_sample_count + 1);
    float step_size_model = length(step_model);

    float scattering_lod = CalcLOD(step_size_model, scattering_volume_size);
    float absorption_lod = CalcLOD(step_size_model, absorption_volume_size);

    float3 transmittance = float3(1, 1, 1);
    float3 pos_model = from_model + 0.5f * step_model;
    for (int i = 0; i < shadow_sample_count; ++i) {
        float3 scattering = SampleScatteringTexture(pos_model, scattering_lod);
        float3 absorption = SampleAbsorptionTexture(pos_model, absorption_lod);
        float3 extinction = scattering + absorption;
        transmittance *= exp(-extinction * step_size_model / (1.0 + float(i) * step_size_model * shadow_gain));
        pos_model += step_model;
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

    float3 scattering = SampleScatteringTexture(input.pos_model, 0);
    float3 absorption = SampleAbsorptionTexture(input.pos_model, 0);
    float3 extinction = scattering + absorption;
    float3 albedo = scattering / extinction;

    float3 lumi = float3(0, 0, 0);

    // Loop through directional lights.
    for (int i = 0; i < directional_light_count; ++i) {
        float3 light_vector_norm = -normalize(mul(world_inverse_mat, float4(directional_light_directions[i], 0)).xyz);
        float3 corner = volume_origin + (0.5f + 0.5f * sign(light_vector_norm)) * volume_size;
        float3 light_pos = input.pos_model + 1.1f * length(corner - input.pos_model) * light_vector_norm;
        float3 shadow = RayTransmittance(input.pos_model, light_pos);
        float3 light_color = directional_light_colors[i] * directional_light_intensities[i];
        lumi += albedo * light_color * shadow;
    }

    // Loop through point lights.
    for (int i = 0; i < point_light_count; ++i) {
        float3 light_pos_model = mul(world_inverse_mat, float4(point_light_positions[i], 1)).xyz;
        float3 shadow = RayTransmittance(input.pos_model, light_pos_model);
        float3 light_color = point_light_colors[i] * point_light_intensities[i];
        lumi += albedo * light_color * shadow;
    }

    float3 view_dir_model = normalize(mul(world_inverse_mat, float4(view_dir_world, 0)).xyz);
    float3 dom = DominantAxis(view_dir_model/volume_size, float3(0, 0, 1));
    float slice_thickness = dot(dom, volume_size / float(max_slice_count)) / abs(dot(dom, view_dir_model));

    float3 transmittance = exp(-extinction * slice_thickness);

    float3 emission = SampleEmissionTexture(input.pos_model, 0);
    lumi += emission;

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

    void SlicedDisplay::preDrawCallback(MHWRender::MDrawContext& context, const MHWRender::MRenderItemList& /*renderItemList*/, MHWRender::MShaderInstance* shaderInstance)
    {
        // Collect light data.

        constexpr int MAX_POINT_LIGHTS = 8;
        int point_light_count = 0;
        std::array<float, 3*MAX_POINT_LIGHTS> point_light_positions;
        std::array<float, 3*MAX_POINT_LIGHTS> point_light_colors;
        std::array<float, MAX_POINT_LIGHTS>   point_light_intensities;

        constexpr int MAX_DIRECTIONAL_LIGHTS = 8;
        int directional_light_count = 0;
        std::array<float, 3*MAX_DIRECTIONAL_LIGHTS> directional_light_directions;
        std::array<float, 3*MAX_DIRECTIONAL_LIGHTS> directional_light_colors;
        std::array<float, MAX_DIRECTIONAL_LIGHTS>   directional_light_intensities;

        const auto light_count = context.numberOfActiveLights();
        for (unsigned int i = 0; i < light_count; ++i)
        {
            MIntArray int_array;
            MFloatArray float_array;

            const auto light_params = context.getLightParameterInformation(i);

            // Continue if light is not enabled.
            const auto status = light_params->getParameter(MLightParameterInformation::kLightEnabled, float_array);
            if (status != MS::kSuccess || float_array[0] != 1)
                continue;

            // Get additional info based on light type and save light data.
            const auto light_type = light_params->lightType();
            if (light_type == "pointLight")
            {
                if (point_light_count == MAX_POINT_LIGHTS)
                    continue;

                // Position.
                light_params->getParameter(MLightParameterInformation::kWorldPosition, float_array);
                memcpy(point_light_positions.data() + 3 * point_light_count, &float_array[0], 3 * sizeof(float));

                // Color.
                light_params->getParameter(MLightParameterInformation::kColor, float_array);
                memcpy(point_light_colors.data() + 3 * point_light_count, &float_array[0], 3 * sizeof(float));

                // Intensity.
                light_params->getParameter(MLightParameterInformation::kIntensity, float_array);
                point_light_intensities[point_light_count] = float_array[0];

                ++point_light_count;
            }
            else if (light_type == "directionalLight")
            {
                if (directional_light_count == MAX_DIRECTIONAL_LIGHTS)
                    continue;

                // Direction.
                light_params->getParameter(MLightParameterInformation::kWorldDirection, float_array);
                memcpy(directional_light_directions.data() + 3 * directional_light_count, &float_array[0], 3 * sizeof(float));

                // Color.
                light_params->getParameter(MLightParameterInformation::kColor, float_array);
                memcpy(directional_light_colors.data() + 3 * directional_light_count, &float_array[0], 3 * sizeof(float));

                // Intensity.
                light_params->getParameter(MLightParameterInformation::kIntensity, float_array);
                directional_light_intensities[directional_light_count] = float_array[0];

                ++directional_light_count;
            }
            else
            {
                std::stringstream ss;
                ss << "unsupported light type: " << light_type;
                LOG_ERROR(ss.str());
                continue;
            }
        }

        // Set shader params.

        CHECK_MSTATUS(shaderInstance->setParameter("point_light_count", point_light_count));
        CHECK_MSTATUS(shaderInstance->setArrayParameter("point_light_positions", point_light_positions.data(), point_light_count));
        CHECK_MSTATUS(shaderInstance->setArrayParameter("point_light_colors", point_light_colors.data(), point_light_count));
        CHECK_MSTATUS(shaderInstance->setArrayParameter("point_light_intensities", point_light_intensities.data(), point_light_count));

        CHECK_MSTATUS(shaderInstance->setParameter("directional_light_count", directional_light_count));
        CHECK_MSTATUS(shaderInstance->setArrayParameter("directional_light_directions", directional_light_directions.data(), directional_light_count));
        CHECK_MSTATUS(shaderInstance->setArrayParameter("directional_light_colors", directional_light_colors.data(), directional_light_count));
        CHECK_MSTATUS(shaderInstance->setArrayParameter("directional_light_intensities", directional_light_intensities.data(), directional_light_count));
    }

    SlicedDisplay::SlicedDisplay(MHWRender::MPxSubSceneOverride& parent) : m_parent(parent), m_enabled(false), m_selected(false), m_volume_sampler_state(nullptr)
    {
        const MHWRender::MShaderManager* shader_manager = getShaderManager();
        assert(shader_manager);

        // Load volume shader from effect file.
        m_volume_shader.reset(shader_manager->getEffectsBufferShader(s_effect_code.c_str(), unsigned(s_effect_code.size()), "Main", 0, 0, false, preDrawCallback));
        if (!m_volume_shader) {
            LOG_ERROR("Cannot compile cgfx.");
            return;
        }
        m_volume_shader->setIsTransparent(true);

        // Create sampler state for textures.
        MHWRender::MSamplerStateDesc volume_sampler_state_desc;
        volume_sampler_state_desc.filter = MHWRender::MSamplerState::kMinMagMipLinear;
        volume_sampler_state_desc.addressU = MHWRender::MSamplerState::kTexBorder;
        volume_sampler_state_desc.addressV = MHWRender::MSamplerState::kTexBorder;
        volume_sampler_state_desc.addressW = MHWRender::MSamplerState::kTexBorder;
        m_volume_sampler_state = MHWRender::MStateManager::acquireSamplerState(volume_sampler_state_desc);
        CHECK_MSTATUS(m_volume_shader->setParameter("scattering_sampler", *m_volume_sampler_state));
        CHECK_MSTATUS(m_volume_shader->setParameter("absorption_sampler", *m_volume_sampler_state));
        CHECK_MSTATUS(m_volume_shader->setParameter("emission_sampler",  *m_volume_sampler_state));
    }

    void SlicedDisplay::enable(bool enable)
    {
        m_enabled = enable;
        if (m_bbox_renderable.render_item)
            m_bbox_renderable.render_item->enable(m_enabled && m_selected);
        if (m_slices_renderable.render_item)
            m_slices_renderable.render_item->enable(m_enabled);
    }

    namespace {
        const char *SLICES_RENDER_ITEM_NAME = "vdb_volume_slices";
        const char *SELECTION_BBOX_RENDER_ITEM_NAME = "vdb_volume_slices_bbox";
    } // unnamed namespace

    bool SlicedDisplay::initRenderItems(MHWRender::MSubSceneContainer& container, const VDBSubSceneOverrideData& data)
    {
        if (!container.find(SLICES_RENDER_ITEM_NAME))
        {
            // Slices.

            auto render_item = MHWRender::MRenderItem::Create(
                    SLICES_RENDER_ITEM_NAME,
                    MHWRender::MRenderItem::RenderItemType::MaterialSceneItem,
                    MHWRender::MGeometry::kTriangles);
            render_item->setDrawMode(MHWRender::MGeometry::kAll);
            render_item->castsShadows(false);
            render_item->receivesShadows(false);
            if (!render_item->setShader(m_volume_shader.get())) {
                LOG_ERROR("Could not set shader for volume render item.");
                return false;
            }

            // Create geo buffers.
            // Note: descriptor name (first ctor arg) MUST be "", or setGeometryForRenderItem will return kFailure.
            const MHWRender::MVertexBufferDescriptor pos_desc("", MHWRender::MGeometry::kPosition, MHWRender::MGeometry::kFloat, 3);
            m_slices_renderable.position_buffer.reset(new MHWRender::MVertexBuffer(pos_desc));
            m_slices_renderable.index_buffer.reset(new MHWRender::MIndexBuffer(MHWRender::MGeometry::kUnsignedInt32));

            // Add render item to subscene container.
            if (!container.add(render_item)) {
                LOG_ERROR("Could not add m_slices_renderable render item.");
                return false;
            }
            m_slices_renderable.render_item = render_item;

            updateSliceGeo(data);
        }

        if (!container.find(SELECTION_BBOX_RENDER_ITEM_NAME))
        {
            // Selection bbox.

            const MHWRender::MShaderManager* shader_manager = getShaderManager();
            if (!shader_manager)
                return false;
            static auto shader = shader_manager->getStockShader(MHWRender::MShaderManager::k3dSolidShader);
            if (!shader) {
                LOG_ERROR("Couldn't get stock shader: k3dSolidShader.");
                return false;
            }

            auto render_item = MHWRender::MRenderItem::Create(
                    SELECTION_BBOX_RENDER_ITEM_NAME,
                    MHWRender::MRenderItem::RenderItemType::DecorationItem,
                    MHWRender::MGeometry::kLines);
            if (!render_item) {
                LOG_ERROR("Failed to create bbox render item.");
                return false;
            }
            render_item->setDrawMode(MHWRender::MGeometry::kAll);
            render_item->depthPriority(MHWRender::MRenderItem::sActiveWireDepthPriority);
            if (!render_item->setShader(shader)) {
                LOG_ERROR("Failed to set shader for bbox render item.");
                return false;
            }

            // Create geo buffers.
            // Note: descriptor name (first ctor arg) MUST be "", or setGeometryForRenderItem will return kFailure.
            const MHWRender::MVertexBufferDescriptor pos_desc("", MHWRender::MGeometry::kPosition, MHWRender::MGeometry::kFloat, 3);
            m_bbox_renderable.position_buffer.reset(new MHWRender::MVertexBuffer(pos_desc));
            m_bbox_renderable.vertex_buffer_array.clear();
            CHECK_MSTATUS(m_bbox_renderable.vertex_buffer_array.addBuffer("pos_model", m_bbox_renderable.position_buffer.get()));

            m_bbox_renderable.index_buffer.reset(new MHWRender::MIndexBuffer(MHWRender::MGeometry::kUnsignedInt32));
            constexpr auto index_count = 2 * 12;
            static const uint32_t BOX_WIREFRAME_INDICES[] = { 0, 1, 1, 3, 3, 2, 2, 0, 4, 5, 5, 7, 7, 6, 6, 4, 0, 4, 1, 5, 3, 7, 2, 6 };
            CHECK_MSTATUS(m_bbox_renderable.index_buffer->update(BOX_WIREFRAME_INDICES, 0, index_count, true));

            // Add render item to subscene container.
            if (!container.add(render_item)) {
                LOG_ERROR("Could not add bbox render item.");
                m_bbox_renderable.position_buffer.reset();
                m_bbox_renderable.index_buffer.reset();
                m_bbox_renderable.vertex_buffer_array.clear();
                return false;
            }
            m_bbox_renderable.render_item = render_item;

            updateBBox(data.bbox);
        }

        return true;
    }

    void SlicedDisplay::updateSliceGeo(const VDBSubSceneOverrideData& data)
    {
        // - Vertices
        const auto vertex_count = data.max_slice_count * 4;
        MFloatVector* positions = reinterpret_cast<MFloatVector*>(m_slices_renderable.position_buffer->acquire(vertex_count, true));
        for (int i = 0; i < data.max_slice_count; ++i) {
            const float z = i * 1.0f / (data.max_slice_count - 1.0f);
            positions[4 * i + 0] = MFloatVector(0.0, 0.0, z);
            positions[4 * i + 1] = MFloatVector(1.0, 0.0, z);
            positions[4 * i + 2] = MFloatVector(0.0, 1.0, z);
            positions[4 * i + 3] = MFloatVector(1.0, 1.0, z);
        }
        m_slices_renderable.position_buffer->commit(positions);

        // - Indices
        const auto index_count = data.max_slice_count * 6;
        unsigned int* indices = reinterpret_cast<unsigned int*>(m_slices_renderable.index_buffer->acquire(index_count, true));
        for (int i = 0; i < data.max_slice_count; ++i) {
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
        m_slices_renderable.update(m_parent, data.bbox);

        CHECK_MSTATUS(m_volume_shader->setParameter("max_slice_count", data.max_slice_count));
    }

    void SlicedDisplay::updateBBox(const MBoundingBox& bbox)
    {
        const auto extents = bbox.max() - bbox.min();
        constexpr auto vertex_count = 8;
        MFloatVector* positions = static_cast<MFloatVector*>(m_bbox_renderable.position_buffer->acquire(vertex_count, true));
        for (unsigned int i = 0; i < vertex_count; ++i) {
            positions[i] = bbox.min() + MVector((i & 1) * extents.x, ((i & 2) >> 1) * extents.y, ((i & 4) >> 2) * extents.z);
        }
        m_bbox_renderable.position_buffer->commit(positions);

        m_bbox_renderable.update(m_parent, bbox);

        CHECK_MSTATUS(m_volume_shader->setParameter("volume_origin", bbox.min()));
        CHECK_MSTATUS(m_volume_shader->setParameter("volume_size", extents));
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
        if (!m_volume_shader)
            return false;

        initRenderItems(container, data);
        if (!m_bbox_renderable || !m_slices_renderable)
            return false;

        // Handle selection.
        m_selected = data.is_selected;
        m_bbox_renderable.render_item->enable(m_enabled && m_selected);
        m_slices_renderable.render_item->enable(m_enabled);

        // Set wireframe color.
        const auto& color = data.wireframe_color;
        const float color_as_array[] = { color.r, color.g, color.b, color.a };
        CHECK_MSTATUS(m_bbox_renderable.render_item->getShader()->setParameter("solidColor", color_as_array));

        // Set world matrix.
        m_slices_renderable.render_item->setMatrix(&data.world_matrix);
        m_bbox_renderable.render_item->setMatrix(&data.world_matrix);

        // Update shader params.
        CHECK_MSTATUS(m_volume_shader->setParameter("scattering_color", data.scattering_color));
        CHECK_MSTATUS(m_volume_shader->setParameter("absorption_color", data.attenuation_color));
        CHECK_MSTATUS(m_volume_shader->setParameter("emission_color", data.emission_color));
        CHECK_MSTATUS(m_volume_shader->setParameter("shadow_gain", data.sliced_display_shader_params.shadow_gain));
        CHECK_MSTATUS(m_volume_shader->setParameter("shadow_sample_count", data.sliced_display_shader_params.shadow_sample_count));

        // === Update channels. ===

        typedef VDBSubSceneOverrideData::ChangeSet ChangeSet;

        // Bail if nothing has changed which would affect the volume textures.
        if (data.change_set <= ChangeSet::GENERIC_ATTRIBUTE)
            return true;

        const bool vdb_file_changed        = hasChange(data.change_set, ChangeSet::VDB_FILE);
        const bool max_slice_count_changed = hasChange(data.change_set, ChangeSet::MAX_SLICE_COUNT);
        const bool scattering_grid_changed = hasChange(data.change_set, ChangeSet::SCATTERING_GRID);
        const bool absorption_grid_changed = hasChange(data.change_set, ChangeSet::ATTENUATION_GRID);
        const bool emission_grid_changed   = hasChange(data.change_set, ChangeSet::EMISSION_GRID);

        // Needs to come before updateBBox: updateBBox updates the slices renderable too,
        // but at first run updateSliceGeo populates geo buffers.
        if (max_slice_count_changed)
            updateSliceGeo(data);

        // Update file-level bbox and invalidate channel cache if a new file is loaded.
        if (vdb_file_changed) {
            updateBBox(data.bbox);
            m_slices_renderable.update(m_parent, data.bbox);
            m_channel_cache.clear();
        }

        // Should come after channel cache is cleared.
        if (max_slice_count_changed)
            // Resample all channels.
            for (auto& channel : m_channel_cache)
                channel.second.lock()->sample(m_volume_sampler, data.max_slice_count);

        // Update channels.

        auto handle_channel_change = [&](const std::string& param_group, const std::string& channel_name, VolumeChannel::Ptr& channel) {
            if (channel_name.empty())
            {
                channel.reset();
                return;
            }

            auto it = m_channel_cache.find(channel_name);
            if (it == m_channel_cache.end())
            {
                // Create and cache channel.
                channel.reset(new VolumeChannel(data.vdb_file, channel_name));
                if (!channel->isValid())
                    return;
                channel->sample(m_volume_sampler, data.max_slice_count);
                m_channel_cache.insert({ channel_name, channel });
            }
            else
            {
                // Use cached channel.
                channel = it->second.lock();
            }

            const auto grid_bbox_is = getIndexSpaceBoundingBox(channel->grid.get());
            const auto grid_bbox_ws = channel->grid->transform().indexToWorld(grid_bbox_is);
            CHECK_MSTATUS(m_volume_shader->setParameter(format("^1s_volume_size", param_group), mayavecFromVec3f(grid_bbox_ws.extents())));
            CHECK_MSTATUS(m_volume_shader->setParameter(format("^1s_volume_origin", param_group), mayavecFromVec3f(grid_bbox_ws.min())));
        };

        if (vdb_file_changed || scattering_grid_changed)
            handle_channel_change("scattering", data.scattering_channel, m_scattering_channel);

        if (vdb_file_changed || absorption_grid_changed)
            handle_channel_change("absorption", data.attenuation_channel, m_absorption_channel);

        if (vdb_file_changed || emission_grid_changed)
            handle_channel_change("emission", data.emission_channel, m_emission_channel);

        // Clean-up channel cache.
        for (auto it = m_channel_cache.cbegin(); it != m_channel_cache.cend(); ++it)
            if (it->second.expired())
                it = m_channel_cache.erase(it);

        // Set channel params.
        MHWRender::MTextureAssignment texture_assignment;

        const bool has_scattering_texture = m_scattering_channel && m_scattering_channel->isValid();
        CHECK_MSTATUS(m_volume_shader->setParameter("use_scattering_texture", has_scattering_texture));
        if (has_scattering_texture)
        {
            texture_assignment.texture = m_scattering_channel->volume_texture.texture.get();
            CHECK_MSTATUS(m_volume_shader->setParameter("scattering_texture", texture_assignment));
            CHECK_MSTATUS(m_volume_shader->setParameter("scattering_value_range", mayavecFromFloatRange(m_scattering_channel->volume_texture.value_range)));
        }

        const bool has_absorption_texture = m_absorption_channel && m_absorption_channel->isValid();
        CHECK_MSTATUS(m_volume_shader->setParameter("use_absorption_texture", has_absorption_texture));
        if (has_absorption_texture)
        {
            texture_assignment.texture = m_absorption_channel->volume_texture.texture.get();
            CHECK_MSTATUS(m_volume_shader->setParameter("absorption_texture", texture_assignment));
            CHECK_MSTATUS(m_volume_shader->setParameter("absorption_value_range", mayavecFromFloatRange(m_absorption_channel->volume_texture.value_range)));
        }

        const bool has_emission_texture = m_emission_channel && m_emission_channel->isValid();
        CHECK_MSTATUS(m_volume_shader->setParameter("use_emission_texture", has_emission_texture));
        if (has_emission_texture)
        {
            texture_assignment.texture = m_emission_channel->volume_texture.texture.get();
            CHECK_MSTATUS(m_volume_shader->setParameter("emission_texture", texture_assignment));
            CHECK_MSTATUS(m_volume_shader->setParameter("emission_value_range", mayavecFromFloatRange(m_emission_channel->volume_texture.value_range)));
        }

        return true;
    }

    void Renderable::update(MHWRender::MPxSubSceneOverride& subscene_override, const MBoundingBox& bbox)
    {
        // Note: render item has to be added to the MSubSceneContainer before calling setGeometryForRenderItem.
        CHECK_MSTATUS(subscene_override.setGeometryForRenderItem(*render_item, vertex_buffer_array, *index_buffer, &bbox));
    }
}
