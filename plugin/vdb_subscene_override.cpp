#include "vdb_subscene_override.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <map>
#include <new>
#include <random>
#include <unordered_map>

#include <maya/MArgList.h>
#include <maya/MArgParser.h>
#include <maya/MDrawContext.h>
#include <maya/MGlobal.h>
#include <maya/MHwGeometryUtilities.h>
#include <maya/MSelectionList.h>
#include <maya/MShaderManager.h>

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>

#include <Cg/cg.h>

#include "paths.h"
#include "blackbody.h"
#include "progress_bar.h"
#include "volume_sampler.h"
#include "vdb_maya_utils.hpp"
#include "vdb_visualizer_data.h"

namespace {
    const MHWRender::MShaderManager* get_shader_manager()
    {
        auto renderer = MHWRender::MRenderer::theRenderer();
        if (renderer == nullptr)
            return nullptr;

        return renderer->getShaderManager();
    }

    MHWRender::MTextureManager* get_texture_manager()
    {
        auto renderer = MHWRender::MRenderer::theRenderer();
        if (renderer == nullptr)
            return nullptr;

        return renderer->getTextureManager();
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
    inline void log_error(const char *msg, const char *file_name, int line_no)
    {
#if _DEBUG
        std::cerr << "[openvdb_render] error: " << file_name << ": line " << line_no << ": " << msg << std::endl;
#else
        std::cerr << "[openvdb_render] error: " << msg << std::endl;
#endif
    }
    inline void log_error(const MString& msg, const char *file_name, int line_no)
    {
        log_error(msg.asChar(), file_name, line_no);
    }

} // unnamed namespace

// === Renderable ==========================================================

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

void Renderable::update(MHWRender::MPxSubSceneOverride& subscene_override, const MBoundingBox& bbox)
{
    // Note: render item has to be added to the MSubSceneContainer before calling setGeometryForRenderItem.
    CHECK_MSTATUS(subscene_override.setGeometryForRenderItem(*render_item, vertex_buffer_array, *index_buffer, &bbox));
}

// === VDBVolumeSpec =======================================================

struct VDBVolumeSpec
{
    std::string vdb_file_name;
    std::string vdb_file_tag; // To get different hash for different file versions.
    std::string vdb_grid_name;
    openvdb::Coord texture_size;

    VDBVolumeSpec() {}
    VDBVolumeSpec(const std::string& vdb_file_name_, const std::string& vdb_file_tag_, const std::string& vdb_grid_name_, openvdb::Coord texture_size_)
        : vdb_file_name(vdb_file_name_), vdb_file_tag(vdb_file_tag_), vdb_grid_name(vdb_grid_name_), texture_size(texture_size_) {}
};

namespace
{

    template<typename T> void hash_combine(size_t& seed, T const& v)
    {
        seed ^= std::hash<T>{}(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }

} // unnamed namespace

namespace std
{
    template <> struct hash<VDBVolumeSpec>
    {
        typedef VDBVolumeSpec argument_type;
        typedef std::size_t result_type;
        result_type operator()(argument_type const& spec) const
        {
            result_type res = std::hash<std::string>{}(spec.vdb_file_name);
            hash_combine(res, spec.vdb_file_tag);
            hash_combine(res, spec.vdb_grid_name);
            hash_combine(res, spec.texture_size.x());
            hash_combine(res, spec.texture_size.y());
            hash_combine(res, spec.texture_size.z());
            return res;
        }
    };

    bool operator==(const VDBVolumeSpec& lhs, const VDBVolumeSpec& rhs)
    {
        return lhs.vdb_file_name == rhs.vdb_file_name &&
               lhs.vdb_file_tag == rhs.vdb_file_tag &&
               lhs.vdb_grid_name == rhs.vdb_grid_name &&
               lhs.texture_size == rhs.texture_size;
    }
}

// === SamplerState ========================================================

class SamplerState
{
public:
    SamplerState(MHWRender::MSamplerState::TextureFilter filter,
        MHWRender::MSamplerState::TextureAddress address,
        int min_lod = 0, int max_lod = 16)
        : m_sampler_state(createSamplerState(filter, address, min_lod, max_lod))
    {
    }

    void assign(MHWRender::MShaderInstance *shader, const MString& param_name) const
    {
        CHECK_MSTATUS(shader->setParameter(param_name, *m_sampler_state));
    }

private:
    SamplerStatePtr m_sampler_state;

    static SamplerStatePtr createSamplerState(MHWRender::MSamplerState::TextureFilter filter,
        MHWRender::MSamplerState::TextureAddress address,
        int min_lod, int max_lod);
};

SamplerStatePtr SamplerState::createSamplerState(MHWRender::MSamplerState::TextureFilter filter,
    MHWRender::MSamplerState::TextureAddress address,
    int min_lod, int max_lod)
{
    MHWRender::MSamplerStateDesc desc;
    desc.filter = filter;
    desc.addressU = address;
    desc.addressV = address;
    desc.addressW = address;
    desc.minLOD = min_lod;
    desc.maxLOD = max_lod;
    desc.mipLODBias = 0;
    memset(desc.borderColor, 0, 4 * sizeof(float));
    return SamplerStatePtr(MHWRender::MStateManager::acquireSamplerState(desc));
}

// === VolumeTexture =======================================================

struct VolumeTexture
{
    TexturePtr texture_ptr;
    MFloatVector value_range;
    MFloatVector volume_size;
    MFloatVector volume_origin;
    openvdb::Coord extents;

    VolumeTexture() : texture_ptr(nullptr) {}
    VolumeTexture(const VolumeTexture&) = delete;
    VolumeTexture& operator=(const VolumeTexture&) = delete;
    VolumeTexture(VolumeTexture&&) = default;
    VolumeTexture& operator=(VolumeTexture&&) = default;

    template <typename RealType>
    void acquireBuffer(const openvdb::Coord& texture_extents, const VolumeBufferHandle<RealType>& volume_buffer);
    void clear() { texture_ptr.reset(); }
    bool isValid() const { return texture_ptr.get() != nullptr; }

    void assign(MHWRender::MShaderInstance* shader_instance, const MString& param)
    {
        m_texture_assignment.texture = texture_ptr.get();
        CHECK_MSTATUS(shader_instance->setParameter(param, m_texture_assignment));
    }

private:
    MHWRender::MTextureAssignment m_texture_assignment;
    static std::vector<float> s_staging;
};
std::vector<float> VolumeTexture::s_staging;

namespace {

    template <typename RealType>
    MFloatVector inline mayavecFromArray2(RealType *a)
    {
        return { a[0], a[1] };
    }

    template <typename RealType>
    MFloatVector inline mayavecFromArray3(RealType *a)
    {
        return { a[0], a[1], a[2] };
    }

    MFloatVector inline mayavecFromVec2f(const openvdb::Vec2f& vec)
    {
        return { vec.x(), vec.y() };
    }

    MFloatVector inline mayavecFromVec3f(const openvdb::Vec3f& vec)
    {
        return { vec.x(), vec.y(), vec.z() };
    }

    size_t voxel_count(const openvdb::Coord& extents)
    {
        return extents.x() * extents.y() * extents.z();
    }

} // unnamed namespace

template <typename RealType>
void VolumeTexture::acquireBuffer(const openvdb::Coord& texture_extents, const VolumeBufferHandle<RealType>& volume_buffer)
{
    if (!volume_buffer)
    {
        clear();
        return;
    }

    // Set metadata.
    value_range = mayavecFromArray2(volume_buffer.header->value_range);
    volume_size = mayavecFromArray3(volume_buffer.header->size);
    volume_origin = mayavecFromArray3(volume_buffer.header->origin);

    void* buffer = nullptr;
    if (std::is_same<RealType, float>::value)
    {
        buffer = volume_buffer.voxel_array;
    }
    else
    {
        // Convert voxels to float.
        // Note: uploading 'half' voxel data (i.e. using raster type kR16_FLOAT) is
        // SLOWER by a factor of ~2 compared to uploading floats. I have no idea why.
        // If you know the answer, please explain it to me (zoltan.gilian@gmail.com).
        s_staging.resize(voxel_count(texture_extents));
        typedef tbb::blocked_range<size_t> tbb_range;
        tbb::parallel_for(tbb_range(0, s_staging.size()),
            [input = reinterpret_cast<RealType*>(volume_buffer.voxel_array), output = s_staging.data()]
        (const tbb_range& range){
            for (size_t i = range.begin(); i < range.end(); ++i)
                output[i] = static_cast<float>(input[i]);
        });
        buffer = s_staging.data();
    }

    // If texture size didn't change, texture data can be updated in place,
    // providing there is an actual texture owned by this instance.
    if (extents == texture_extents && texture_ptr.get() != nullptr)
    {
        texture_ptr->update(buffer, true);
        return;
    }

    // Otherwise a new texture has to be created.
    MHWRender::MTextureDescription texture_desc;
    texture_desc.fWidth = texture_extents.x();
    texture_desc.fHeight = texture_extents.y();
    texture_desc.fDepth = texture_extents.z();
    texture_desc.fBytesPerRow = sizeof(float) * texture_desc.fWidth;
    texture_desc.fBytesPerSlice = texture_desc.fBytesPerRow * texture_desc.fHeight;
    texture_desc.fMipmaps = 0;
    texture_desc.fArraySlices = 1;
    texture_desc.fFormat = MHWRender::kR32_FLOAT;
    texture_desc.fTextureType = MHWRender::kVolumeTexture;
    texture_desc.fEnvMapType = MHWRender::kEnvNone;
    texture_ptr.reset(get_texture_manager()->acquireTexture("", texture_desc, buffer, true));

    extents = texture_extents;
}

// === VolumeCache =========================================================

class VolumeCache
{
public:
    static VolumeCache& instance();

    void getVolume(const VDBVolumeSpec& spec, VolumeTexture& output);

    enum class VoxelType { FLOAT, HALF };
    VoxelType getVoxelType() const { return m_voxel_type; }
    void setVoxelType(VoxelType voxel_type);

    void setMemoryLimitBytes(size_t mem_limit_bytes);
    size_t getMemoryLimitBytes() const { return m_mem_limit_bytes; }
    size_t getAllocatedBytes() const { return m_buffer.size(); }

private:
    static size_t s_refcount;

public:
    // Each SubSceneOverride should call registerUsage in its ctor and unregisterUsage in its dtor
    // so that the cache can be cleared if e.g. a new scene is created.
    static void registerUsage() { ++s_refcount; }
    static void unregisterUsage() { --s_refcount; if (s_refcount == 0) instance().clear(); }

private:
    VoxelType m_voxel_type;

    struct BufferRange
    {
        size_t begin;
        size_t end;
        BufferRange(size_t begin_, size_t end_) : begin(begin_), end(end_) {}
    };

    size_t m_mem_limit_bytes;

    // Objects are stored contiguously in a vector.
    std::vector<uint8_t> m_buffer;
    // FIFO cache eviction policy is used for its simplicity. The m_buffer_head
    // wraps around effectively creating a ring buffer, but objects are always
    // laid out linearly; m_buffer_head will wrap prematurely if the next object
    // doesn't fit into m_buffer.
    size_t m_buffer_head;
    // Associates VDBVolumeSpec values to ranges in m_buffer.
    // This makes the cache addressable by VDBVolumeSpec.
    typedef std::unordered_map<VDBVolumeSpec, BufferRange> BufferMap;
    BufferMap m_buffer_map;
    // m_allocation_map associates buffer range offsets ('begin') to buffer
    // map keys, so that old allocations overlapping a new one
    // can be deleted. The map needs to be ordered.
    std::map<size_t, VDBVolumeSpec> m_allocation_map;

    VolumeCache();
    openvdb::FloatGrid::ConstPtr loadGrid(const VDBVolumeSpec& spec);
    void clear();
    void clearRange(const BufferRange& range);
    void growBuffer(size_t minimum_buffer_size_bytes);

    template <typename RealType>
    void getVolume(const VDBVolumeSpec& spec, VolumeTexture& output);
    template <typename RealType>
    VolumeBufferHandle<RealType> allocate(const VDBVolumeSpec& spec);
    template <typename RealType>
    bool sampleGrid(const VDBVolumeSpec& spec, const openvdb::FloatGrid& grid, const VolumeBufferHandle<RealType>& output);

    static const size_t DEFAULT_LIMIT_BYTES;
    static const size_t DEFAULT_SIZE_BYTES;
    static const size_t GROW_AMOUNT_BYTES;
};

size_t VolumeCache::s_refcount = 0;

namespace {

    class VDBFile
    {
    public:
        VDBFile(const std::string& file_name) : m_vdb_file(file_name) { m_vdb_file.open(false); }
        ~VDBFile() { m_vdb_file.close(); }
        operator bool() const { return m_vdb_file.isOpen(); }
        openvdb::FloatGrid::ConstPtr loadFloatGrid(const std::string& grid_name);

    private:
        openvdb::io::File m_vdb_file;
    };

    openvdb::FloatGrid::ConstPtr VDBFile::loadFloatGrid(const std::string& grid_name)
    {
        if (!m_vdb_file.isOpen())
            return nullptr;

        openvdb::GridBase::ConstPtr grid_base_ptr;
        try
        {
            grid_base_ptr = m_vdb_file.readGrid(grid_name);
        }
        catch (const openvdb::Exception&)
        {
            return nullptr;
        }

        auto grid_ptr = openvdb::gridConstPtrCast<openvdb::FloatGrid>(grid_base_ptr);
        if (!grid_ptr)
        {
            LOG_ERROR(format("Grid '^1s' is not a FloatGrid.", grid_name));
            return nullptr;
        }

        return grid_ptr;
    }

    constexpr size_t KILOBYTE = 1024;
    constexpr size_t MEGABYTE = 1024 * KILOBYTE;
    constexpr size_t GIGABYTE = 1024 * MEGABYTE;

} // unnamed namespace

const size_t VolumeCache::DEFAULT_LIMIT_BYTES = 2 * GIGABYTE;
const size_t VolumeCache::DEFAULT_SIZE_BYTES = 256 * MEGABYTE;
const size_t VolumeCache::GROW_AMOUNT_BYTES = 256 * MEGABYTE;

VolumeCache& VolumeCache::instance()
{
    static VolumeCache volume_cache;
    return volume_cache;
}

void VolumeCache::setVoxelType(VoxelType voxel_type)
{
    if (m_voxel_type != voxel_type)
        clear();

    m_voxel_type = voxel_type;
}

VolumeCache::VolumeCache() : m_voxel_type(VoxelType::HALF), m_mem_limit_bytes(DEFAULT_LIMIT_BYTES), m_buffer_head(0)
{
    // Don't allocate anything in the ctor to avoid unnecessary consumption of memory (e.g. batch mode).
}

openvdb::FloatGrid::ConstPtr VolumeCache::loadGrid(const VDBVolumeSpec& spec)
{
    // Open VDB file or bail.
    auto vdb_file = VDBFile(spec.vdb_file_name);
    if (!vdb_file)
        return nullptr;

    return vdb_file.loadFloatGrid(spec.vdb_grid_name);
}

template <typename RealType>
bool VolumeCache::sampleGrid(const VDBVolumeSpec& spec, const openvdb::FloatGrid& grid, const VolumeBufferHandle<RealType>& output)
{
    // Set up volume sampler.
    ProgressBar progress_bar(format("vdb_visualizer: sampling grid ^1s", spec.vdb_grid_name));
    VolumeSampler<RealType> volume_sampler;
    volume_sampler.attachProgressBar(&progress_bar);

    // Sample grid.
    return volume_sampler.sampleGrid(grid, spec.texture_size, output);
}

void VolumeCache::getVolume(const VDBVolumeSpec& spec, VolumeTexture& output)
{
    if (m_voxel_type == VoxelType::HALF)
        getVolume<half>(spec, output);
    else if (m_voxel_type == VoxelType::FLOAT)
        getVolume<float>(spec, output);
}

template <typename RealType>
void VolumeCache::getVolume(const VDBVolumeSpec& spec, VolumeTexture& output)
{
    // Check if in cache.
    auto it = m_buffer_map.find(spec);
    if (it != m_buffer_map.end())
    {
        // Load from cache.
        output.acquireBuffer(spec.texture_size, VolumeBufferHandle<RealType>(m_buffer.data() + it->second.begin));
        return;
    }

    // Not in cache; try to load the grid.
    const auto grid = loadGrid(spec);
    if (!grid)
    {
        output.clear();
        return;
    }

    // Allocate space, return if not succesful.
    VolumeBufferHandle<RealType> volume_buffer = allocate<RealType>(spec);
    if (!volume_buffer)
    {
        output.clear();
        return;
    }

    // Sample grid, return if not succesful (i.e. user cancelled the sampling procedure).
    const auto status = sampleGrid<RealType>(spec, *grid, volume_buffer);
    if (!status)
    {
        output.clear();

        // Delete this allocation.

        const auto buffer_range_it = m_buffer_map.find(spec);
        if (buffer_range_it == m_buffer_map.end())
            // If caching is disabled, the allocation hasn't been recorded, so just return.
            return;

        m_allocation_map.erase(buffer_range_it->second.begin);
        m_buffer_map.erase(buffer_range_it);
        return;
    }

    // Update texture.
    output.acquireBuffer<RealType>(spec.texture_size, volume_buffer);

    // Clear buffer if caching is disabled.
    if (m_mem_limit_bytes == 0)
        m_buffer.clear();
}

void VolumeCache::setMemoryLimitBytes(size_t mem_limit_bytes)
{
    m_mem_limit_bytes = mem_limit_bytes;

    // Shrink buffer if requested.
    if (m_mem_limit_bytes < m_buffer.size())
    {
        clearRange(BufferRange(m_mem_limit_bytes, m_buffer.size()));
        m_buffer.resize(mem_limit_bytes);
        m_buffer.shrink_to_fit();
    }

    // Reset head if current positions will become invalid.
    if (m_buffer_head >= m_mem_limit_bytes)
        m_buffer_head = 0;
}

template <typename RealType>
VolumeBufferHandle<RealType> VolumeCache::allocate(const VDBVolumeSpec& spec)
{
    const size_t item_size_bytes = sizeof(VolumeBufferHeader<RealType>) + voxel_count(spec.texture_size) * sizeof(RealType);

    if (m_mem_limit_bytes == 0)
    {
        // Caching is disabled. Use buffer for this request, but perform no accounting.
        m_buffer.resize(item_size_bytes);
        m_buffer.shrink_to_fit();
        return { m_buffer.data() };
    }
    else if (item_size_bytes > m_mem_limit_bytes)
        return {};

    // If this allocation would exceed the memory limit, reset the head index.
    // Otherwise grow the buffer if needed.
    {
        const auto allocation_end = m_buffer_head + item_size_bytes;
        if (m_buffer.size() < allocation_end)
            growBuffer(allocation_end);
    }

    // Allocate buffer range.
    const size_t buffer_begin = m_buffer_head;
    const size_t buffer_end = buffer_begin + item_size_bytes;
    const auto buffer_range = BufferRange(buffer_begin, buffer_end);
    m_buffer_head += item_size_bytes;

    // Throw away old allocations which overlap [buffer_begin, buffer_end).
    clearRange(buffer_range);

    // Update maps.
    m_buffer_map.insert(std::make_pair(spec, buffer_range)).first;
    m_allocation_map.insert(std::make_pair(buffer_begin, spec));
    return { m_buffer.data() + buffer_begin };
}

void VolumeCache::growBuffer(size_t minimum_buffer_size_bytes)
{
    // Reset head if requested size exceeds the limit.
    if (minimum_buffer_size_bytes > m_mem_limit_bytes)
    {
        m_buffer_head = 0;
        return;
    }

    size_t new_size_bytes = m_buffer.size() + GROW_AMOUNT_BYTES;
    new_size_bytes = std::max(minimum_buffer_size_bytes, new_size_bytes);
    new_size_bytes = std::min(m_mem_limit_bytes, new_size_bytes);

    // Bail if buffer size is large enough already.
    if (new_size_bytes <= m_buffer.size())
        return;

    // Try to resize or reset head if failed.
    try
    {
        m_buffer.resize(new_size_bytes);
    }
    catch (const std::bad_alloc&)
    {
        m_buffer_head = 0;
    }
}

void VolumeCache::clear()
{
    m_buffer_head = 0;
    m_buffer.clear();
    m_buffer.shrink_to_fit();
    m_buffer_map.clear();
    m_allocation_map.clear();
}

void VolumeCache::clearRange(const BufferRange& range_to_clear)
{
    // Find first range which begins strictly after the given range.begin.
    auto erase_it = m_allocation_map.upper_bound(range_to_clear.begin);
    // Start erasing at the previous range if it contains range_to_clear.begin.
    if (erase_it != m_allocation_map.begin())
    {
        auto prev_it = erase_it;
        --prev_it;
        const auto prev_range_it = m_buffer_map.find(prev_it->second);
        assert(prev_range_it != m_buffer_map.end());
        const auto& prev_range = prev_range_it->second;
        if (prev_range.begin <= range_to_clear.begin && range_to_clear.begin < prev_range.end)
            erase_it = prev_it;
    }
    // The first range to keep is the first one which begins on or after range.end.
    auto erase_end = m_allocation_map.lower_bound(range_to_clear.end);

    while (erase_it != erase_end)
    {
        m_buffer_map.erase(erase_it->second);
        erase_it = m_allocation_map.erase(erase_it);
    }
}

// === VolumeParam =========================================================

class VolumeParam
{
public:
    VolumeParam(const char *shader_param_prefix = nullptr, MHWRender::MShaderInstance *shader_instance = nullptr);
    void setParamPrefix(const char *param_prefix);
    void setShaderInstance(MHWRender::MShaderInstance *shader_instance) { m_shader_instance = shader_instance; }
    void setVolume(const VDBVolumeSpec& volume_spec);

private:
    MString use_texture_param;
    MString texture_param;
    MString value_range_param;
    MString volume_size_param;
    MString volume_origin_param;
    MHWRender::MShaderInstance *m_shader_instance;
    VolumeTexture m_volume_texture;

    void assign();
};

VolumeParam::VolumeParam(const char *shader_param_prefix, MHWRender::MShaderInstance *shader_instance)
{
    setParamPrefix(shader_param_prefix);
    setShaderInstance(shader_instance);
}

void VolumeParam::setParamPrefix(const char *param_prefix)
{
    if (!param_prefix)
        return;
    use_texture_param = format("use_^1s_texture", param_prefix);
    texture_param = format("^1s_texture", param_prefix);
    value_range_param = format("^1s_value_range", param_prefix);
    volume_size_param = format("^1s_volume_size", param_prefix);
    volume_origin_param = format("^1s_volume_origin", param_prefix);
}

void VolumeParam::setVolume(const VDBVolumeSpec& volume_spec)
{
    VolumeCache::instance().getVolume(volume_spec, m_volume_texture);
    assign();
}

void VolumeParam::assign()
{
    const bool use_texture = m_volume_texture.isValid();
    CHECK_MSTATUS(m_shader_instance->setParameter(use_texture_param, use_texture));
    if (!use_texture)
        return;

    m_volume_texture.assign(m_shader_instance, texture_param);
    CHECK_MSTATUS(m_shader_instance->setParameter(value_range_param, m_volume_texture.value_range));
    CHECK_MSTATUS(m_shader_instance->setParameter(volume_size_param, m_volume_texture.volume_size));
    CHECK_MSTATUS(m_shader_instance->setParameter(volume_origin_param, m_volume_texture.volume_origin));
}

// === RampTextureBase =====================================================

class RampTextureBase
{
public:
    void assignSamplerToShader(MHWRender::MShaderInstance* shader_instance, const MString& sampler_param);
    void assignTextureToShader(MHWRender::MShaderInstance* shader_instance, const MString& texture_param) const;
    void assignDomainToShader(MHWRender::MShaderInstance* shader_instance, const MString& domain_param) const;

    operator bool() const { return m_texture.get() != nullptr; }

protected:
    RampTextureBase(const unsigned int resolution, const MHWRender::MRasterFormat raster_format, const unsigned int bytes_per_pixel);

    unsigned int m_resolution;
    std::vector<uint8_t> m_staging;

    TexturePtr m_texture;
    openvdb::Vec2f m_domain;
    const SamplerState m_ramp_sampler_state;
};

RampTextureBase::RampTextureBase(const unsigned int resolution, const MHWRender::MRasterFormat raster_format, const unsigned int bytes_per_pixel)
    : m_resolution(resolution), m_staging(bytes_per_pixel * resolution, 0),
    m_ramp_sampler_state(MHWRender::MSamplerState::kMinMagMipLinear, MHWRender::MSamplerState::kTexClamp, 0, 0)
{
    MHWRender::MTextureDescription ramp_desc;
    ramp_desc.fWidth = resolution;
    ramp_desc.fHeight = 1;
    ramp_desc.fDepth = 1;
    ramp_desc.fBytesPerRow = bytes_per_pixel * ramp_desc.fWidth;
    ramp_desc.fBytesPerSlice = ramp_desc.fBytesPerRow;
    ramp_desc.fMipmaps = 1;
    ramp_desc.fArraySlices = 1;
    ramp_desc.fFormat = raster_format;
    ramp_desc.fTextureType = MHWRender::kImage1D;
    ramp_desc.fEnvMapType = MHWRender::kEnvNone;
    m_texture.reset(get_texture_manager()->acquireTexture("", ramp_desc, m_staging.data(), false));
}

void RampTextureBase::assignSamplerToShader(MHWRender::MShaderInstance* shader_instance, const MString& sampler_param)
{
    m_ramp_sampler_state.assign(shader_instance, sampler_param);
}

void RampTextureBase::assignTextureToShader(MHWRender::MShaderInstance* shader_instance, const MString& texture_param) const
{
    MHWRender::MTextureAssignment ta;
    ta.texture = m_texture.get();
    CHECK_MSTATUS(shader_instance->setParameter(texture_param, ta));
}

void RampTextureBase::assignDomainToShader(MHWRender::MShaderInstance* shader_instance, const MString& domain_param) const
{
    CHECK_MSTATUS(shader_instance->setParameter(domain_param, mayavecFromVec2f(m_domain)));
}

// === FloatRampTexture ====================================================

class FloatRampTexture : public RampTextureBase
{
public:
    FloatRampTexture(const unsigned int resolution) : RampTextureBase(resolution, MHWRender::MRasterFormat::kR8_UNORM, 1) {}

    void updateFromData(const float *data);
    void updateFromGradient(const Gradient& gradient);

private:
};

void FloatRampTexture::updateFromData(const float *floats)
{
    for (unsigned int i = 0; i < m_resolution; ++i)
        m_staging[i] = uint8_t(floats[i] * 255);
    m_texture->update(m_staging.data(), false);
}

void FloatRampTexture::updateFromGradient(const Gradient& gradient)
{
    if (!m_texture)
        return;

    const auto& sample_vector = gradient.getFloatRamp();
    if (sample_vector.size() < m_resolution)
        return;

    updateFromData(sample_vector.data());
    m_domain = openvdb::Vec2d(gradient.getInputMin(), gradient.getInputMax());
}

// === RGBRampTexture ======================================================

class RGBRampTexture : public RampTextureBase
{
public:
    RGBRampTexture(const unsigned int resolution) : RampTextureBase(resolution, MHWRender::MRasterFormat::kR8G8B8X8, 4) {}

    void updateFromData(const MFloatVector *colors, const float normalizer=1);
    void updateFromGradient(const Gradient& gradient);

private:
};

void RGBRampTexture::updateFromData(const MFloatVector *colors, const float normalizer)
{
    for (unsigned int i = 0; i < m_resolution; ++i)
    {
        auto srgb_color = SRGBFromLinear(colors[i] / normalizer);
        m_staging[4 * i + 0] = uint8_t(srgb_color.x * 255);
        m_staging[4 * i + 1] = uint8_t(srgb_color.y * 255);
        m_staging[4 * i + 2] = uint8_t(srgb_color.z * 255);
        m_staging[4 * i + 3] = 0;
    }
    m_texture->update(m_staging.data(), false);
}

void RGBRampTexture::updateFromGradient(const Gradient& gradient)
{
    if (!m_texture)
        return;

    const auto& sample_vector = gradient.getRgbRamp();
    if (sample_vector.size() < m_resolution)
        return;

    updateFromData(sample_vector.data());
    m_domain = openvdb::Vec2d(gradient.getInputMin(), gradient.getInputMax());
}

// === BlackbodyLUT ========================================================

struct BlackbodyLUT
{
    BlackbodyLUT();
    RGBRampTexture lut;
};

BlackbodyLUT::BlackbodyLUT() : lut(Blackbody::TABLE_SIZE)
{
    lut.updateFromData(Blackbody::LUT, Blackbody::LUT_NORMALIZER);
}

// === VolumeShader ========================================================

class VolumeShader
{
public:
    VolumeShader();
    ~VolumeShader();
    MHWRender::MShaderInstance *get() { return m_shader_ptr; }
    const MHWRender::MShaderInstance *get() const { return m_shader_ptr; }
    MHWRender::MShaderInstance *operator*() { return get(); }
    const MHWRender::MShaderInstance *operator*() const { return get(); }
    MHWRender::MShaderInstance *operator->() { return get(); }
    const MHWRender::MShaderInstance *operator->() const { return get(); }

    operator bool() const { return get() != nullptr; }

private:
    ShaderPtr m_shader_clone_owner;
    MHWRender::MShaderInstance *m_shader_ptr;

    static void loadShader();
    static void preDrawCallback(MHWRender::MDrawContext& context,
                                const MHWRender::MRenderItemList& renderItemList,
                                MHWRender::MShaderInstance* shaderInstance);

    static ShaderPtr s_volume_shader;
    static size_t s_refcount;
    static const unsigned int MAX_LIGHT_COUNT;
    static const std::string VOLUME_EFFECT_CODE;
};

ShaderPtr VolumeShader::s_volume_shader = nullptr;
size_t VolumeShader::s_refcount = 0;

const unsigned int VolumeShader::MAX_LIGHT_COUNT = 16;

namespace {

    template <typename T>
    MHWRender::MShaderCompileMacro makeMacroDef(const MString& name, const T& value)
    {
        return { name, format("^1s", value) };
    };

    template <>
    MHWRender::MShaderCompileMacro makeMacroDef<float>(const MString& name, const float& value)
    {
        return { name, format("^1sf", value) };
    }

} // unnamed namespace

VolumeShader::VolumeShader()
{
    if (s_refcount == 0)
    {
        // Load the shader. 
        loadShader();
        // The first instance uses the static MShaderInstance.
        m_shader_ptr = s_volume_shader.get();
    }
    else if (s_volume_shader.get() != nullptr)
    {
        // Not the first instance, clone the static shader.
        m_shader_clone_owner.reset(s_volume_shader->clone());
        m_shader_ptr = m_shader_clone_owner.get();
    }
    ++s_refcount;
}

VolumeShader::~VolumeShader()
{
    --s_refcount;
    if (s_refcount == 0)
        // Last instance; release the static MShaderInstance.
        s_volume_shader.reset();
}

void VolumeShader::loadShader()
{
    // Only OpenGL core profile is supported.
    if (MHWRender::MRenderer::theRenderer()->drawAPI() != MHWRender::kOpenGLCoreProfile) {
        MGlobal::displayError("[openvdb_render] the visualizer node only works "
                              "with OpenGL Core Profile rendering engine.");
        return;
    }

    const MHWRender::MShaderManager* shader_manager = get_shader_manager();
    if (!shader_manager)
        return;

    MHWRender::MShaderCompileMacro macros[] = { makeMacroDef("BLACKBODY_LUT_MIN_TEMP",   Blackbody::TEMPERATURE_MIN),
                                                makeMacroDef("BLACKBODY_LUT_MAX_TEMP",   Blackbody::TEMPERATURE_MAX),
                                                makeMacroDef("BLACKBODY_LUT_NORMALIZER", Blackbody::LUT_NORMALIZER),
                                                makeMacroDef("MAX_LIGHT_COUNT",          MAX_LIGHT_COUNT) };
    constexpr int macro_count = sizeof(macros) / sizeof(MHWRender::MShaderCompileMacro);
    const auto effect_file = Paths::getVolumeEffectFile();
    s_volume_shader.reset(shader_manager->getEffectsFileShader(
            effect_file.c_str(),
            "Main", macros, macro_count, /*useEffectCache=*/false, preDrawCallback));
    if (!s_volume_shader) {
        return;
    }

    s_volume_shader->setIsTransparent(true);
}

namespace {

    typedef std::array<float, 3> Float3;

    template <typename ParamSpec>
    void getLightParam(MHWRender::MLightParameterInformation* light_params, const ParamSpec& param_spec, float& output) {
        MFloatArray float_array;
        if (light_params->getParameter(param_spec, float_array) == MStatus::kSuccess)
            output = float_array[0];
    };

    template <typename ParamSpec>
    void getLightParam(MHWRender::MLightParameterInformation* light_params, const ParamSpec& param_spec, Float3& output) {
        MFloatArray float_array;
        if (light_params->getParameter(param_spec, float_array) == MStatus::kSuccess)
            memcpy(&output, &float_array[0], 3 * sizeof(float));
    };

    float radiansFromDegrees(float degrees)
    {
        return float(M_PI) * degrees / 180.0f;
    }

    template <int N>
    union Float3Array
    {
        std::array<float, 3 * N> float_array;
        std::array<Float3, N> float3_array;
    };

} // unnamed namespace

void VolumeShader::preDrawCallback(MHWRender::MDrawContext& context, const MHWRender::MRenderItemList& /*renderItemList*/, MHWRender::MShaderInstance* shader_instance)
{
    // Check for errors.
    if (shader_instance->bind(context) != MStatus::kSuccess) {
        const auto shader_manager = get_shader_manager();
        std::cerr << "[openvdb_render] Failed to compile volume shader: "
            << shader_manager->getLastError() << std::endl;
        return;
    }

    // Set blend state.
    {
        auto state_manager = context.getStateManager();
        MHWRender::MBlendStateDesc blend_state_desc;
        blend_state_desc.setDefaults();
        for (int i = 0; i < (blend_state_desc.independentBlendEnable ? MHWRender::MBlendState::kMaxTargets : 1); ++i) {
            blend_state_desc.targetBlends[i].blendEnable = true;
            blend_state_desc.targetBlends[i].sourceBlend = MHWRender::MBlendState::kOne;
            blend_state_desc.targetBlends[i].destinationBlend = MHWRender::MBlendState::kInvSourceAlpha;
            blend_state_desc.targetBlends[i].blendOperation = MHWRender::MBlendState::kAdd;
            blend_state_desc.targetBlends[i].alphaSourceBlend = MHWRender::MBlendState::kOne;
            blend_state_desc.targetBlends[i].alphaDestinationBlend = MHWRender::MBlendState::kInvSourceAlpha;
            blend_state_desc.targetBlends[i].alphaBlendOperation = MHWRender::MBlendState::kAdd;
        }

        blend_state_desc.blendFactor[0] = 1.0f;
        blend_state_desc.blendFactor[1] = 1.0f;
        blend_state_desc.blendFactor[2] = 1.0f;
        blend_state_desc.blendFactor[3] = 1.0f;

        const MHWRender::MBlendState* blend_state = state_manager->acquireBlendState(blend_state_desc);
        CHECK_MSTATUS(state_manager->setBlendState(blend_state));
    }

    // Set view position.
    {
        MStatus status;
        const auto view_pos = context.getTuple(MHWRender::MFrameContext::kViewPosition, &status);
        CHECK_MSTATUS(status);
        CHECK_MSTATUS(shader_instance->setParameter("view_pos_world", MFloatVector(float(view_pos[0]), float(view_pos[1]), float(view_pos[2]))));
    }

    // Set view direction in model space and dominant axis index (0: x, 1: y, 2: z).
    {
        MStatus status;
        const auto view_dir = context.getTuple(MHWRender::MFrameContext::kViewDirection, &status);
        CHECK_MSTATUS(status);
        const auto world_inverse_mat = context.getMatrix(MHWRender::MFrameContext::kWorldInverseMtx, &status);
        CHECK_MSTATUS(status);
        auto view_dir_vec = MVector(view_dir[0], view_dir[1], view_dir[2]);
        view_dir_vec *= world_inverse_mat;
        view_dir_vec.normalize();
        const auto view_dir_model = MFloatVector(view_dir_vec[0], view_dir_vec[1], view_dir_vec[2]);
        CHECK_MSTATUS(shader_instance->setParameter("view_dir_model", view_dir_model));

        const float ax = std::fabs(view_dir_model.x), ay = std::fabs(view_dir_model.y), az = std::fabs(view_dir_model.z);
        int dominant_axis = 0;
        if (ay > ax && ay > az)
            dominant_axis = 1;
        else if (az > ax && az > ay)
            dominant_axis = 2;
        CHECK_MSTATUS(shader_instance->setParameter("dominant_axis", dominant_axis));
    }

    // Collect light data.

    constexpr int LIGHT_FLAG_POINT_LIGHT       = 0;
    constexpr int LIGHT_FLAG_DIRECTIONAL_LIGHT = 1;
    constexpr int LIGHT_FLAG_SPOTLIGHT         = 2;
    constexpr int LIGHT_FLAG_CAST_SHADOWS      = 8;

    Float3Array<MAX_LIGHT_COUNT> light_position;
    Float3Array<MAX_LIGHT_COUNT> light_direction;
    Float3Array<MAX_LIGHT_COUNT> light_color;
    Float3Array<MAX_LIGHT_COUNT> light_shadow_color;
    std::array<int,   MAX_LIGHT_COUNT> light_flags;
    std::array<float, MAX_LIGHT_COUNT> light_intensity;
    std::array<float, MAX_LIGHT_COUNT> light_decay_exponent;
    std::array<float, MAX_LIGHT_COUNT> light_cutoff_costheta1;
    std::array<float, MAX_LIGHT_COUNT> light_cutoff_costheta2;
    std::array<float, MAX_LIGHT_COUNT> light_dropoff;

    using MHWRender::MLightParameterInformation;

    const auto light_count_total = std::min(context.numberOfActiveLights(), MAX_LIGHT_COUNT);
    int shader_light_count = 0;
    for (unsigned int i = 0; i < light_count_total; ++i)
    {
        MIntArray int_array;
        MFloatArray float_array;

        const auto light_params = context.getLightParameterInformation(i);

        // Proceed to next light if this light is not enabled.
        {
            const auto status = light_params->getParameter(MLightParameterInformation::kLightEnabled, float_array);
            if (status != MS::kSuccess || float_array[0] != 1)
                continue;
        }

        light_flags[shader_light_count] = 0;

        // Light type.
        const auto light_type = light_params->lightType();
        if (light_type == "pointLight")
            light_flags[shader_light_count] |= LIGHT_FLAG_POINT_LIGHT;
        else if (light_type == "directionalLight")
            light_flags[shader_light_count] |= LIGHT_FLAG_DIRECTIONAL_LIGHT;
        else if (light_type == "spotLight")
            light_flags[shader_light_count] |= LIGHT_FLAG_SPOTLIGHT;
        else
        {
            LOG_ERROR(format("unsupported light type ^1s", light_type).asChar());
            continue;
        }

        // Position.
        getLightParam(light_params, MLightParameterInformation::kWorldPosition, light_position.float3_array[shader_light_count]);

        // Direction.
        getLightParam(light_params, MLightParameterInformation::kWorldDirection, light_direction.float3_array[shader_light_count]);

        // Color.
        getLightParam(light_params, MLightParameterInformation::kColor, light_color.float3_array[shader_light_count]);

        // Intensity.
        getLightParam(light_params, MLightParameterInformation::kIntensity, light_intensity[shader_light_count]);

        // Shadow.
        {
            // Cast shadows if the ShadowOn param is on, or true by default if there is no such param.
            auto status = light_params->getParameter(MLightParameterInformation::kShadowOn, int_array);
            if (status == MStatus::kSuccess && int_array[0] == 1) {
                light_flags[shader_light_count] |= LIGHT_FLAG_CAST_SHADOWS;
            }
#if MAYA_API_VERSION < 201600
            // ShadowOn attrib of point lights is buggy in Maya 2015, cast shadows by default.
            if (light_type == "pointLight") {
                light_flags[shader_light_count] |= LIGHT_FLAG_CAST_SHADOWS;
            }
#endif
        }

        // Shadow color.
        getLightParam(light_params, MLightParameterInformation::kShadowColor, light_shadow_color.float3_array[shader_light_count]);

        // Decay rate.
        getLightParam(light_params, MLightParameterInformation::kDecayRate, light_decay_exponent[shader_light_count]);

        // Spot light dropoff.
        getLightParam(light_params, MLightParameterInformation::kDropoff, light_dropoff[shader_light_count]);

        // Spot light cutoffs.
        if (light_params->getParameter(MLightParameterInformation::kCosConeAngle, float_array) == MStatus::kSuccess)
        {
            light_cutoff_costheta1[shader_light_count] = float_array[0];
            light_cutoff_costheta2[shader_light_count] = float_array[1];
        }

        shader_light_count++;
    }

    // Convert colors to linear color space.
    LinearFromSRGB(light_color.float_array.data(), 3 * shader_light_count);
    LinearFromSRGB(light_shadow_color.float_array.data(), 3 * shader_light_count);

    // Set shader params.
    CHECK_MSTATUS(shader_instance->setParameter("light_count", shader_light_count));
    CHECK_MSTATUS(shader_instance->setArrayParameter("light_flags", light_flags.data(), shader_light_count));
    CHECK_MSTATUS(shader_instance->setArrayParameter("light_position", light_position.float_array.data(), shader_light_count));
    CHECK_MSTATUS(shader_instance->setArrayParameter("light_direction", light_direction.float_array.data(), shader_light_count));
    CHECK_MSTATUS(shader_instance->setArrayParameter("light_color", light_color.float_array.data(), shader_light_count));
    CHECK_MSTATUS(shader_instance->setArrayParameter("light_intensity", light_intensity.data(), shader_light_count));
    CHECK_MSTATUS(shader_instance->setArrayParameter("light_shadow_color", light_shadow_color.float_array.data(), shader_light_count));
    CHECK_MSTATUS(shader_instance->setArrayParameter("light_decay_exponent", light_decay_exponent.data(), shader_light_count));
    CHECK_MSTATUS(shader_instance->setArrayParameter("light_cutoff_costheta1", light_cutoff_costheta1.data(), shader_light_count));
    CHECK_MSTATUS(shader_instance->setArrayParameter("light_cutoff_costheta2", light_cutoff_costheta2.data(), shader_light_count));
    CHECK_MSTATUS(shader_instance->setArrayParameter("light_dropoff", light_dropoff.data(), shader_light_count));
}


// === SlicedDisplay =======================================================

class SlicedDisplay
{
public:
    SlicedDisplay(MHWRender::MPxSubSceneOverride& parent);
    ~SlicedDisplay();
    bool update(MHWRender::MSubSceneContainer& container, const VDBSubSceneOverrideData& data);
    void enable(bool enable);

private:
    bool initRenderItems(MHWRender::MSubSceneContainer& container, const VDBSubSceneOverrideData& data);

    void updateBBox(const MBoundingBox& bbox);
    void updateSliceGeo(const VDBSubSceneOverrideData& data);

    MHWRender::MPxSubSceneOverride& m_parent;

    VolumeShader m_volume_shader;

    VolumeParam m_density_channel;
    VolumeParam m_scattering_channel;
    VolumeParam m_emission_channel;
    VolumeParam m_transparency_channel;
    VolumeParam m_temperature_channel;

    // Must be the same as Gradient resolution.
    static const unsigned int RAMP_RESOLUTION;
    FloatRampTexture m_density_ramp;
    RGBRampTexture m_scattering_ramp;
    RGBRampTexture m_emission_ramp;

    BlackbodyLUT m_blackbody_lut;

    Renderable m_slices_renderable;
    Renderable m_bbox_renderable;
    SamplerState m_volume_sampler_state;

    bool m_enabled;
    bool m_selected;
    bool m_force_shadows_off;
};
const unsigned int SlicedDisplay::RAMP_RESOLUTION = 128;

// === VDBSubSceneOverrideData ====================================================

enum class ChangeSet : unsigned int {
    NO_CHANGES = 0,
    GENERIC_ATTRIBUTE =    1<<0,
    GRADIENT =             1<<1,
    VDB_FILE =             1<<2,
    MAX_SLICE_COUNT =      1<<3,
    SCATTERING_CHANNEL =   1<<4,
    ATTENUATION_CHANNEL =  1<<5,
    EMISSION_CHANNEL =     1<<6,
    DENSITY_CHANNEL =      1<<7,
    TRANSPARENCY_CHANNEL = 1<<8,
    TEMPERATURE_CHANNEL =  1<<9,
};

struct VDBSubSceneOverrideData : public VDBVisualizerData {
    MMatrix world_matrix;
    MColor wireframe_color;
    bool is_selected, is_visible;

    openvdb::GridBase::ConstPtr scattering_grid;
    openvdb::GridBase::ConstPtr attenuation_grid;
    openvdb::GridBase::ConstPtr emission_grid;

    ChangeSet change_set;

    VDBSubSceneOverrideData();
    ~VDBSubSceneOverrideData();
    void clear();
    bool update(const VDBVisualizerData* data, const MObject& obj);
};

inline ChangeSet& operator|=(ChangeSet& lhs, ChangeSet rhs)
{
    return lhs = ChangeSet(unsigned(lhs) | unsigned(rhs));
}
inline ChangeSet operator|(ChangeSet lhs, ChangeSet rhs)
{
    lhs |= rhs;
    return lhs;
}
inline ChangeSet& operator&=(ChangeSet& lhs, ChangeSet rhs)
{
    return lhs = ChangeSet(unsigned(lhs) & unsigned(rhs));
}
inline ChangeSet operator&(ChangeSet lhs, ChangeSet rhs)
{
    lhs &= rhs;
    return lhs;
}
inline bool hasChange(ChangeSet change_set, ChangeSet mask)
{
    return (change_set & mask) != ChangeSet::NO_CHANGES;
}
inline ChangeSet operator~(ChangeSet cs)
{
    return ChangeSet(~unsigned(cs));
}

VDBSubSceneOverrideData::VDBSubSceneOverrideData() : is_selected(false), change_set(ChangeSet::NO_CHANGES)
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

namespace {

    ChangeSet setup_channel(ChannelParams& target, const ChannelParams& source, ChangeSet channel_mask)
    {
        if (target.name != source.name)
        {
            target = source;
            return channel_mask;
        }
        else if (target.gradient != source.gradient)
        {
            target = source;
            return ChangeSet::GRADIENT;
        }
        else if (target.color_source != source.color_source || target.intensity != source.intensity || target.color != source.color)
        {
            target = source;
            return ChangeSet::GENERIC_ATTRIBUTE;
        }
        else
        {
            return ChangeSet::NO_CHANGES;
        }
    }

} // unnamed namespace

bool VDBSubSceneOverrideData::update(const VDBVisualizerData* data, const MObject& obj)
{
    // TODO: we can limit some of the comparisons to the display mode
    // ie, we don't need to compare certain things if we are using the bounding
    // box mode
    change_set = ChangeSet::NO_CHANGES;

    MDagPath dg = MDagPath::getAPathTo(obj);
    const MMatrix inc_world_matrix = dg.inclusiveMatrix();
    change_set |= setup_parameter(world_matrix, inc_world_matrix, ChangeSet::GENERIC_ATTRIBUTE);
    change_set |= setup_parameter(wireframe_color, MHWRender::MGeometryUtilities::wireframeColor(dg), ChangeSet::GENERIC_ATTRIBUTE);
    change_set |= setup_parameter(is_selected, isPathSelected(dg), ChangeSet::GENERIC_ATTRIBUTE);

    if (dg.isVisible() && !is_visible)
        // Update everything if becoming visible.
        change_set = ~ChangeSet::NO_CHANGES;
    is_visible = dg.isVisible();
    //change_set |= setup_parameter(is_visible, dg.isVisible(), ChangeSet::GENERIC_ATTRIBUTE);

    if (data == nullptr || update_trigger == data->update_trigger)
        return change_set != ChangeSet::NO_CHANGES;

    update_trigger = data->update_trigger;

    // A change in display mode triggers complete recalculation; mark everything as changed.
    change_set |= setup_parameter(display_mode, data->display_mode, ~ChangeSet::NO_CHANGES);

    change_set |= setup_parameter(bbox, data->bbox, ChangeSet::GENERIC_ATTRIBUTE);

    change_set |= setup_channel(density_channel, data->density_channel, ChangeSet::DENSITY_CHANNEL);
    change_set |= setup_channel(scattering_channel, data->scattering_channel, ChangeSet::SCATTERING_CHANNEL);
    change_set |= setup_channel(attenuation_channel, data->attenuation_channel, ChangeSet::ATTENUATION_CHANNEL);
    change_set |= setup_channel(emission_channel, data->emission_channel, ChangeSet::EMISSION_CHANNEL);
    change_set |= setup_channel(transparency_channel, data->transparency_channel, ChangeSet::TRANSPARENCY_CHANNEL);
    change_set |= setup_channel(temperature_channel, data->temperature_channel, ChangeSet::TEMPERATURE_CHANNEL);
    change_set |= setup_parameter(anisotropy, data->anisotropy, ChangeSet::GENERIC_ATTRIBUTE);
    change_set |= setup_parameter(blackbody_intensity, data->blackbody_intensity, ChangeSet::GENERIC_ATTRIBUTE);
    change_set |= setup_parameter(emission_mode, data->emission_mode, ChangeSet::GENERIC_ATTRIBUTE);

    change_set |= setup_parameter(vdb_path, data->vdb_path, ChangeSet::VDB_FILE);
    vdb_file = data->vdb_file;
    if (!vdb_file)
        clear();

    if (display_mode == VDBDisplayMode::DISPLAY_SLICES)
    {
        change_set |= setup_parameter(max_slice_count, data->max_slice_count, ChangeSet::MAX_SLICE_COUNT);
        change_set |= setup_parameter(shadow_sample_count, data->shadow_sample_count, ChangeSet::GENERIC_ATTRIBUTE);
        change_set |= setup_parameter(shadow_gain, data->shadow_gain, ChangeSet::GENERIC_ATTRIBUTE);
    }

    return change_set != ChangeSet::NO_CHANGES;
}

// === VDBSubSceneOverride implementation ==================================

void VDBSubSceneOverride::shader_instance_deleter::operator()(MHWRender::MShaderInstance* p)
{
    auto shmgr = get_shader_manager();
    if (shmgr != nullptr)
        shmgr->releaseShader(p);
}

MString VDBSubSceneOverride::registrantId("VDBVisualizerSubSceneOverride");

MHWRender::MPxSubSceneOverride* VDBSubSceneOverride::creator(const MObject& obj)
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
}

VDBSubSceneOverride::~VDBSubSceneOverride()
{
}

MHWRender::DrawAPI VDBSubSceneOverride::supportedDrawAPIs() const
{
#if MAYA_API_VERSION >= 201600
    return MHWRender::kOpenGLCoreProfile;
#else
    return MHWRender::kOpenGL;
#endif
}

void VDBSubSceneOverride::update(MHWRender::MSubSceneContainer& container, const MHWRender::MFrameContext& /*frameContext*/)
{
    VDBSubSceneOverrideData* data = p_data.get();

    const MHWRender::MShaderManager* shader_manager = get_shader_manager();
    if (shader_manager == nullptr)
        return;

    MHWRender::MRenderItem* bounding_box = container.find("bounding_box");
    if (bounding_box == nullptr)
    {
        bounding_box = MHWRender::MRenderItem::Create("bounding_box",
            MHWRender::MRenderItem::NonMaterialSceneItem,
            MHWRender::MGeometry::kLines);
        bounding_box->enable(false);
        bounding_box->setDrawMode(MHWRender::MGeometry::kAll);
        bounding_box->depthPriority(MHWRender::MRenderItem::sDormantWireDepthPriority);

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

    bounding_box->setMatrix(&data->world_matrix);


    MHWRender::MRenderItem* selection_bounding_box = container.find("selection_bounding_box");
    if (selection_bounding_box == nullptr)
    {
        selection_bounding_box = MHWRender::MRenderItem::Create("selection_bounding_box",
            MHWRender::MRenderItem::NonMaterialSceneItem,
            MHWRender::MGeometry::kTriangles);
        selection_bounding_box->enable(true);
        selection_bounding_box->setDrawMode(MHWRender::MGeometry::kSelectionOnly);
        selection_bounding_box->depthPriority(MHWRender::MRenderItem::sSelectionDepthPriority);

        MHWRender::MShaderInstance* shader = shader_manager->getStockShader(
            MHWRender::MShaderManager::k3dSolidShader, nullptr, nullptr);
        if (shader)
        {
            selection_bounding_box->setShader(shader);
        }

        container.add(selection_bounding_box);
    }

    selection_bounding_box->setMatrix(&data->world_matrix);


    if (data->change_set == ChangeSet::NO_CHANGES)
        return;

    const bool file_exists = data->vdb_file != nullptr;

    const static MHWRender::MVertexBufferDescriptor position_buffer_desc("", MHWRender::MGeometry::kPosition, MHWRender::MGeometry::kFloat, 3);
    const static MHWRender::MVertexBufferDescriptor color_buffer_desc("", MHWRender::MGeometry::kColor, MHWRender::MGeometry::kFloat, 4);

    // Box vertex positions.
    p_bbox_position.reset(new MHWRender::MVertexBuffer(position_buffer_desc));
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
    }

    // Selection bbox.
    {
        p_selection_bbox_indices.reset(new MHWRender::MIndexBuffer(MHWRender::MGeometry::kUnsignedInt32));
        set_bbox_indices_triangles(1, p_selection_bbox_indices.get());

        MHWRender::MVertexBufferArray vertex_buffers;
        vertex_buffers.addBuffer("", p_bbox_position.get());
        setGeometryForRenderItem(*selection_bounding_box, vertex_buffers, *p_selection_bbox_indices.get(), &data->bbox);
        selection_bounding_box->enable(data->is_visible);
    }

    if (!file_exists || data->display_mode <= DISPLAY_GRID_BBOX)
    {
        bounding_box->enable(data->is_visible);
        m_sliced_display->enable(false);

        MHWRender::MVertexBufferArray vertex_buffers;
        p_bbox_indices.reset(new MHWRender::MIndexBuffer(MHWRender::MGeometry::kUnsignedInt32));

        if ((data->display_mode == DISPLAY_AXIS_ALIGNED_BBOX) || !file_exists)
        {
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
                    p_bbox_position.reset(new MHWRender::MVertexBuffer(position_buffer_desc));
                    p_bbox_indices.reset(new MHWRender::MIndexBuffer(MHWRender::MGeometry::kUnsignedInt32));
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
    else if (data->display_mode == DISPLAY_SLICES)
    {
        bounding_box->enable(false);

        m_sliced_display->enable(data->is_visible);

        if (data->is_visible)
            m_sliced_display->update(container, *data);
    }
    else
    {
        m_sliced_display->enable(false);
    }

    data->change_set = ChangeSet::NO_CHANGES;
}

bool VDBSubSceneOverride::requiresUpdate(const MHWRender::MSubSceneContainer& /*container*/, const MHWRender::MFrameContext& /*frameContext*/) const
{
    return p_data->update(p_vdb_visualizer->get_update(), m_object);
}

// === Sliced display mode implementation ===================================

SlicedDisplay::SlicedDisplay(MHWRender::MPxSubSceneOverride& parent)
    : m_parent(parent),
    m_density_channel("density"), m_scattering_channel("scattering"), m_emission_channel("emission"), m_transparency_channel("transparency"), m_temperature_channel("temperature"),
    m_density_ramp(RAMP_RESOLUTION), m_scattering_ramp(RAMP_RESOLUTION), m_emission_ramp(RAMP_RESOLUTION),
    m_volume_sampler_state(MHWRender::MSamplerState::kMinMagMipLinear, MHWRender::MSamplerState::kTexBorder),
    m_enabled(false), m_selected(false), m_force_shadows_off(false)
{
    if (!m_volume_shader)
        return;

    VolumeCache::registerUsage();

    m_density_channel.setShaderInstance(m_volume_shader.get());
    m_scattering_channel.setShaderInstance(m_volume_shader.get());
    m_emission_channel.setShaderInstance(m_volume_shader.get());
    m_transparency_channel.setShaderInstance(m_volume_shader.get());
    m_temperature_channel.setShaderInstance(m_volume_shader.get());

    // Create sampler state for textures.
    for (MString param : { "density_sampler", "scattering_sampler", "transparency_sampler", "emission_sampler", "temperature_sampler" })
        m_volume_sampler_state.assign(m_volume_shader.get(), param);

    // Assign ramp sampler states.
    m_density_ramp.assignSamplerToShader(m_volume_shader.get(), "density_ramp_sampler");
    m_scattering_ramp.assignSamplerToShader(m_volume_shader.get(), "scattering_ramp_sampler");
    m_emission_ramp.assignSamplerToShader(m_volume_shader.get(), "emission_ramp_sampler");

    // Set up blackbody LUT texture.
    m_blackbody_lut.lut.assignSamplerToShader(m_volume_shader.get(), "blackbody_lut_sampler");
    m_blackbody_lut.lut.assignTextureToShader(m_volume_shader.get(), "blackbody_lut_texture");

    // Shadows don't display correctly on Radeon HD 7950 under Linux,
    // so switch them off if the OpenGL vendor string doesn't contain
    // the string "nvidia" as a temporary (TM) workaround (TM).
    const std::string gl_vendor = reinterpret_cast<const char*>(glGetString(GL_VENDOR));
    static const std::string pattern = "nvidia";
    auto it = std::search(
        gl_vendor.begin(), gl_vendor.end(), pattern.begin(), pattern.end(),
        [](char a, char b) { return std::tolower(a, std::locale()) == std::tolower(b, std::locale()); });
    m_force_shadows_off = (it == gl_vendor.end());
}

SlicedDisplay::~SlicedDisplay()
{
    VolumeCache::unregisterUsage();
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

        const MHWRender::MShaderManager* shader_manager = get_shader_manager();
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
        const auto z = float(i);
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
    constexpr size_t vertex_count = 8;
    MFloatVector* positions = static_cast<MFloatVector*>(m_bbox_renderable.position_buffer->acquire(vertex_count, true));
    for (size_t i = 0; i < vertex_count; ++i) {
        positions[i] = bbox.min() + MVector((i & 1) * extents.x, ((i & 2) >> 1) * extents.y, ((i & 4) >> 2) * extents.z);
    }
    m_bbox_renderable.position_buffer->commit(positions);

    m_bbox_renderable.update(m_parent, bbox);

    CHECK_MSTATUS(m_volume_shader->setParameter("volume_origin", bbox.min()));
    CHECK_MSTATUS(m_volume_shader->setParameter("volume_size", extents));
}

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
    CHECK_MSTATUS(m_volume_shader->setParameter("density", data.density_channel.intensity));
    CHECK_MSTATUS(m_volume_shader->setParameter("density_source", int(data.density_channel.color_source)));
    CHECK_MSTATUS(m_volume_shader->setParameter("scattering_intensity", data.scattering_channel.intensity));
    CHECK_MSTATUS(m_volume_shader->setParameter("scattering_color", data.scattering_channel.color));
    CHECK_MSTATUS(m_volume_shader->setParameter("scattering_color_source", int(data.scattering_channel.color_source)));
    CHECK_MSTATUS(m_volume_shader->setParameter("scattering_anisotropy", data.anisotropy));
    CHECK_MSTATUS(m_volume_shader->setParameter("transparency", data.transparency_channel.color));
    CHECK_MSTATUS(m_volume_shader->setParameter("emission_intensity", data.emission_channel.intensity));
    CHECK_MSTATUS(m_volume_shader->setParameter("emission_color", data.emission_channel.color));
    CHECK_MSTATUS(m_volume_shader->setParameter("emission_color_source", int(data.emission_channel.color_source)));
    CHECK_MSTATUS(m_volume_shader->setParameter("emission_mode", int(data.emission_mode)));
    CHECK_MSTATUS(m_volume_shader->setParameter("temperature", data.temperature_channel.intensity));
    CHECK_MSTATUS(m_volume_shader->setParameter("blackbody_intensity", data.blackbody_intensity));
    CHECK_MSTATUS(m_volume_shader->setParameter("shadow_gain", data.shadow_gain));
    CHECK_MSTATUS(m_volume_shader->setParameter("shadow_sample_count", m_force_shadows_off ? 0 : data.shadow_sample_count));

    // === Update channels. ===

    // Bail if nothing has changed which affects the volume textures or the gradient ramps.
    if (data.change_set == ChangeSet::GENERIC_ATTRIBUTE)
        return true;

    // Update ramps.
    m_density_ramp.updateFromGradient(data.density_channel.gradient);
    m_density_ramp.assignTextureToShader(m_volume_shader.get(), "density_ramp_texture");
    m_density_ramp.assignDomainToShader(m_volume_shader.get(), "density_ramp_domain");
    m_scattering_ramp.updateFromGradient(data.scattering_channel.gradient);
    m_scattering_ramp.assignTextureToShader(m_volume_shader.get(), "scattering_ramp_texture");
    m_scattering_ramp.assignDomainToShader(m_volume_shader.get(), "scattering_ramp_domain");
    m_emission_ramp.updateFromGradient(data.emission_channel.gradient);
    m_emission_ramp.assignTextureToShader(m_volume_shader.get(), "emission_ramp_texture");
    m_emission_ramp.assignDomainToShader(m_volume_shader.get(), "emission_ramp_domain");

    // Bail if nothing has changed which affects the volume textures.
    if (data.change_set == ChangeSet::GRADIENT)
        return true;

    // Update slice geometry if number of slices has changed.
    if (hasChange(data.change_set, ChangeSet::MAX_SLICE_COUNT))
        updateSliceGeo(data);

    // Update file-level bbox.
    updateBBox(data.bbox);

    // Update volumes.
    const auto extents = openvdb::Coord(data.max_slice_count, data.max_slice_count, data.max_slice_count);
    m_density_channel.setVolume({ data.vdb_path, data.vdb_file->getUniqueTag(), data.density_channel.name, extents });
    m_scattering_channel.setVolume({ data.vdb_path, data.vdb_file->getUniqueTag(), data.scattering_channel.name, extents });
    m_emission_channel.setVolume({ data.vdb_path, data.vdb_file->getUniqueTag(), data.emission_channel.name, extents });
    m_transparency_channel.setVolume({ data.vdb_path, data.vdb_file->getUniqueTag(), data.transparency_channel.name, extents });
    m_temperature_channel.setVolume({ data.vdb_path, data.vdb_file->getUniqueTag(), data.temperature_channel.name, extents });

    return true;
}

MSyntax VDBVolumeCacheCmd::create_syntax()
{
    MSyntax syntax;
    syntax.enableQuery();
    syntax.enableEdit();
    syntax.addFlag("h", "help", MSyntax::kNoArg);
    syntax.addFlag("l", "limit", MSyntax::kLong);
    syntax.makeFlagQueryWithFullArgs("limit", true);
    syntax.addFlag("vt", "voxelType", MSyntax::kString);
    syntax.makeFlagQueryWithFullArgs("voxelType", true);
    return syntax;
}

namespace
{
    MString getVoxelTypeString()
    {
        const auto voxel_type = VolumeCache::instance().getVoxelType();
        if (voxel_type == VolumeCache::VoxelType::HALF)
            return "half";
        else if (voxel_type == VolumeCache::VoxelType::FLOAT)
            return "float";
        else
            return "unknown";
    }

} // unnamed namespace

MStatus VDBVolumeCacheCmd::doIt(const MArgList& args)
{
    MStatus status;
    MArgParser parser(syntax(), args);

    const auto display_error = [](const MString& message)
    {
        MGlobal::displayError(format("[openvdb_render] command ^1s: ^2s", COMMAND_STRING, message));
    };

    if (parser.isEdit())
    {
        if (parser.isFlagSet("limit"))
        {
            // Set volume cache limit to the given value in gigabytes.
            const int new_limit_gigabytes = parser.flagArgumentInt("limit", 0, &status);
            if (status != MStatus::kSuccess || new_limit_gigabytes < 0)
            {
                display_error("In edit mode argument to 'limit' has to be a non-negative integer representing gigabytes.");
                return MS::kFailure;
            }

            VolumeCache::instance().setMemoryLimitBytes(size_t(new_limit_gigabytes) << 30);
        }

        if (parser.isFlagSet("voxelType"))
        {
            const auto voxel_type_str = parser.flagArgumentString("voxelType", 0, &status);
            if (status != MStatus::kSuccess)
            {
                display_error("In edit mode the 'voxelType' flag requires a string argument, either 'half' or 'float'.");
                return MS::kFailure;
            }

            if (voxel_type_str == "half")
                VolumeCache::instance().setVoxelType(VolumeCache::VoxelType::HALF);
            else if (voxel_type_str == "float")
                VolumeCache::instance().setVoxelType(VolumeCache::VoxelType::FLOAT);
            else
            {
                display_error("In edit mode argument to 'voxelType' has to be either 'half' or 'float'.");
                return MS::kFailure;
            }
        }

        return MS::kSuccess;
    }
    else if (parser.isQuery())
    {
        if (parser.isFlagSet("limit"))
        {
            // Return volume cache limit in gigabytes.
            const size_t limit_bytes = VolumeCache::instance().getMemoryLimitBytes();
            MPxCommand::setResult(unsigned(limit_bytes / (1 << 30)));
            return MS::kSuccess;
        }
        else if (parser.isFlagSet("voxelType"))
        {
            // Return the voxel type as string.
            MPxCommand::setResult(getVoxelTypeString());
            return MS::kSuccess;
        }

        display_error("In query mode either 'limit' or 'voxelType' flag has to be specified.");
        return MS::kFailure;
    }

    // Neither edit nor query mode: display info.

    if (parser.isFlagSet("voxelType"))
    {
        // Display voxel type.
        MGlobal::displayInfo(format("[openvdb_render] volume cache voxel type is '^1s'.", getVoxelTypeString()));
        return MS::kSuccess;
    }
    else if (parser.isFlagSet("limit"))
    {
        // Display allocated bytes and cache limit.
        const auto pretty_string_size = [](size_t size) -> std::string
        {
            static const std::array<char, 3> prefixes = { { 'G', 'M', 'K' } };

            std::stringstream ss;
            ss << std::setprecision(2) << std::setiosflags(std::ios_base::fixed);

            for (size_t i = 0; i < prefixes.size(); ++i)
            {
                size_t prefix_size = 1LL << (10 * (prefixes.size() - i));
                if (size < prefix_size)
                    continue;

                ss << double(size) / double(prefix_size) << prefixes[i] << "B";
                return ss.str();
            }

            ss << double(size) << "B";
            return ss.str();
        };

        const size_t limit = VolumeCache::instance().getMemoryLimitBytes();
        if (limit == 0)
        {
            MGlobal::displayInfo("[openvdb_render] Volume caching is off.");
            return MS::kSuccess;
        }

        const size_t alloc = VolumeCache::instance().getAllocatedBytes();
        MGlobal::displayInfo(format("[openvdb_render] Volume cache allocated/total: ^1s/^2s.",
            pretty_string_size(alloc),
            pretty_string_size(limit)));
        return MS::kSuccess;
    }

    // Default: display help.
    MGlobal::displayInfo(format("[openvdb_render] Usage: ^1s [-h|-help] [-q|-query|-e|-edit] [-vt|-voxelType [\"half\"|\"float\"]] [-l|-limit [<limit_in_gigabytes>]]", COMMAND_STRING));
    return MS::kSuccess;
}

const char *VDBVolumeCacheCmd::COMMAND_STRING = "vdb_visualizer_volume_cache";
