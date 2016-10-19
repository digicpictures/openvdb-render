#pragma once

#include <cstdlib>
#include <algorithm>
#include <limits>
#include <string>
#include <sstream>
#include <utility>
#include <maya/MFloatVector.h>
#include <openvdb/Types.h>

openvdb::CoordBBox inline
getIndexSpaceBoundingBox(const openvdb::GridBase* grid)
{
    const auto file_bbox_min = openvdb::Coord(grid->metaValue<openvdb::Vec3i>("file_bbox_min"));
    if (file_bbox_min.x() == std::numeric_limits<int>::max() ||
        file_bbox_min.y() == std::numeric_limits<int>::max() ||
        file_bbox_min.z() == std::numeric_limits<int>::max()) {
        return {};
    }
    const auto file_bbox_max = openvdb::Coord(grid->metaValue<openvdb::Vec3i>("file_bbox_max"));
    if (file_bbox_max.x() == std::numeric_limits<int>::min() ||
        file_bbox_max.y() == std::numeric_limits<int>::min() ||
        file_bbox_max.z() == std::numeric_limits<int>::min()) {
        return {};
    }

    return { file_bbox_min, file_bbox_max };
}


template <typename T>
struct ValueRange
{
    T min, max;
    ValueRange() noexcept : min(std::numeric_limits<T>::max()), max(std::numeric_limits<T>::min()) {}
    ValueRange(T min_ , T max_) noexcept : min(min_), max(max_) {}
    void update(T value) { min = std::min(min, value); max = std::max(max, value); }
    void update(const ValueRange<T>& value_range) { min = std::min(min, value_range.min); max = std::max(max, value_range.max); }
    ValueRange merge(const ValueRange<T>& value_range) { return { std::min(min, value_range.min), std::max(max, value_range.max) }; }
};
using FloatRange = ValueRange<float>;


template <typename T>
class CachedData {
public:
    template <typename... Args>
    CachedData(Args&&... args) : m_data(std::forward<Args>(args)...), m_dirty(false) {}
    const T& peek() const { return m_data; }
    bool isDirty() const { return m_dirty; }
    T& getDirtyRef() { m_dirty = true; return m_data; }
    const T& clearAndGet() const { m_dirty = false; return m_data; }

private:
    T m_data;
    mutable bool m_dirty;
};

class MayaPathSpec
{
public:
    MayaPathSpec(const std::string& path_spec) : m_path_spec(path_spec)
    {
        m_field_begin = m_path_spec.find(FIELD_CHAR, 0);
        const auto field_end = m_path_spec.find_first_not_of(FIELD_CHAR, m_field_begin + 1);
        m_field_length = field_end - m_field_begin;
    }
    MayaPathSpec(const MayaPathSpec&) = default;
    MayaPathSpec(MayaPathSpec&&) = default;
    MayaPathSpec& operator=(const MayaPathSpec&) = default;
    MayaPathSpec& operator=(MayaPathSpec&&) = default;

    bool hasFrameField() const { return m_field_begin != std::string::npos; }

    std::string getPath() const { return m_path_spec; }
    std::string getPath(int frame_num) const
    {
        assert(hasFrameField());

        std::stringstream ss;
        ss.fill('0');
        ss.width(m_field_length);
        ss << frame_num;

        std::string path = m_path_spec;
        path.replace(m_field_begin, m_field_length, ss.str());
        return path;
    }

private:
    std::string m_path_spec;
    size_t m_field_begin, m_field_length;

    static constexpr char FIELD_CHAR = '#';
};

#define LOG_ERROR(msg) ::log_error(msg, __FILE__, __LINE__)
inline void log_error(const std::string& msg, const char *file_name, int line_no)
{
#if _DEBUG
    std::cerr << "openvdb_render error: " << file_name << ": line " << line_no << ": " << msg << std::endl;
#else
    std::cerr << "openvdb_render error: " << msg << std::endl;
#endif
}
