#pragma once

#include <string>
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

MFloatVector inline
mayavecFromVec3f(const openvdb::Vec3f& vec)
{
    return { vec.x(), vec.y(), vec.z() };
}

template <typename VecT>
typename VecT::value_type maxComponent(const VecT& v)
{
    return std::max(std::max(v.x(), v.y()), v.z());
}

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
