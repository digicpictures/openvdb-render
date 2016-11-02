#pragma once

#include <array>

#include <openvdb/openvdb.h>

#include <maya/MBoundingBox.h>
#include <maya/MMatrix.h>

inline bool read_grid_transformed_bbox_wire(openvdb::GridBase::ConstPtr grid, std::array<MFloatVector, 8>& vertices)
{
    const openvdb::Vec3i file_bbox_min = grid->metaValue<openvdb::Vec3i>("file_bbox_min");
    if (file_bbox_min.x() == std::numeric_limits<int>::max() ||
        file_bbox_min.y() == std::numeric_limits<int>::max() ||
        file_bbox_min.z() == std::numeric_limits<int>::max())
        return false;
    const openvdb::Vec3i file_bbox_max = grid->metaValue<openvdb::Vec3i>("file_bbox_max");
    if (file_bbox_max.x() == std::numeric_limits<int>::min() ||
        file_bbox_max.y() == std::numeric_limits<int>::min() ||
        file_bbox_max.z() == std::numeric_limits<int>::min())
        return false;
    const openvdb::math::Transform& transform = grid->transform();

    openvdb::Vec3d pnt = transform.indexToWorld(openvdb::Vec3i(file_bbox_min.x(), file_bbox_min.y(), file_bbox_min.z()));
    vertices[0] = MFloatVector(static_cast<float>(pnt.x()), static_cast<float>(pnt.y()), static_cast<float>(pnt.z()));

    pnt = transform.indexToWorld(openvdb::Vec3i(file_bbox_min.x(), file_bbox_max.y(), file_bbox_min.z()));
    vertices[1] = MFloatVector(static_cast<float>(pnt.x()), static_cast<float>(pnt.y()), static_cast<float>(pnt.z()));

    pnt = transform.indexToWorld(openvdb::Vec3i(file_bbox_min.x(), file_bbox_max.y(), file_bbox_max.z()));
    vertices[2] = MFloatVector(static_cast<float>(pnt.x()), static_cast<float>(pnt.y()), static_cast<float>(pnt.z()));

    pnt = transform.indexToWorld(openvdb::Vec3i(file_bbox_min.x(), file_bbox_min.y(), file_bbox_max.z()));
    vertices[3] = MFloatVector(static_cast<float>(pnt.x()), static_cast<float>(pnt.y()), static_cast<float>(pnt.z()));

    pnt = transform.indexToWorld(openvdb::Vec3i(file_bbox_max.x(), file_bbox_min.y(), file_bbox_min.z()));
    vertices[4] = MFloatVector(static_cast<float>(pnt.x()), static_cast<float>(pnt.y()), static_cast<float>(pnt.z()));

    pnt = transform.indexToWorld(openvdb::Vec3i(file_bbox_max.x(), file_bbox_max.y(), file_bbox_min.z()));
    vertices[5] = MFloatVector(static_cast<float>(pnt.x()), static_cast<float>(pnt.y()), static_cast<float>(pnt.z()));

    pnt = transform.indexToWorld(openvdb::Vec3i(file_bbox_max.x(), file_bbox_max.y(), file_bbox_max.z()));
    vertices[6] = MFloatVector(static_cast<float>(pnt.x()), static_cast<float>(pnt.y()), static_cast<float>(pnt.z()));

    pnt = transform.indexToWorld(openvdb::Vec3i(file_bbox_max.x(), file_bbox_min.y(), file_bbox_max.z()));
    vertices[7] = MFloatVector(static_cast<float>(pnt.x()), static_cast<float>(pnt.y()), static_cast<float>(pnt.z()));

    return true;
}

inline bool read_transformed_bounding_box(openvdb::GridBase::ConstPtr grid, MBoundingBox& bbox)
{
    const openvdb::Vec3i file_bbox_min = grid->metaValue<openvdb::Vec3i>("file_bbox_min");
    if (file_bbox_min.x() == std::numeric_limits<int>::max() ||
        file_bbox_min.y() == std::numeric_limits<int>::max() ||
        file_bbox_min.z() == std::numeric_limits<int>::max())
        return false;
    const openvdb::Vec3i file_bbox_max = grid->metaValue<openvdb::Vec3i>("file_bbox_max");
    if (file_bbox_max.x() == std::numeric_limits<int>::min() ||
        file_bbox_max.y() == std::numeric_limits<int>::min() ||
        file_bbox_max.z() == std::numeric_limits<int>::min())
        return false;
    const openvdb::math::Transform& transform = grid->transform();

    openvdb::Vec3d pnt = transform.indexToWorld(openvdb::Vec3i(file_bbox_min.x(), file_bbox_min.y(), file_bbox_min.z()));
    bbox.expand(MPoint(pnt.x(), pnt.y(), pnt.z(), 1.0));

    pnt = transform.indexToWorld(openvdb::Vec3i(file_bbox_max.x(), file_bbox_min.y(), file_bbox_min.z()));
    bbox.expand(MPoint(pnt.x(), pnt.y(), pnt.z(), 1.0));

    pnt = transform.indexToWorld(openvdb::Vec3i(file_bbox_min.x(), file_bbox_max.y(), file_bbox_min.z()));
    bbox.expand(MPoint(pnt.x(), pnt.y(), pnt.z(), 1.0));

    pnt = transform.indexToWorld(openvdb::Vec3i(file_bbox_min.x(), file_bbox_min.y(), file_bbox_max.z()));
    bbox.expand(MPoint(pnt.x(), pnt.y(), pnt.z(), 1.0));

    pnt = transform.indexToWorld(openvdb::Vec3i(file_bbox_max.x(), file_bbox_max.y(), file_bbox_min.z()));
    bbox.expand(MPoint(pnt.x(), pnt.y(), pnt.z(), 1.0));

    pnt = transform.indexToWorld(openvdb::Vec3i(file_bbox_max.x(), file_bbox_min.y(), file_bbox_max.z()));
    bbox.expand(MPoint(pnt.x(), pnt.y(), pnt.z(), 1.0));

    pnt = transform.indexToWorld(openvdb::Vec3i(file_bbox_min.x(), file_bbox_max.y(), file_bbox_max.z()));
    bbox.expand(MPoint(pnt.x(), pnt.y(), pnt.z(), 1.0));

    pnt = transform.indexToWorld(openvdb::Vec3i(file_bbox_max.x(), file_bbox_max.y(), file_bbox_max.z()));
    bbox.expand(MPoint(pnt.x(), pnt.y(), pnt.z(), 1.0));

    return true;
}

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

