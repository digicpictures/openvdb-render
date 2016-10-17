#pragma once

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
