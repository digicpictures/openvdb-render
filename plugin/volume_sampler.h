#pragma once

#include <new>

#include <openvdb/openvdb.h>
#include <openvdb/tools/MultiResGrid.h>

#include "vdb_subscene_utils.hpp"

// The default ctor shouldn't initialize anything so this struct can be
// placement new'd onto a buffer without changing the contents.
template <typename RealType>
struct VolumeBufferHeader
{
    RealType value_range[2];
    RealType size[3];
    RealType origin[3];
};

template <typename RealType>
struct VolumeBufferHandle
{
    typedef VolumeBufferHeader<RealType> HeaderType;

    HeaderType *header;
    RealType *voxel_array;

    VolumeBufferHandle() : header(nullptr), voxel_array(nullptr) {}
    VolumeBufferHandle(void *raw_buffer) :
            header(new (raw_buffer) HeaderType),
            voxel_array(reinterpret_cast<RealType*>(raw_buffer) + sizeof(HeaderType) / sizeof(RealType))
    {}

    VolumeBufferHandle(const VolumeBufferHandle&) = delete;
    VolumeBufferHandle& operator=(const VolumeBufferHandle&) = delete;
    VolumeBufferHandle(VolumeBufferHandle&&) = default;
    VolumeBufferHandle& operator=(VolumeBufferHandle&&) = default;

    operator bool() const { return header && voxel_array; }
};

class ProgressBar;

// Sample a 3D domain (e.g. openvdb grid) using a sampling function (e.g. box filter),
// and create a Maya texture using the samples.
template <typename SampleType>
class VolumeSampler
{
public:
    typedef VolumeBufferHandle<SampleType> VolumeBufferHandleType;

    VolumeSampler() : m_progress_bar(nullptr) {}
    void attachProgressBar(ProgressBar *progress_bar) { m_progress_bar = progress_bar; }

    // Sample grid. Use MultiResGrid filtering only if texture_extents is coarser than the grid, otherwise use box filter.
    bool sampleGrid(const openvdb::FloatGrid& grid, const openvdb::Coord& texture_extents, const VolumeBufferHandleType& output);
    // Sample a single grid using simple box filter.
    bool sampleGridWithBoxFilter(const openvdb::FloatGrid& grid, const openvdb::Coord& texture_extents, const VolumeBufferHandleType& output);
    // Sample a MultiResGrid. Filtering is done by the built-in sampling mechanism of MultiResGrid.
    bool sampleMultiResGrid(const openvdb::tools::MultiResGrid<openvdb::FloatTree>& multires, const openvdb::Coord& texture_extents, const VolumeBufferHandleType& output);

private:
    ProgressBar *m_progress_bar;
};
