#pragma once

#include <new>

#include <openvdb/Types.h>
#include <openvdb/tools/MultiResGrid.h>

#include "vdb_subscene_utils.hpp"

// The default ctor shouldn't initialize anything so this struct can be
// placement new'd onto a buffer without changing the contents.
struct VolumeBufferHeader
{
    float value_range[2];
    float size[3];
    float origin[3];
};

struct VolumeBufferHandle
{
    static const size_t VOXEL_ARRAY_OFFSET;
    VolumeBufferHeader *header;
    float *voxel_array;

    VolumeBufferHandle() : header(nullptr), voxel_array(nullptr) {}
    VolumeBufferHandle(float *raw_buffer) :
            header(new (raw_buffer) VolumeBufferHeader),
            voxel_array(raw_buffer + VOXEL_ARRAY_OFFSET)
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
class VolumeSampler
{
public:
    VolumeSampler() : m_progress_bar(nullptr) {}
    void attachProgressBar(ProgressBar *progress_bar) { m_progress_bar = progress_bar; }

    // Sample grid. Use MultiResGrid filtering only if texture_extents is coarser than the grid, otherwise use box filter.
    void sampleGrid(const openvdb::FloatGrid& grid, const openvdb::Coord& texture_extents, const VolumeBufferHandle& output);
    // Sample a single grid using simple box filter.
    void sampleGridWithBoxFilter(const openvdb::FloatGrid& grid, const openvdb::Coord& texture_extents, const VolumeBufferHandle& output);
    // Sample a MultiResGrid. Filtering is done by the built-in sampling mechanism of MultiResGrid.
    void sampleMultiResGrid(const openvdb::tools::MultiResGrid<openvdb::FloatTree>& multires, const openvdb::Coord& texture_extents, const VolumeBufferHandle& output);

private:
    ProgressBar *m_progress_bar;
};
