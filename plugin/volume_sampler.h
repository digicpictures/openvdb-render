#pragma once

#include <type_traits>
#include <vector>

#include <openvdb/openvdb.h>
#include <openvdb/tools/MultiResGrid.h>

namespace MHWRender {
    class MTexture;
    class MTextureManager;
}

struct Texture
{
    MHWRender::MTexture* texture;
    float min_value;
    float max_value;
};

class VolumeSampler
{
public:
    VolumeSampler(MHWRender::MTextureManager* texture_manager) : m_texture_manager(texture_manager) { initSampleJitter(); }

    // Sample a single grid. Filtering is done by getting multiple (jittered) samples.
    Texture sampleGrid(const openvdb::FloatGrid& grid, const openvdb::Coord& texture_extents);
    // Sample a multi res grid. Filtering is done by the built-in sampling mechanism of MultiResGrid.
    Texture sampleMultiResGrid(const openvdb::tools::MultiResGrid<openvdb::FloatTree>& multires, const openvdb::Coord& texture_extents);

private:
    MHWRender::MTexture* acquireVolumeTexture(const openvdb::Coord& texture_extents, const float* pixel_data);

    //std::vector<openvdb::Vec3d> createSampleJitter(int n_samp);

    std::vector<float>
    samplingLoop(const openvdb::CoordBBox& input_domain, const openvdb::Coord& output_extents,
                 std::function<float(openvdb::Vec3d)> sampling_func);

    MHWRender::MTextureManager* m_texture_manager;


    static constexpr int N_FILT_SAMP = 1;
    void initSampleJitter();
    std::vector<openvdb::Vec3d> m_sample_jitter;
};
