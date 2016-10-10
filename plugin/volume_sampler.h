#pragma once

#include <openvdb/openvdb.h>
#include <openvdb/tools/MultiResGrid.h>

namespace MHWRender {
    class MTexture;
    class MTextureManager;
}
MHWRender::MTexture* volumeTextureFromGrid(const openvdb::tools::MultiResGrid<openvdb::FloatTree>& grid, const openvdb::Coord& slice_counts, MHWRender::MTextureManager* texture_manager);
