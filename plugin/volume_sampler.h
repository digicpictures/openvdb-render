#pragma once

#include <openvdb/openvdb.h>

namespace MHWRender {
    class MTexture;
    class MTextureManager;
}
MHWRender::MTexture* volumeTextureFromGrid(const openvdb::FloatGrid* grid, const openvdb::Coord& slice_counts, MHWRender::MTextureManager* texture_manager);
