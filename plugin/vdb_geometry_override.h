#pragma once

#include <maya/MPxGeometryOverride.h>

#include <memory>

#include "vdb_visualizer.h"

struct VDBGeometryOverrideData;
class VDBGeometryOverride : public MHWRender::MPxGeometryOverride
{
public:
    VDBGeometryOverride(const MObject& obj);
    ~VDBGeometryOverride();

    static MHWRender::MPxGeometryOverride* creator(const MObject& obj);

    void updateDG();

    void updateRenderItems(
        const MDagPath& path,
        MHWRender::MRenderItemList& list);

    void populateGeometry(
        const MHWRender::MGeometryRequirements& requirements,
        const MHWRender::MRenderItemList& renderItems,
        MHWRender::MGeometry& data);
    void cleanUp();

    static MString registrantId;
private:
    VDBVisualizerShape* p_vdb_visualizer;

    std::unique_ptr<VDBGeometryOverrideData> p_data;
};
