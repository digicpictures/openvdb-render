#pragma once

#include <maya/MPxDrawOverride.h>

class VDBDrawOverride : public MHWRender::MPxDrawOverride
{

public:
    static MHWRender::MPxDrawOverride* creator(const MObject& obj);

    class UserData;
    static void drawCb(const MHWRender::MDrawContext& drawContext, const MUserData* data);

    virtual ~VDBDrawOverride() override {}

    VDBDrawOverride(const VDBDrawOverride& obj) = delete;
    const VDBDrawOverride& operator=(const VDBDrawOverride& obj) = delete;

    virtual MHWRender::DrawAPI supportedDrawAPIs() const override
    { return MHWRender::kOpenGL; }

    virtual bool isBounded(
        const MDagPath& objPath,
        const MDagPath& cameraPath) const override
    { return true; }

    virtual MBoundingBox boundingBox(
        const MDagPath& objPath,
        const MDagPath& cameraPath) const override;

    virtual bool disableInternalBoundingBoxDraw() const override
    { return false; };

    virtual MUserData* prepareForDraw(
        const MDagPath& objPath,
        const MDagPath& cameraPath,
        const MHWRender::MFrameContext& frameContext,
        MUserData* oldData) override;

private:
    VDBDrawOverride(const MObject& obj) : MHWRender::MPxDrawOverride(obj, &drawCb) {}
};