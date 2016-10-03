#include "vdb_draw_override.h"

#include <maya/MDagPath.h>
#include <maya/MDrawContext.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MGLFunctionTable.h>
#include <maya/MHardwareRenderer.h>
#include <maya/MUserData.h>

#include "vdb_visualizer.h"

class GLFunctionTable
{
public:
    GLFunctionTable() : m_glft(nullptr) {}

    MGLFunctionTable* operator->() const
    {
        if (!m_glft) {
            MHardwareRenderer* renderer = MHardwareRenderer::theRenderer();
            if (!renderer) {
                return nullptr;
            }
            m_glft = renderer->glFunctionTable();
        }
        return m_glft;
    }

private:
    mutable MGLFunctionTable* m_glft;

};

GLFunctionTable g_glft;
typedef VDBVisualizerShape ShapeType;

// === UserData ===========================================

class VDBDrawOverride::UserData : public MUserData
{
public:
    UserData(const ShapeType* shape) : MUserData(false), m_shape(shape) {}
    virtual ~UserData() override {}

    void draw(const MHWRender::MDrawContext& context) const;

private:
    const ShapeType* m_shape;
};

void VDBDrawOverride::UserData::draw(const MHWRender::MDrawContext& context) const
{
    using namespace MHWRender;

    MRenderer* renderer = MRenderer::theRenderer();
    if (!renderer) return;

    // Only OpenGL is supported.
    if (!renderer->drawAPIIsOpenGL()) {
        return;
    }

    MStateManager* stateMgr = context.getStateManager();
    if (!stateMgr) return;

    const unsigned int displayStyle = context.getDisplayStyle();
    if (displayStyle == 0) return;

    if (displayStyle & MFrameContext::kXray) {
        // Viewport 2.0 will call draw() twice when drawing transparent objects
        // (X-Ray mode). We skip the first draw() call.
        const MRasterizerState* rasterState = stateMgr->getRasterizerState();
        if (rasterState && rasterState->desc().cullMode == MRasterizerState::kCullFront) {
            return;
        }
    }

    MStatus status;
    const MMatrix wvp_matrix = context.getMatrix(MFrameContext::kWorldViewProjMtx, &status);
    if (status != MStatus::kSuccess) return;

    //const MStringArray& pass_semantics = context.passS
    //std::find_if()
}

// === VDBDrawOverride ====================================

MHWRender::MPxDrawOverride* VDBDrawOverride::creator(const MObject& obj)
{
    return new VDBDrawOverride(obj);
}

void VDBDrawOverride::drawCb(const MHWRender::MDrawContext& drawContext, const MUserData* mdata)
{
    const UserData* data = dynamic_cast<const UserData*>(mdata);
    if (!data) {
        return;
    }
    data->draw(drawContext);
}

MBoundingBox VDBDrawOverride::boundingBox(const MDagPath& objPath, const MDagPath& cameraPath) const
{
    // Try to get bounding box data from the VDB shape node at objPath.

    MStatus status;
    MFnDependencyNode node(objPath.node(), &status);
    if (!status) return MBoundingBox();

    const ShapeType* shapeNode = dynamic_cast<ShapeType*>(node.userNode());
    if (!shapeNode) return MBoundingBox();

    return shapeNode->boundingBox();
}

MUserData* VDBDrawOverride::prepareForDraw(const MDagPath& objPath, const MDagPath& cameraPath, const MHWRender::MFrameContext& frameContext, MUserData* oldData)
{
    MObject object = objPath.node();

    // Get user data from prev call or create a new one.
    UserData* data = dynamic_cast<UserData*>(oldData);
    if (!data) {
        MStatus status;
        MFnDependencyNode node(object, &status);
        if (!status) return nullptr;

        data = new UserData(dynamic_cast<ShapeType*>(node.userNode()));
    }

    // TODO: set data required for rendering.

    return data;
}
