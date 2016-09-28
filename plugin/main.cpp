#include <GL/glew.h>

#include <boost/filesystem.hpp>

#include <maya/MFnPlugin.h>
#include <maya/MGlobal.h>
#include <maya/MDrawRegistry.h>
#include <maya/MShaderManager.h>

#include "vdb_geometry_override.h"
#include "vdb_visualizer.h"
#include "vdb_query.h"
#include "vdb_sampler.h"
#include "vdb_shader.h"
#include "vdb_simple_shader.h"


MHWRender::MShaderInstance *volume_shader = nullptr;
const MHWRender::MShaderManager* shader_manager = nullptr;
void loadVolumeShader()
{
    MHWRender::MRenderer* renderer = MHWRender::MRenderer::theRenderer();
    if (renderer == nullptr)
        return;

    shader_manager = renderer->getShaderManager();
    if (shader_manager == nullptr)
        return;

    using boost::filesystem::path;
    // TODO: use correct effect file path
    path effect_file = path(__FILE__).parent_path() / "volume.cgfx";
    volume_shader = shader_manager->getEffectsFileShader(effect_file.c_str(), "Main", 0, 0, false);
}
void unloadVolumeShader()
{
    shader_manager->releaseShader(volume_shader);
}

__declspec(dllexport) MStatus initializePlugin(MObject obj)
{
    const bool is_interactive = MGlobal::mayaState() == MGlobal::kInteractive;
    MStatus status = MS::kFailure;

    if (is_interactive)
    {
        /*if (glewInit() != GLEW_OK)
        {
            status.perror("[openvdb] Error initializing glew.");
            return status;
        }

        if (!glewIsSupported("GL_EXT_direct_state_access"))
        {
            status.perror("[openvdb] Direct State Access is not available, update your drivers or use a newer GPU!");
            return status;
        }*/
    }

    MFnPlugin plugin(obj, "Luma Pictures", "0.0.1", "Any");

    status = plugin.registerShape(VDBVisualizerShape::typeName, VDBVisualizerShape::typeId,
                                  VDBVisualizerShape::creator, VDBVisualizerShape::initialize,
                                  VDBVisualizerShapeUI::creator, &VDBVisualizerShape::drawDbClassification);

    if (!status)
    {
        status.perror("[openvdb] Error registering the VDBVisualizer Node.");
        return status;
    }

    status = plugin.registerNode(VDBSamplerNode::s_type_name, VDBSamplerNode::s_type_id,
                                 VDBSamplerNode::creator, VDBSamplerNode::initialize, MPxNode::kDependNode, &VDBSamplerNode::s_classification);

    if (!status)
    {
        status.perror("[openvdb] Error registering the VDBSampler Node.");
        return status;
    }

    status = plugin.registerNode(VDBShaderNode::s_type_name, VDBShaderNode::s_type_id,
                                 VDBShaderNode::creator, VDBShaderNode::initialize, MPxNode::kDependNode, &VDBShaderNode::s_classification);

    if (!status)
    {
        status.perror("[openvdb] Error registering the VDBShader Node.");
        return status;
    }

    status = plugin.registerNode(VDBSimpleShaderNode::s_type_name, VDBSimpleShaderNode::s_type_id,
                                 VDBSimpleShaderNode::creator, VDBSimpleShaderNode::initialize, MPxNode::kDependNode, &VDBSimpleShaderNode::s_classification);

    if (!status)
    {
        status.perror("[openvdb] Error registering the VDBSimpleShader Node.");
        return status;
    }

    status = MHWRender::MDrawRegistry::registerGeometryOverrideCreator(
            VDBVisualizerShape::drawDbClassification,
            MHWRender::VDBGeometryOverride::registrantId,
            MHWRender::VDBGeometryOverride::creator
    );

    if (!status)
    {
        status.perror("[openvdb] Error registering the VDBVisualizer Geometry Override.");
        return status;
    }

    openvdb::initialize();

    status = plugin.registerCommand("vdb_query", VDBQueryCmd::creator, VDBQueryCmd::create_syntax);

    if (!status)
    {
        status.perror("[openvdb] Error registering the VDBQuery Command.");
        return status;
    }

    if (is_interactive)
        MGlobal::executePythonCommand("import AEvdb_visualizerTemplate; import AEvdb_samplerTemplate; import AEvdb_shaderTemplate");

    // Load volume shader.
    loadVolumeShader();

    return status;
}

__declspec(dllexport) MStatus uninitializePlugin(MObject obj)
{
    MStatus status = MS::kSuccess;

    MFnPlugin plugin(obj);

    status = plugin.deregisterNode(VDBVisualizerShape::typeId);

    if (!status)
    {
        status.perror("[openvdb] Error deregistering the VDBVisualizer Node.");
        return status;
    }

    status = plugin.deregisterNode(VDBSamplerNode::s_type_id);

    if (!status)
    {
        status.perror("[openvdb] Error deregistering the VDBSampler Node.");
        return status;
    }

    status = plugin.deregisterNode(VDBShaderNode::s_type_id);

    if (!status)
    {
        status.perror("[openvdb] Error deregistering the VDBShader Node.");
        return status;
    }

    status = plugin.deregisterNode(VDBSimpleShaderNode::s_type_id);

    if (!status)
    {
        status.perror("[openvdb] Error deregistering the VDBSimpleShader Node.");
        return status;
    }

    status = MHWRender::MDrawRegistry::deregisterGeometryOverrideCreator(
            VDBVisualizerShape::drawDbClassification,
            MHWRender::VDBGeometryOverride::registrantId
    );

    if (!status)
    {
        status.perror("[openvdb] Error deregistering the VDBVisualizer Geometry Override.");
        return status;
    }

    status = plugin.deregisterCommand("vdb_query");

    if (!status)
    {
        status.perror("[openvdb] Error deregistering the VDBQuery Command.");
        return status;
    }

    // Unload volume shader.
    unloadVolumeShader();

    openvdb::uninitialize();

    return status;
}
