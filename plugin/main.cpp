#include <GL/glew.h>

#include <maya/MFnPlugin.h>
#include <maya/MGlobal.h>
#include <maya/MDrawRegistry.h>
#include <maya/MShaderManager.h>

#include "vdb_subscene_override.h"
#include "vdb_visualizer.h"
#include "vdb_query.h"
#include "vdb_sampler.h"
#include "vdb_shader.h"
#include "vdb_simple_shader.h"
#include "paths.h"

PLUGIN_EXPORT MStatus initializePlugin(MObject obj)
{
    const bool is_interactive = MGlobal::mayaState() == MGlobal::kInteractive;
    MStatus status = MS::kFailure;

    MFnPlugin plugin(obj, "Luma Pictures", "0.0.1", "Any");

    Paths::init(plugin.loadPath());

    // Register VDBVisualizer shape node.
    status = plugin.registerShape(VDBVisualizerShape::typeName, VDBVisualizerShape::typeId,
                                  VDBVisualizerShape::creator, VDBVisualizerShape::initialize,
                                  VDBVisualizerShapeUI::creator, &VDBVisualizerShape::drawDbClassification);
    if (!status)
    {
        status.perror("[openvdb] Error registering the VDBVisualizer Node.");
        return status;
    }

    // Register VDBSampler node.
    status = plugin.registerNode(VDBSamplerNode::s_type_name, VDBSamplerNode::s_type_id,
                                 VDBSamplerNode::creator, VDBSamplerNode::initialize, MPxNode::kDependNode, &VDBSamplerNode::s_classification);
    if (!status)
    {
        status.perror("[openvdb] Error registering the VDBSampler Node.");
        return status;
    }

    // Register VDBShader node.
    status = plugin.registerNode(VDBShaderNode::s_type_name, VDBShaderNode::s_type_id,
                                 VDBShaderNode::creator, VDBShaderNode::initialize, MPxNode::kDependNode, &VDBShaderNode::s_classification);
    if (!status)
    {
        status.perror("[openvdb] Error registering the VDBShader Node.");
        return status;
    }

    // Register VDBSimpleShader node.
    status = plugin.registerNode(VDBSimpleShaderNode::s_type_name, VDBSimpleShaderNode::s_type_id,
                                 VDBSimpleShaderNode::creator, VDBSimpleShaderNode::initialize, MPxNode::kDependNode, &VDBSimpleShaderNode::s_classification);
    if (!status)
    {
        status.perror("[openvdb] Error registering the VDBSimpleShader Node.");
        return status;
    }

    // Register subscene override for shape node.
    status = MHWRender::MDrawRegistry::registerSubSceneOverrideCreator(
        VDBVisualizerShape::drawDbClassification,
        MHWRender::VDBSubSceneOverride::registrantId,
        MHWRender::VDBSubSceneOverride::creator
    );
    if (!status)
    {
        status.perror("[openvdb] Error registering the VDBVisualizer Sub Scene Override.");
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

    return status;
}

PLUGIN_EXPORT MStatus uninitializePlugin(MObject obj)
{
    MStatus status = MS::kSuccess;

    MFnPlugin plugin(obj);

    // Deregister VDBVisualizer shape node.
    status = plugin.deregisterNode(VDBVisualizerShape::typeId);
    if (!status)
    {
        status.perror("[openvdb] Error deregistering the VDBVisualizer Node.");
        return status;
    }

    // Deregister VDBSampler node.
    status = plugin.deregisterNode(VDBSamplerNode::s_type_id);
    if (!status)
    {
        status.perror("[openvdb] Error deregistering the VDBSampler Node.");
        return status;
    }

    // Deregister shader node.
    status = plugin.deregisterNode(VDBShaderNode::s_type_id);
    if (!status)
    {
        status.perror("[openvdb] Error deregistering the VDBShader Node.");
        return status;
    }

    // Deregister simple shader node.
    status = plugin.deregisterNode(VDBSimpleShaderNode::s_type_id);
    if (!status)
    {
        status.perror("[openvdb] Error deregistering the VDBSimpleShader Node.");
        return status;
    }

    // Deregister subscene override
    status = MHWRender::MDrawRegistry::deregisterSubSceneOverrideCreator(
            VDBVisualizerShape::drawDbClassification,
            MHWRender::VDBSubSceneOverride::registrantId
    );
    if (!status)
    {
        status.perror("[openvdb] Error deregistering the VDBVisualizer Sub Scene Override.");
        return status;
    }

    status = plugin.deregisterCommand("vdb_query");
    if (!status)
    {
        status.perror("[openvdb] Error deregistering the VDBQuery Command.");
        return status;
    }

    openvdb::uninitialize();

    return status;
}
