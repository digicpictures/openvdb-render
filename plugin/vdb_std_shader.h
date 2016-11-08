#pragma once

#include <maya/MPxNode.h>
#include <maya/MString.h>
#include <maya/MTypeId.h>
#include <maya/MPlug.h>

#include "vdb_sampler.h"

struct VDBShaderParams;
struct VDBVolumeStandardShaderParams
{
    VDBVolumeStandardShaderParams();

    void create_params(VDBShaderParams *shared_params);
    void affect_output(MObject& out_object);
    bool check_plug(const MPlug& plug);

    MObject density;
    MObject density_channel;
    VDBGradientParams density_gradient;

    MObject scattering_intensity;
    MObject scattering_color;
    MObject scattering_channel;
    MObject anisotropy;
    VDBGradientParams scattering_gradient;

    MObject transparency;
    MObject transparency_channel;
    VDBGradientParams transparency_gradient;

    MObject emission_intensity;
    MObject emission_color;
    MObject emission_channel;
    MObject emission_mode;
    VDBGradientParams emission_gradient;

    MObject temperature;
    MObject temperature_channel;
    VDBGradientParams temperature_gradient;
};
