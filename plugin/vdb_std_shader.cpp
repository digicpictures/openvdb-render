#include "vdb_std_shader.h"

#include <maya/MFnNumericAttribute.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnStringData.h>

#include "vdb_shader.h"

VDBVolumeStandardShaderParams::VDBVolumeStandardShaderParams()
    : density_gradient("density"), scattering_gradient("scattering"), transparency_gradient("transparency"),
    emission_gradient("emission"), temperature_gradient("temperature")
{
}

void VDBVolumeStandardShaderParams::create_params(VDBShaderParams *shared_params)
{
    MFnEnumAttribute eAttr;
    MFnNumericAttribute nAttr;
    MFnTypedAttribute tAttr;
    MFnStringData sData;

    // === Params unique to volume standard shading mode. ===

    // Density.

    density = nAttr.create("density", "density", MFnNumericData::kFloat);
    nAttr.setDefault(1.0f);
    nAttr.setMin(0.0f);
    nAttr.setSoftMax(1.0f);
    nAttr.setChannelBox(true);
    MPxNode::addAttribute(density);

    density_channel = tAttr.create("densityChannel", "density_channel", MFnData::kString);
    tAttr.setDefault(sData.create("density"));
    MPxNode::addAttribute(density_channel);

    density_gradient.create_params();

    // Transparency.

    transparency = nAttr.createColor("transparency", "transparency");
    nAttr.setDefault(1.0, 1.0, 1.0);
    MPxNode::addAttribute(transparency);

    transparency_channel = tAttr.create("transparencyChannel", "transparency_channel", MFnData::kString);
    tAttr.setDefault(sData.create("density"));
    MPxNode::addAttribute(transparency_channel);

    transparency_gradient.create_params();

    // Emission mode.

    emission_mode = eAttr.create("emissionMode", "emission_mode");
    eAttr.addField("Uniform", 0);
    eAttr.addField("Density", 1);
    eAttr.addField("Channel", 2);
    eAttr.addField("Blackbody", 3);
    eAttr.addField("Density and blackbody", 4);
    eAttr.setDefault(0);
    MPxNode::addAttribute(emission_mode);

    // Temperature.

    temperature = nAttr.create("temperature", "temperature", MFnNumericData::kFloat);
    nAttr.setDefault(5000.0f);
    nAttr.setMin(0.0f);
    nAttr.setSoftMax(20000.0f);
    nAttr.setChannelBox(true);
    MPxNode::addAttribute(temperature);

    temperature_channel = tAttr.create("temperatureChannel", "temperature_channel", MFnData::kString);
    tAttr.setDefault(sData.create("temperature"));
    MPxNode::addAttribute(temperature_channel);

    temperature_gradient.create_params();

    // === Params sharable with other shading modes. ===

    // Get sharable params from pointer if supplied...
    if (shared_params)
    {
        scattering_intensity = shared_params->scattering_intensity;
        scattering_color = shared_params->scattering_color;
        scattering_channel = shared_params->scattering_channel;
        anisotropy = shared_params->anisotropy;
        scattering_gradient = shared_params->scattering_gradient;

        emission_intensity = shared_params->emission_intensity;
        emission_color = shared_params->emission_color;
        emission_channel = shared_params->emission_channel;
        emission_gradient = shared_params->emission_gradient;
        return;
    }
    // ... oherwise create sharable params here.

    // Scattering.

    scattering_intensity = nAttr.create("scatteringIntensity", "scattering_intensity", MFnNumericData::kFloat);
    nAttr.setDefault(1.0f);
    nAttr.setMin(0.0f);
    nAttr.setSoftMax(1.0f);
    nAttr.setChannelBox(true);
    MPxNode::addAttribute(scattering_intensity);

    scattering_color = nAttr.createColor("scatteringColor", "scattering_color");
    nAttr.setDefault(1.0, 1.0, 1.0);
    MPxNode::addAttribute(scattering_color);

    scattering_channel = tAttr.create("scatteringChannel", "scattering_channel", MFnData::kString);
    tAttr.setDefault(sData.create("density"));
    MPxNode::addAttribute(scattering_channel);

    anisotropy = nAttr.create("anisotropy", "anisotropy", MFnNumericData::kFloat);
    nAttr.setDefault(0.0f);
    nAttr.setMin(-1.0f);
    nAttr.setMax(1.0f);
    MPxNode::addAttribute(anisotropy);

    scattering_gradient.create_params();

    // Emission.

    emission_intensity = nAttr.create("emissionIntensity", "emission_intensity", MFnNumericData::kFloat);
    nAttr.setDefault(1.0f);
    nAttr.setMin(0.0f);
    nAttr.setSoftMax(1.0f);
    nAttr.setChannelBox(true);
    MPxNode::addAttribute(emission_intensity);

    emission_color = nAttr.createColor("emissionColor", "emission_color");
    nAttr.setDefault(1.0, 1.0, 1.0);
    MPxNode::addAttribute(emission_color);

    emission_channel = tAttr.create("emissionChannel", "emission_channel", MFnData::kString);
    tAttr.setDefault(sData.create("heat"));
    MPxNode::addAttribute(emission_channel);

    emission_gradient.create_params();
}

void VDBVolumeStandardShaderParams::affect_output(MObject& out_object)
{
    MPxNode::attributeAffects(density, out_object);
    MPxNode::attributeAffects(density_channel, out_object);
    density_gradient.affect_output(out_object);

    MPxNode::attributeAffects(scattering_intensity, out_object);
    MPxNode::attributeAffects(scattering_color, out_object);
    MPxNode::attributeAffects(scattering_channel, out_object);
    MPxNode::attributeAffects(anisotropy, out_object);
    scattering_gradient.affect_output(out_object);

    MPxNode::attributeAffects(transparency, out_object);
    MPxNode::attributeAffects(transparency_channel, out_object);
    transparency_gradient.affect_output(out_object);

    MPxNode::attributeAffects(emission_intensity, out_object);
    MPxNode::attributeAffects(emission_color, out_object);
    MPxNode::attributeAffects(emission_channel, out_object);
    MPxNode::attributeAffects(emission_mode, out_object);
    emission_gradient.affect_output(out_object);

    MPxNode::attributeAffects(temperature_channel, out_object);
    temperature_gradient.affect_output(out_object);
}

bool VDBVolumeStandardShaderParams::check_plug(const MPlug& plug)
{
    return plug == density ||
           plug == density_channel ||
           density_gradient.check_plug(plug) ||
           plug == scattering_intensity ||
           plug == scattering_color ||
           plug == scattering_channel ||
           plug == anisotropy ||
           scattering_gradient.check_plug(plug) ||
           plug == transparency ||
           plug == transparency_channel ||
           transparency_gradient.check_plug(plug) ||
           plug == emission_intensity ||
           plug == emission_color ||
           plug == emission_channel ||
           plug == emission_mode ||
           emission_gradient.check_plug(plug) ||
           plug == temperature ||
           plug == temperature_channel ||
           temperature_gradient.check_plug(plug);
}
