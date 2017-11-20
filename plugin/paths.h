#pragma once

#include <string>
#include <maya/MString.h>

class Paths
{
public:
    Paths() = delete;
    Paths(const Paths&) = delete;
    Paths(const Paths&&) = delete;
    Paths& operator=(const Paths&) = delete;
    Paths& operator=(const Paths&&) = delete;

    static void init(const std::string& plugin_path) { s_plugin_path = plugin_path; }
    static void init(const MString& plugin_path) { s_plugin_path = plugin_path.asChar(); }
    static std::string getVolumeEffectFile() { return s_plugin_path + '/' + SHADER_DIR + '/' + VOLUME_SHADER_FILE; }

private:
    static constexpr const char* SHADER_DIR = "shader";
#if MAYA_API_VERSION < 201600
    static constexpr const char* VOLUME_SHADER_FILE = "volume.cgfx";
#else
    static constexpr const char* VOLUME_SHADER_FILE = "volume.ogsfx";
#endif

    static std::string s_plugin_path;
};
