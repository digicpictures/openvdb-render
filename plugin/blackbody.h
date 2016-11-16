#include <maya/MFloatVector.h>

namespace Blackbody
{
    constexpr float TEMPERATURE_MIN = 800.0f;
    constexpr float TEMPERATURE_MAX = 12000.0f;
    constexpr int TABLE_SIZE = 256;

    extern const float LUT_NORMALIZER;
    extern const MFloatVector LUT[TABLE_SIZE];
};