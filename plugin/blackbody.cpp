#include "blackbody.h"

#include <cmath>
#include <maya/MColor.h>

typedef MFloatVector AtPoint2;
typedef MColor AtColor;
typedef MColor AtRGB;

// Copy-pasted from Arnold shaders.

/// A structure representing a color coordinate system
struct ColorSystem
{
    ColorSystem(const std::string& name_, const AtPoint2& r, const AtPoint2& g, const AtPoint2& b, const AtPoint2& w):
        name(name_), red(r), green(g), blue(b), white(w)
    {}

    std::string name;   /// name of the system, e.g. rec709, ACES etc
    AtPoint2 red;       /// chromaticity of the red primary
    AtPoint2 green;     /// chromaticity of the green primary
    AtPoint2 blue;      /// chromticity of the blue primary
    AtPoint2 white;     /// chromaticity of the white point
};

/// White points for common illuminants
extern  AtPoint2 illumC;
extern  AtPoint2 illumD65;
extern  AtPoint2 illumE;

/// Color system definitions for common color spaces
extern  ColorSystem CsNTSC;
extern  ColorSystem CsSMPTE;
extern  ColorSystem CsHDTV;
extern  ColorSystem CsCIE;
extern  ColorSystem CsRec709;

/// CIE spectral color matching functions from 380 to 780 nm, in 5nm increments
extern float cieMatch[81][3];

/// Function to convert a spectral intensity to CIE XYZ. The template argument SpecFunc should be a functor that
/// descrbes the desired spectral function. The functor should have the member:
/// float operator(const float wavelength)
/// which should return the spectral intensity of the function at the given wavelength
template <class SpecFunc>
AtColor spectrumToXyz(const SpecFunc& sf)
{
    AtColor result;
    float lambda = 380.0f;
    float X=0.0f, Y=0.0f, Z=0.0f;
    for (int i = 0; lambda < 780.1; i++, lambda += 5)
    {
        float Me;

        Me = sf(lambda);
        X += Me * cieMatch[i][0];
        Y += Me * cieMatch[i][1];
        Z += Me * cieMatch[i][2];
    }
    float XYZ = (X + Y + Z);

    result.r = X / XYZ ;
    result.g = Y / XYZ ;
    result.b = Z / XYZ ;

    return result;
}

inline AtPoint2 point2(float x, float y)
{
    AtPoint2 r;
    r.x = x;
    r.y = y;
    return r;
}

AtPoint2 illumD65(point2(0.3127f, 0.3291f));
ColorSystem CsRec709("Rec709", point2(0.64f, 0.33f), point2(0.30f, 0.60f), point2(0.15f, 0.06f), illumD65);

float cieMatch[81][3] =
{
    {0.0014f,0.0000f,0.0065f}, {0.0022f,0.0001f,0.0105f}, {0.0042f,0.0001f,0.0201f},
    {0.0076f,0.0002f,0.0362f}, {0.0143f,0.0004f,0.0679f}, {0.0232f,0.0006f,0.1102f},
    {0.0435f,0.0012f,0.2074f}, {0.0776f,0.0022f,0.3713f}, {0.1344f,0.0040f,0.6456f},
    {0.2148f,0.0073f,1.0391f}, {0.2839f,0.0116f,1.3856f}, {0.3285f,0.0168f,1.6230f},
    {0.3483f,0.0230f,1.7471f}, {0.3481f,0.0298f,1.7826f}, {0.3362f,0.0380f,1.7721f},
    {0.3187f,0.0480f,1.7441f}, {0.2908f,0.0600f,1.6692f}, {0.2511f,0.0739f,1.5281f},
    {0.1954f,0.0910f,1.2876f}, {0.1421f,0.1126f,1.0419f}, {0.0956f,0.1390f,0.8130f},
    {0.0580f,0.1693f,0.6162f}, {0.0320f,0.2080f,0.4652f}, {0.0147f,0.2586f,0.3533f},
    {0.0049f,0.3230f,0.2720f}, {0.0024f,0.4073f,0.2123f}, {0.0093f,0.5030f,0.1582f},
    {0.0291f,0.6082f,0.1117f}, {0.0633f,0.7100f,0.0782f}, {0.1096f,0.7932f,0.0573f},
    {0.1655f,0.8620f,0.0422f}, {0.2257f,0.9149f,0.0298f}, {0.2904f,0.9540f,0.0203f},
    {0.3597f,0.9803f,0.0134f}, {0.4334f,0.9950f,0.0087f}, {0.5121f,1.0000f,0.0057f},
    {0.5945f,0.9950f,0.0039f}, {0.6784f,0.9786f,0.0027f}, {0.7621f,0.9520f,0.0021f},
    {0.8425f,0.9154f,0.0018f}, {0.9163f,0.8700f,0.0017f}, {0.9786f,0.8163f,0.0014f},
    {1.0263f,0.7570f,0.0011f}, {1.0567f,0.6949f,0.0010f}, {1.0622f,0.6310f,0.0008f},
    {1.0456f,0.5668f,0.0006f}, {1.0026f,0.5030f,0.0003f}, {0.9384f,0.4412f,0.0002f},
    {0.8544f,0.3810f,0.0002f}, {0.7514f,0.3210f,0.0001f}, {0.6424f,0.2650f,0.0000f},
    {0.5419f,0.2170f,0.0000f}, {0.4479f,0.1750f,0.0000f}, {0.3608f,0.1382f,0.0000f},
    {0.2835f,0.1070f,0.0000f}, {0.2187f,0.0816f,0.0000f}, {0.1649f,0.0610f,0.0000f},
    {0.1212f,0.0446f,0.0000f}, {0.0874f,0.0320f,0.0000f}, {0.0636f,0.0232f,0.0000f},
    {0.0468f,0.0170f,0.0000f}, {0.0329f,0.0119f,0.0000f}, {0.0227f,0.0082f,0.0000f},
    {0.0158f,0.0057f,0.0000f}, {0.0114f,0.0041f,0.0000f}, {0.0081f,0.0029f,0.0000f},
    {0.0058f,0.0021f,0.0000f}, {0.0041f,0.0015f,0.0000f}, {0.0029f,0.0010f,0.0000f},
    {0.0020f,0.0007f,0.0000f}, {0.0014f,0.0005f,0.0000f}, {0.0010f,0.0004f,0.0000f},
    {0.0007f,0.0002f,0.0000f}, {0.0005f,0.0002f,0.0000f}, {0.0003f,0.0001f,0.0000f},
    {0.0002f,0.0001f,0.0000f}, {0.0002f,0.0001f,0.0000f}, {0.0001f,0.0000f,0.0000f},
    {0.0001f,0.0000f,0.0000f}, {0.0001f,0.0000f,0.0000f}, {0.0000f,0.0000f,0.0000f}
};

AtRGB xyzToRgb(const ColorSystem& cs, const AtColor& xyz)
{
    float xr, yr, zr, xg, yg, zg, xb, yb, zb;
    float xw, yw, zw;
    float rx, ry, rz, gx, gy, gz, bx, by, bz;
    float rw, gw, bw;
    AtRGB rgb;

    xr = cs.red.x;    yr = cs.red.y;    zr = 1 - (xr + yr);
    xg = cs.green.x;  yg = cs.green.y;  zg = 1 - (xg + yg);
    xb = cs.blue.x;   yb = cs.blue.y;   zb = 1 - (xb + yb);

    xw = cs.white.x;  yw = cs.white.y;  zw = 1 - (xw + yw);

    // xyz -> rgb matrix, before scaling to white
    rx = (yg * zb) - (yb * zg);  ry = (xb * zg) - (xg * zb);  rz = (xg * yb) - (xb * yg);
    gx = (yb * zr) - (yr * zb);  gy = (xr * zb) - (xb * zr);  gz = (xb * yr) - (xr * yb);
    bx = (yr * zg) - (yg * zr);  by = (xg * zr) - (xr * zg);  bz = (xr * yg) - (xg * yr);

    // White scaling factors.
    // Dividing by yw scales the white luminance to unity, as conventional
    rw = ((rx * xw) + (ry * yw) + (rz * zw)) / yw;
    gw = ((gx * xw) + (gy * yw) + (gz * zw)) / yw;
    bw = ((bx * xw) + (by * yw) + (bz * zw)) / yw;

    // xyz -> rgb matrix, correctly scaled to white
    rx = rx / rw;  ry = ry / rw;  rz = rz / rw;
    gx = gx / gw;  gy = gy / gw;  gz = gz / gw;
    bx = bx / bw;  by = by / bw;  bz = bz / bw;

    // rgb of the desired point
    rgb.r = (rx * xyz.r) + (ry * xyz.g) + (rz * xyz.b);
    rgb.g = (gx * xyz.r) + (gy * xyz.g) + (gz * xyz.b);
    rgb.b = (bx * xyz.r) + (by * xyz.g) + (bz * xyz.b);

    return rgb;
}

struct BlackbodySpectrum
{
    BlackbodySpectrum(const float temperature)
        : temp(temperature)
    {}

    float operator()(float wavelength) const
    {
        double lambda = wavelength * 1e-9;
        return float((3.74183e-16 * pow(lambda, -5.0)) / (exp(1.4388e-2 / (lambda * temp)) - 1.0));
        // 3.74183e-16 * pow(lambda, -5.0) == 3.74183e-16 * pow(wavelength, -5.0) * 1e45    <--- TODO
    }

    double temp;
};

MFloatVector blackbodyColorRGB(const float temperature)
{
    const AtColor xyz = spectrumToXyz(BlackbodySpectrum(temperature));
    AtRGB rgb = xyzToRgb(CsRec709, xyz);
    rgb.r = std::max(0.f, rgb.r);
    rgb.g = std::max(0.f, rgb.g);
    rgb.b = std::max(0.f, rgb.b);
    return { rgb.r, rgb.g, rgb.b };
}