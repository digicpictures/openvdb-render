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
extern double cieMatch[81][3];

/// Function to convert a spectral intensity to CIE XYZ. The template argument SpecFunc should be a functor that
/// descrbes the desired spectral function. The functor should have the member:
/// float operator(const float wavelength)
/// which should return the spectral intensity of the function at the given wavelength
template <class SpecFunc>
AtColor spectrumToXyz(const SpecFunc& sf)
{
    double lambda = 380.0f;
    double X=0.0f, Y=0.0f, Z=0.0f;
    double lambda_step = 5;
    double weight = lambda_step * 1e-9;
    for (int i = 0; lambda < 780.1; i++, lambda += lambda_step)
    {
        double Me;

        Me = sf(lambda);
        X += Me * cieMatch[i][0] * weight;
        Y += Me * cieMatch[i][1] * weight;
        Z += Me * cieMatch[i][2] * weight;
    }

    return { float(X), float(Y), float(Z) };
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

double cieMatch[81][3] =
{
    {0.0014,0.0000,0.0065}, {0.0022,0.0001,0.0105}, {0.0042,0.0001,0.0201},
    {0.0076,0.0002,0.0362}, {0.0143,0.0004,0.0679}, {0.0232,0.0006,0.1102},
    {0.0435,0.0012,0.2074}, {0.0776,0.0022,0.3713}, {0.1344,0.0040,0.6456},
    {0.2148,0.0073,1.0391}, {0.2839,0.0116,1.3856}, {0.3285,0.0168,1.6230},
    {0.3483,0.0230,1.7471}, {0.3481,0.0298,1.7826}, {0.3362,0.0380,1.7721},
    {0.3187,0.0480,1.7441}, {0.2908,0.0600,1.6692}, {0.2511,0.0739,1.5281},
    {0.1954,0.0910,1.2876}, {0.1421,0.1126,1.0419}, {0.0956,0.1390,0.8130},
    {0.0580,0.1693,0.6162}, {0.0320,0.2080,0.4652}, {0.0147,0.2586,0.3533},
    {0.0049,0.3230,0.2720}, {0.0024,0.4073,0.2123}, {0.0093,0.5030,0.1582},
    {0.0291,0.6082,0.1117}, {0.0633,0.7100,0.0782}, {0.1096,0.7932,0.0573},
    {0.1655,0.8620,0.0422}, {0.2257,0.9149,0.0298}, {0.2904,0.9540,0.0203},
    {0.3597,0.9803,0.0134}, {0.4334,0.9950,0.0087}, {0.5121,1.0000,0.0057},
    {0.5945,0.9950,0.0039}, {0.6784,0.9786,0.0027}, {0.7621,0.9520,0.0021},
    {0.8425,0.9154,0.0018}, {0.9163,0.8700,0.0017}, {0.9786,0.8163,0.0014},
    {1.0263,0.7570,0.0011}, {1.0567,0.6949,0.0010}, {1.0622,0.6310,0.0008},
    {1.0456,0.5668,0.0006}, {1.0026,0.5030,0.0003}, {0.9384,0.4412,0.0002},
    {0.8544,0.3810,0.0002}, {0.7514,0.3210,0.0001}, {0.6424,0.2650,0.0000},
    {0.5419,0.2170,0.0000}, {0.4479,0.1750,0.0000}, {0.3608,0.1382,0.0000},
    {0.2835,0.1070,0.0000}, {0.2187,0.0816,0.0000}, {0.1649,0.0610,0.0000},
    {0.1212,0.0446,0.0000}, {0.0874,0.0320,0.0000}, {0.0636,0.0232,0.0000},
    {0.0468,0.0170,0.0000}, {0.0329,0.0119,0.0000}, {0.0227,0.0082,0.0000},
    {0.0158,0.0057,0.0000}, {0.0114,0.0041,0.0000}, {0.0081,0.0029,0.0000},
    {0.0058,0.0021,0.0000}, {0.0041,0.0015,0.0000}, {0.0029,0.0010,0.0000},
    {0.0020,0.0007,0.0000}, {0.0014,0.0005,0.0000}, {0.0010,0.0004,0.0000},
    {0.0007,0.0002,0.0000}, {0.0005,0.0002,0.0000}, {0.0003,0.0001,0.0000},
    {0.0002,0.0001,0.0000}, {0.0002,0.0001,0.0000}, {0.0001,0.0000,0.0000},
    {0.0001,0.0000,0.0000}, {0.0001,0.0000,0.0000}, {0.0000,0.0000,0.0000}
};

AtRGB xyzToRgb(const ColorSystem& cs, const AtColor& xyz)
{
    double xr, yr, zr, xg, yg, zg, xb, yb, zb;
    double xw, yw, zw;
    double rx, ry, rz, gx, gy, gz, bx, by, bz;
    double rw, gw, bw;
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

AtRGB xyzToLinearsRGB(const AtColor& xyz)
{
    AtRGB rgb;
    rgb.r =  3.2406 * xyz.r - 1.5372 * xyz.g - 0.4986 * xyz.b;
    rgb.g = -0.9689 * xyz.r + 1.8758 * xyz.g + 0.0415 * xyz.b;
    rgb.b =  0.0557 * xyz.r - 0.2040 * xyz.g + 1.0570 * xyz.b;
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
        return float((1.19031e-16 * pow(lambda, -5.0)) / (exp(1.4388e-2 / (lambda * temp)) - 1.0));
    }

    double temp;
};

MFloatVector blackbodyColorRGB(const float temperature)
{
    const AtColor xyz = spectrumToXyz(BlackbodySpectrum(temperature));
    AtRGB rgb = xyzToLinearsRGB(xyz);
    rgb.r = std::max(0.f, rgb.r);
    rgb.g = std::max(0.f, rgb.g);
    rgb.b = std::max(0.f, rgb.b);

    return { rgb.r, rgb.g, rgb.b };
}