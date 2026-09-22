#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include "Timer.hpp"
#include "types.hpp"

#include <cmath>

namespace vaspml::constants
{

// Pi variants
constexpr Real PI = 3.14159265358979323846264338327950288419716939937510;
constexpr Real PI2 = 2 * PI;
constexpr Real PI3 = 3 * PI;
constexpr Real PI4 = 4 * PI;
constexpr Real SQRT_PI = 1.77245385090551602729816748334114518279754945612239;

// Euler's number
constexpr Real E = 2.71828182845904523536028747135266249775724709369995;

// Other constants
constexpr Real SQRT2 = 1.41421356237309504880168872420969807856967187537695;
constexpr Real EPS_TOL = 1.0e-10;

/// conversion factor from angstroem to bohr
constexpr Real AUTOA = 0.529177249;
constexpr Real RYTOEV = 13.605826;
constexpr Real EVTOJ = 1.60217733e-19;
constexpr Real AMTOKG = 1.660540e-27;

constexpr Real EUNIT = 2.0 * RYTOEV;
constexpr Real MUNIT =
    1.054571726e-034 * 1.054571726e-034 / 2.0 / (RYTOEV * EVTOJ * AUTOA * AUTOA * 1.0e-20) / AMTOKG;
constexpr Real FUNIT = 2.0 * RYTOEV / AUTOA;
constexpr Real SUNIT =
    2.0 * RYTOEV / ((AUTOA * 1.0e-10) * (AUTOA * 1.0e-10) * (AUTOA * 1.0e-10)) * EVTOJ / 1.0e8;

// key list which is needed in Descriptor collector and Kernel routines
// for extracting needed descriptors from descriptor map
const Vec1String descriptorKeyList = {"SHS2-2-body", "SHS3-3-body"};

} // namespace vaspml::constants

#endif
