#include "BasisFunctions.hpp"
#include "BasisFunctionsAngular.hpp"
#include "BasisFunctionsRadialSpline.hpp"
#include "ParallelEnvironemt.hpp"
#include "SmartEnum.hpp"
#include "Tutor.hpp"
#include "cutoff.hpp"

#include <cmath>
#include <iostream>

using namespace vaspml;

using CT = vaspml::math::CutoffType;

BasisFunctions::BasisFunctions() :
    cutoffType(CT::NONE),
    cutoff(0.0),
    widthBroadening(0.0),
    functionType(BasisFunctionType::none)
{}

BasisFunctions::BasisFunctions(const CT                cutoffType_in,
                               const Real              cutoff_in,
                               const Real              widthBroadening_in,
                               const BasisFunctionType type) :
    cutoffType(cutoffType_in),
    cutoff(cutoff_in),
    functionType(type)
{
    set_widthBroadening(widthBroadening_in);
}

void BasisFunctions::printValues()
{
    std::cout << "cutoffType      = " << cutoffType << std::endl;
    std::cout << "cutoff          = " << cutoff << std::endl;
    std::cout << "widthBroadening = " << get_widthBroadening() << std::endl;

    return;
}

// Setters
void BasisFunctions::set_cutoffType(const CT in)
{
    cutoffType = in;
    return;
}

void BasisFunctions::set_cutoff(const Real in)
{
    cutoff = in;

    return;
}

void BasisFunctions::set_widthBroadening(const Real in)
{
    widthBroadening = 0.5 / (in * in);

    return;
}

// Getters
CT BasisFunctions::get_cutoffType() const
{
    return cutoffType;
}

Real BasisFunctions::get_cutoff() const
{
    return cutoff;
}

Real BasisFunctions::get_widthBroadening() const
{
    return std::sqrt(0.5 / widthBroadening);
}

std::size_t BasisFunctions::get_maxOrder(void) const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        return static_cast<const BasisFunctionsRadialSpline*>(this)->get_maxOrder();
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        return static_cast<const BasisFunctionsAngular*>(this)->get_maxOrder();
    }
    else
    {
        String functionName = __func__;
        global_scope::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
    throw;
}

Real BasisFunctions::get_gridSpacing(void) const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        return static_cast<const BasisFunctionsRadialSpline*>(this)->get_gridSpacing();
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        return static_cast<const BasisFunctionsAngular*>(this)->get_gridSpacing();
    }
    else
    {
        String functionName = __func__;
        global_scope::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
    throw;
}

Int BasisFunctions::get_nGrid(void) const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        return static_cast<const BasisFunctionsRadialSpline*>(this)->get_nGrid();
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        return static_cast<const BasisFunctionsAngular*>(this)->get_nGrid();
    }
    else
    {
        String functionName = __func__;
        global_scope::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
    throw;
}

Int BasisFunctions::get_nRootsBessel(void) const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        return static_cast<const BasisFunctionsRadialSpline*>(this)->get_nRootsBessel();
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        return static_cast<const BasisFunctionsAngular*>(this)->get_nRootsBessel();
    }
    else
    {
        String functionName = __func__;
        global_scope::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
    throw;
}

std::size_t BasisFunctions::get_totalNumberBasisFunctions(void) const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        return static_cast<const BasisFunctionsRadialSpline*>(this)
            ->get_totalNumberBasisFunctions();
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        return static_cast<const BasisFunctionsAngular*>(this)->get_totalNumberBasisFunctions();
    }
    else
    {
        String functionName = __func__;
        global_scope::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
    throw;
}

std::size_t BasisFunctions::get_ldim(void) const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        return static_cast<const BasisFunctionsRadialSpline*>(this)->get_ldim();
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        return static_cast<const BasisFunctionsAngular*>(this)->get_ldim();
    }
    else
    {
        String functionName = __func__;
        global_scope::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
    throw;
}

const std::vector<std::size_t>& BasisFunctions::get_nValidRoots(void) const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        return static_cast<const BasisFunctionsRadialSpline*>(this)->get_nValidRoots();
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        return static_cast<const BasisFunctionsAngular*>(this)->get_nValidRoots();
    }
    else
    {
        String functionName = __func__;
        global_scope::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
    throw;
}

void BasisFunctions::set_nGrid(const std::size_t in)
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        static_cast<BasisFunctionsRadialSpline*>(this)->set_nGrid(in);
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        static_cast<BasisFunctionsAngular*>(this)->set_nGrid(in);
    }
    else
    {
        String functionName = __func__;
        global_scope::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
}

void BasisFunctions::set_maxOrderBessel(const std::size_t in)
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        static_cast<BasisFunctionsRadialSpline*>(this)->set_maxOrderBessel(in);
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        static_cast<BasisFunctionsAngular*>(this)->set_maxOrderBessel(in);
    }
    else
    {
        String functionName = __func__;
        global_scope::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
}

void BasisFunctions::set_nRootsBessel(const std::size_t in)
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        static_cast<BasisFunctionsRadialSpline*>(this)->set_nRootsBessel(in);
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        static_cast<BasisFunctionsAngular*>(this)->set_nRootsBessel(in);
    }
    else
    {
        String functionName = __func__;
        global_scope::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
}

void BasisFunctions::printValuesRadial(void) const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        static_cast<const BasisFunctionsRadialSpline*>(this)->printValuesRadial();
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        static_cast<const BasisFunctionsAngular*>(this)->printValuesRadial();
    }
    else
    {
        String functionName = __func__;
        global_scope::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
}

void BasisFunctions::update(void)
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        return static_cast<BasisFunctionsRadialSpline*>(this)->update();
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        return static_cast<BasisFunctionsAngular*>(this)->update();
    }
    else
    {
        String functionName = __func__;
        global_scope::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
}

Vec2Real BasisFunctions::interpolate(const Real x) const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        static_cast<const BasisFunctionsRadialSpline*>(this)->interpolate(x);
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        static_cast<const BasisFunctionsAngular*>(this)->interpolate(x);
    }
    else
    {
        String functionName = __func__;
        global_scope::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
    throw;
}

void BasisFunctions::interpolate(const Real x, Vec2Real& result) const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        static_cast<const BasisFunctionsRadialSpline*>(this)->interpolate(x, result);
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        static_cast<const BasisFunctionsAngular*>(this)->interpolate(x, result);
    }
    else
    {
        String functionName = __func__;
        global_scope::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
}

void BasisFunctions::interpolate(const Real x, Vec2Real& result, Vec2Real& result_derivative) const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        static_cast<const BasisFunctionsRadialSpline*>(this)->interpolate(x,
                                                                          result,
                                                                          result_derivative);
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        static_cast<const BasisFunctionsAngular*>(this)->interpolate(x, result, result_derivative);
    }
    else
    {
        String functionName = __func__;
        global_scope::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
}

VASPML_NV_HOST_DEVICE
void BasisFunctions::interpolate(const Real x, Vec1Real& result, Vec1Real& result_derivative) const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        static_cast<const BasisFunctionsRadialSpline*>(this)->interpolate(x,
                                                                          result,
                                                                          result_derivative);
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        static_cast<const BasisFunctionsAngular*>(this)->interpolate(x, result, result_derivative);
    }
#ifndef USE_NVIDIA_GPU
    else
    {
        String functionName = __func__;
        global_scope::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
#endif
}

VASPML_NV_HOST_DEVICE
void BasisFunctions::computeAngularBasis(const Vec1Real& rhat,
                                         const Vec1Real& rnorm,
                                         Vec1Real&       ylm,
                                         Vec1Real&       ylmd) const
{
    if (functionType == BasisFunctionType::bodyOrder2)
    {
        static_cast<const BasisFunctionsRadialSpline*>(this)->computeAngularBasis(rhat,
                                                                                  rnorm,
                                                                                  ylm,
                                                                                  ylmd);
    }
    else if (functionType == BasisFunctionType::bodyOrder3)
    {
        static_cast<const BasisFunctionsAngular*>(this)->computeAngularBasis(rhat,
                                                                             rnorm,
                                                                             ylm,
                                                                             ylmd);
    }
#ifndef USE_NVIDIA_GPU
    else
    {
        String functionName = __func__;
        global_scope::tutor.bug("ERROR: " + functionName + " function switch not implemented");
    }
#endif
}
