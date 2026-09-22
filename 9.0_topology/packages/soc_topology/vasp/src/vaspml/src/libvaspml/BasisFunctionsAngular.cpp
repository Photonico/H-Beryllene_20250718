#include "BasisFunctionsAngular.hpp"

#include "ParallelEnvironemt.hpp"
#include "SphericalHarmonics.hpp"
#include "constants.hpp"

#include <cmath>

using namespace vaspml;
using CT = vaspml::math::CutoffType;

BasisFunctionsAngular::BasisFunctionsAngular() : BasisFunctionsRadialSpline()
{
    lmax = 0;
}

BasisFunctionsAngular::BasisFunctionsAngular(CT                cutOffType_in,
                                             Real              cutOff_in,
                                             Real              widthBroadening_in,
                                             std::size_t       nGrid_in,
                                             std::size_t       maxOrderBessel_in,
                                             std::size_t       nRootsBessel_in,
                                             BasisFunctionType type,
                                             std::size_t       angularFiltering_in,
                                             Real              filterScale_in) :
    BasisFunctionsRadialSpline(cutOffType_in,
                               cutOff_in,
                               widthBroadening_in,
                               nGrid_in,
                               maxOrderBessel_in,
                               nRootsBessel_in,
                               type)
{
    lmax = maxOrderBessel_in;
    angularFiltering = angularFiltering_in;
    filterScale = filterScale_in;
    angularBasis = std::make_shared<math::SphericalHarmonics>(lmax);
    ldim = angularBasis->get_ldim();
    computeAngularFilteringCoefficient();
}

void BasisFunctionsAngular::set_lmax(std::size_t lmax_in)
{

    lmax = lmax_in;
    angularBasis = std::make_shared<math::SphericalHarmonics>(lmax);
    ldim = angularBasis->get_ldim();
}

VASPML_NV_HOST_DEVICE
void BasisFunctionsAngular::computeAngularBasis(const Vec1Real& normed_xyz,
                                                const Vec1Real& norm,
                                                Vec1Real&       ylm,
                                                Vec1Real&       ylmd) const
{

    angularBasis->computeSphericalHarmonicsAndGradients(normed_xyz, norm, ylm, ylmd);
}

void BasisFunctionsAngular::computeAngularFilteringCoefficient()
{

    Real pi8 = (Real)2 * constants::PI2 * constants::PI2;
    filteringFactors.resize(lmax + 1);
    if (angularFiltering == 1)
    {
        for (std::size_t ll = 0; ll < lmax + 1; ll++)
        {
            filteringFactors[ll] = std::sqrt(pi8 / (Real)(2 * ll + 1)) / (Real)(2 * ll + 1);
        }
    }
    else if (angularFiltering == 2)
    {
        for (std::size_t ll = 0; ll < lmax + 1; ll++)
        {
            Real factor = (Real)(1 + filterScale * ll * (ll + 1) * ll * (ll + 1));
            factor *= factor;
            filteringFactors[ll] = std::sqrt(pi8 / (Real)(2 * ll + 1)) / factor;
        }
    }
    else
    {
        for (std::size_t ll = 0; ll < lmax + 1; ll++)
        {
            filteringFactors[ll] = std::sqrt(pi8 / (Real)(2 * ll + 1));
        }
    }
    return;
}

const Vec1Real& BasisFunctionsAngular::get_filteringFactors() const
{
    return filteringFactors;
}

std::size_t BasisFunctionsAngular::get_ldim() const
{
    return ldim;
}
