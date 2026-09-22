#include "BasisFunctionsRadialSpline.hpp"

#include "ParallelEnvironemt.hpp"
#include "constants.hpp"
#include "cutoff.hpp"
#include "debug.hpp"
#include "math.hpp"
#include "utils.hpp"

#include <cmath>
#include <stdexcept>

using namespace vaspml;
using CT = vaspml::math::CutoffType;

BasisFunctionsRadialSpline::BasisFunctionsRadialSpline() : BasisFunctions()
{
    maxOrderBessel = 0;
    nRootsBessel = 0;
    nGrid = 0;
}

BasisFunctionsRadialSpline::BasisFunctionsRadialSpline(CT                cutoffType_in,
                                                       Real              cutoff_in,
                                                       Real              widthBroadening_in,
                                                       std::size_t       nGrid_in,
                                                       std::size_t       maxOrderBessel_in,
                                                       std::size_t       nRootsBessel_in,
                                                       BasisFunctionType type) :
    BasisFunctions(cutoffType_in, cutoff_in, widthBroadening_in, type)
{
    nGrid = nGrid_in;
    gridSpacing = (cutoff) / (Real)nGrid;
    gridSpacing += gridSpacing / (Real)1;
    gridSpacing += gridSpacing / nGrid;
    maxOrderBessel_store = maxOrderBessel_in;
    maxOrderBessel = maxOrderBessel_store + 1;
    nRootsBessel = nRootsBessel_in;
    update();
}

// getters
Real BasisFunctionsRadialSpline::get_gridSpacing() const
{
    return gridSpacing;
}

Int BasisFunctionsRadialSpline::get_nGrid() const
{
    return nGrid;
}

std::size_t BasisFunctionsRadialSpline::get_maxOrder() const
{
    return maxOrderBessel_store;
}

Int BasisFunctionsRadialSpline::get_nRootsBessel() const
{
    return nRootsBessel;
}

std::size_t BasisFunctionsRadialSpline::get_totalNumberBasisFunctions() const
{
    return totalNumberBasisFunctions;
}

std::size_t BasisFunctionsRadialSpline::get_ldim() const
{
    return 1;
}

const std::vector<std::size_t>& BasisFunctionsRadialSpline::get_nValidRoots() const
{
    return nValidRoots;
}
// ~~~~~~~~~~~~~~~~~~~

// setters
void BasisFunctionsRadialSpline::set_nGrid(const std::size_t in)
{
    nGrid = in;

    return;
}

void BasisFunctionsRadialSpline::set_maxOrderBessel(const std::size_t in)
{
    maxOrderBessel = in + 1;
    maxOrderBessel_store = in;

    return;
}

void BasisFunctionsRadialSpline::set_nRootsBessel(const std::size_t in)
{
    nRootsBessel = in;

    return;
}

void BasisFunctionsRadialSpline::set_gridSpacing()
{
    gridSpacing = cutoff / (Real)nGrid;
    gridSpacing += gridSpacing / nGrid;

    return;
}

// ~~~~~~~~~~~~~~~~~~~

void BasisFunctionsRadialSpline::makeGrid()
{

    grid.resize(nGrid);
    for (std::size_t i = 0; i < nGrid; i++)
    {
        grid[i] = gridSpacing * (Real)(i + 1) + 0.5 * gridSpacing;
    }

    return;
}

void BasisFunctionsRadialSpline::printValuesRadial() const
{
    std::cout << "nGrid          = " << nGrid << std::endl;
    std::cout << "maxOrderBessel = " << maxOrderBessel_store << std::endl;
    std::cout << "nRootsBessel   = " << nRootsBessel << std::endl;

    return;
}

void BasisFunctionsRadialSpline::computeCutoffFunction()
{

    cutoffFunction.resize(nGrid);
    cutoffFunction_derivative.resize(nGrid);

    math::cutoffAndDerivative(grid, cutoffFunction, cutoffFunction_derivative, cutoff, cutoffType);

    return;
}

void BasisFunctionsRadialSpline::update()
{

    Vec3Real sphericalBessel;
    set_gridSpacing();
    makeGrid();
    // compute cutoff function for radial basis
    computeCutoffFunction();

    // temporary array to store roots (zeros) of the basis set
    Vec2Real rootsBesselBasis;
    // compute roots of Bessel functions
    computeSphericalBesselRoots(rootsBesselBasis, maxOrderBessel);
    checkSphericalBesselRootsConvergence(rootsBesselBasis);
    // rescaling roots such that Basis set is zero at cutoff
    rescaleRoots(rootsBesselBasis);
    // compute the spherical Bessel functions, scaled by the roots
    computeSphericalBesselFunctions(sphericalBessel, rootsBesselBasis);
    // normalize spherical Bessel functions
    computeBesselNorms(sphericalBessel);

    Vec3Real projectionRadial;
    Vec2Real projectionRadialDerivativeOrigin;
    computeProjectionRadial(sphericalBessel, projectionRadial);
    // compute derivatives of modified spherical bessel functions
    // at grid origin to compute spline interpolation
    computeProjectionRadialDerivativeOrigin(sphericalBessel, projectionRadialDerivativeOrigin);
    computeSplineCoefficients(projectionRadial, projectionRadialDerivativeOrigin);

    return;
}

void BasisFunctionsRadialSpline::checkSphericalBesselRootsConvergence(Vec2Real& roots)
{

    nValidRoots.resize(maxOrderBessel);
    std::fill(nValidRoots.begin(), nValidRoots.end(), 0);
    totalNumberBasisFunctions = 0;
    for (std::size_t actOrder = 0; actOrder < maxOrderBessel; actOrder++)
    {
        for (std::size_t i = 0; i < roots[actOrder].size(); i++)
        {
            if (roots[actOrder][i] <= rootMax)
            {
                nValidRoots[actOrder]++;
                totalNumberBasisFunctions++;
            }
            else { roots[actOrder][i] = (Real)0; }
        }
    }

    return;
}

void BasisFunctionsRadialSpline::computeSphericalBesselRoots(Vec2Real& rootsBesselBasis,
                                                             const Int maxOrderBessel)
{

    vector_tools::allocate_vector(rootsBesselBasis, maxOrderBessel, nRootsBessel);
    rootMax = 1.0e12;
    for (Int actOrder = 0; actOrder < maxOrderBessel; actOrder++)
    {
        math::sphericalBesselRoots(rootsBesselBasis[actOrder], actOrder);
        rootMax = std::min(rootMax, rootsBesselBasis[actOrder][nRootsBessel - 1]);
    }

    return;
}

void BasisFunctionsRadialSpline::rescaleRoots(Vec2Real& rootsBesselBasis)
{
    for (auto& r : rootsBesselBasis)
    {
        for (auto& x : r) { x /= cutoff; }
    }
    rootMax /= cutoff;

    return;
}

void BasisFunctionsRadialSpline::computeSphericalBesselFunctions(Vec3Real&       bessel,
                                                                 const Vec2Real& rootsBesselBasis)
{

    vector_tools::allocate_vector(bessel, maxOrderBessel);
    for (std::size_t actOrder = 0; actOrder < maxOrderBessel; actOrder++)
    {
        computeSphericalBesselFunctionsSingleOrder(bessel[actOrder],
                                                   rootsBesselBasis[actOrder],
                                                   actOrder,
                                                   nValidRoots[actOrder]);
    }

    return;
}

void BasisFunctionsRadialSpline::computeSphericalBesselFunctionsSingleOrder(
    Vec2Real&         bessel,
    const Vec1Real&   roots,
    const std::size_t actOrderBessel,
    const std::size_t nBasis)
{
    vector_tools::allocate_vector(bessel, nBasis, nGrid);
    for (std::size_t i = 0; i < nBasis; i++)
    {
        computeSphericalBesselFunctionsSingleRoot(bessel[i], roots[i], actOrderBessel);
    }

    return;
}

void BasisFunctionsRadialSpline::computeSphericalBesselFunctionsSingleRoot(Vec1Real& bessel,
                                                                           Real      root,
                                                                           const Int actOrderBessel)
{
    for (std::size_t i = 0; i < bessel.size(); i++)
    {
        bessel[i] = math::sphericalBessel(grid[i] * root, actOrderBessel);
    }

    return;
}

void BasisFunctionsRadialSpline::computeBesselNorms(Vec3Real& bessel)
{

    for (std::size_t actOrder = 0; actOrder < maxOrderBessel; actOrder++)
    {
        computeBesselNormsSingleOrder(bessel[actOrder], nValidRoots[actOrder]);
    }

    return;
}

void BasisFunctionsRadialSpline::computeBesselNormsSingleOrder(Vec2Real& bessel, std::size_t nBasis)
{

    Vec1Real normalization;
    normalization.resize(nBasis);
    for (std::size_t i = 0; i < nBasis; i++)
    {
        normalization[i] = math::trapezSphericalL2Integral(bessel[i], grid, gridSpacing);
        normalization[i] = std::sqrt(normalization[i]);
        math::vectorTimesScalarNoCopy(bessel[i], (Real)1 / normalization[i]);
    }

    return;
}

void BasisFunctionsRadialSpline::computeProjectionRadial(const Vec3Real& sphericalBessel,
                                                         Vec3Real&       projectionRadial)
{

    Real prefactor = gridSpacing * constants::PI4 * pow(widthBroadening / constants::PI, (Real)1.5);

    Vec1Real rrj = grid;
    math::squareVector(rrj);
    math::vectorTimesScalarNoCopy(rrj, widthBroadening);

    Vec2Real product;
    vector_tools::allocate_vector(product, maxOrderBessel, nGrid);

    Vec1Real temp_integral(nGrid);

    // allocate modifiedSphericalBessel maxOrderBessel,nRoots[l],nGrid
    vector_tools::allocate_vector(projectionRadial, maxOrderBessel);
    for (std::size_t actOrder = 0; actOrder < maxOrderBessel; actOrder++)
    {
        vector_tools::allocate_vector(projectionRadial[actOrder], nValidRoots[actOrder]);
        for (std::size_t root = 0; root < nValidRoots[actOrder]; root++)
        {
            vector_tools::allocate_vector(projectionRadial[actOrder][root], nGrid);
        }
    }

    for (std::size_t i = 0; i < nGrid; i++)
    {

        Vec1Real rr = math::vectorTimesScalar(grid, 2.0 * widthBroadening * grid[i]);

        Real rri = widthBroadening * grid[i] * grid[i];

        math::expModSphBessel(rri, rr, rrj, product);
        for (std::size_t actOrder = 0; actOrder < maxOrderBessel; actOrder++)
        {
            for (std::size_t j = 0; j < nGrid; j++)
            {
                temp_integral[j] =
                    prefactor * cutoffFunction[i] * product[actOrder][j] * grid[j] * grid[j];
            }
            for (std::size_t root = 0; root < nValidRoots[actOrder]; root++)
            {
                projectionRadial[actOrder][root][i] =
                    math::dotProduct(temp_integral, sphericalBessel[actOrder][root]);
            }
        }
    }

    return;
}

void BasisFunctionsRadialSpline::computeProjectionRadialDerivativeOrigin(
    const Vec3Real& sphericalBessel,
    Vec2Real&       projectionRadialDerivativeOrigin)
{

    // for more details on this routine one should consider the
    // paper https://doi.org/10.1103/PhysRevB.100.014105 equation 14

    Real prefactor = gridSpacing * constants::PI4 * pow(widthBroadening / constants::PI, 1.5);

    Vec1Real gauss;
    Vec1Real gauss_derivative;

    vector_tools::allocate_vector(gauss, nGrid);
    vector_tools::allocate_vector(gauss_derivative, nGrid);

    Vec2Real modSphBessel;
    Vec2Real modSphBesselDerivative;
    vector_tools::allocate_vector(modSphBessel, maxOrderBessel, nGrid);
    vector_tools::allocate_vector(modSphBesselDerivative, maxOrderBessel, nGrid);

    math::gaussianAndDerivative(grid, gauss, gauss_derivative, widthBroadening);

    // rescaled used to compute the spherical Bessel function and it's derivative
    Vec1Real rescaledGrid = math::vectorTimesScalar(grid, (Real)2 * widthBroadening * grid[0]);

    math::modifiedSphericalBesselAndDerivative(grid, modSphBessel, modSphBesselDerivative);

    // precomputed factors needed for computing the integral in Eq 14
    // of paper https://doi.org/10.1103/PhysRevB.100.014105
    Real factor1 = prefactor * cutoffFunction_derivative[0] * gauss[0];
    Real factor2 = prefactor * cutoffFunction[0] * gauss_derivative[0];
    Real factor3 = prefactor * cutoffFunction[0] * gauss[0] * (Real)2 * widthBroadening;
    vector_tools::allocate_vector(projectionRadialDerivativeOrigin, maxOrderBessel);
    for (std::size_t actOrder = 0; actOrder < maxOrderBessel; actOrder++)
    {
        projectionRadialDerivativeOrigin[actOrder].resize(nValidRoots[actOrder]);
        for (std::size_t root = 0; root < nValidRoots[actOrder]; root++)
        {
            Real value = computeProjectionRadialDerivativeOriginSingle(
                factor1,
                factor2,
                factor3,
                gauss,
                sphericalBessel[actOrder][root],
                modSphBessel[actOrder].begin(),
                modSphBesselDerivative[actOrder].begin());
            projectionRadialDerivativeOrigin[actOrder][root] = value;
        }
    }

    return;
}

Real BasisFunctionsRadialSpline::computeProjectionRadialDerivativeOriginSingle(
    const Real               factor1,
    const Real               factor2,
    const Real               factor3,
    const Vec1Real&          gauss,
    const Vec1Real&          spherical_bessel,
    Vec1Real::const_iterator modifiedBessel,
    Vec1Real::const_iterator modifiedBesselDerivative)
{
    Real integral = 0;
    for (std::size_t ir = 0; ir < nGrid; ir++)
    {
        integral += (gauss[ir] * (*(modifiedBessel + ir)) * (factor1 + factor2)
                     + factor3 * gauss[ir] * (*(modifiedBesselDerivative + ir)) * grid[ir])
                  * grid[ir] * grid[ir] * spherical_bessel[ir];
    }

    return integral;
}

// interpolation part
void BasisFunctionsRadialSpline::computeSplineCoefficients(const Vec3Real& functions,
                                                           const Vec2Real& derivative)
{

    vector_tools::allocate_vector(splines, functions.size());
    for (std::size_t actOrder = 0; actOrder < functions.size(); actOrder++)
    {
        vector_tools::allocate_vector(splines[actOrder], nValidRoots[actOrder]);
        for (std::size_t root = 0; root < functions[actOrder].size(); root++)
        {
            splines[actOrder][root].update(grid,
                                           functions[actOrder][root],
                                           derivative[actOrder][root]);
        }
    }

    return;
}

Vec2Real BasisFunctionsRadialSpline::interpolate(const Real x) const
{

    Vec2Real result;
    vector_tools::allocate_vector(result, splines.size());

    for (std::size_t actOrder = 0; actOrder < splines.size(); actOrder++)
    {
        vector_tools::allocate_vector(result[actOrder], nValidRoots[actOrder]);
        for (std::size_t root = 0; root < splines[actOrder].size(); root++)
        {
            result[actOrder][root] = splines[actOrder][root].interpolateValue(x);
        }
    }

    return result;
}

void BasisFunctionsRadialSpline::interpolate(const Real x, Vec2Real& values) const
{

    for (std::size_t actOrder = 0; actOrder < splines.size(); actOrder++)
    {
        for (std::size_t root = 0; root < splines[actOrder].size(); root++)
        {
            values[actOrder][root] = splines[actOrder][root].interpolateValue(x);
        }
    }

    return;
}

void BasisFunctionsRadialSpline::interpolate(const Real x,
                                             Vec2Real&  values,
                                             Vec2Real&  derivatives) const
{

    vector_tools::allocate_vector(values, splines.size());
    vector_tools::allocate_vector(derivatives, splines.size());
    for (std::size_t actOrder = 0; actOrder < splines.size(); actOrder++)
    {
        values[actOrder].resize(splines[actOrder].size());
        derivatives[actOrder].resize(splines[actOrder].size());
        for (std::size_t root = 0; root < splines[actOrder].size(); root++)
        {
            splines[actOrder][root].computeSplineAndDerivative(x,
                                                               values[actOrder][root],
                                                               derivatives[actOrder][root]);
        }
    }

    return;
}
VASPML_NV_HOST_DEVICE
void BasisFunctionsRadialSpline::interpolate(const Real x,
                                             Vec1Real&  values,
                                             Vec1Real&  derivatives) const
{
#ifndef VASPML_NV_HOST_DEVICE
    VASPML_DEBUG_L1(
        if (totalNumberBasisFunctions != values.size())
        {
            throw std::runtime_error("ERROR in BasisFunctionsRadialSpline::interpolate: Size of "
                                     "supplied function array incorrect");
        }
    );
    VASPML_DEBUG_L1(
        if (totalNumberBasisFunctions != derivatives.size())
        {
            throw std::runtime_error(
              "ERROR in BasisFunctionsRadialSpline::interpolate: Size of supplied derivative"
              "function array incorrect");
        }
    );
#endif

    std::size_t col = 0;
    Real        func;
    Real        func_deri;
    for (std::size_t actOrder = 0; actOrder < splines.size(); actOrder++)
    {
        for (std::size_t root = 0; root < splines[actOrder].size(); root++)
        {
            splines[actOrder][root].computeSplineAndDerivative(x, func, func_deri);
            values[col] = func;
            derivatives[col] = func_deri;
            col++;
        }
    }

    return;
}
VASPML_NV_HOST_DEVICE
void BasisFunctionsRadialSpline::computeAngularBasis(const Vec1Real&,
                                                     const Vec1Real& /* norm */,
                                                     Vec1Real& ylm,
                                                     Vec1Real& ylmd) const
{
    const Real radFactor = (Real)1 / std::sqrt((Real)4 * constants::PI);
    std::fill(ylm.begin(), ylm.end(), radFactor);
    std::fill(ylmd.begin(), ylmd.end(), (Real)0);
    return;
}
