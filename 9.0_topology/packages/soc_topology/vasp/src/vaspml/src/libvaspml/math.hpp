#ifndef MATH_HPP
#define MATH_HPP

#include "Tutor.hpp"
#include "debug.hpp"
#include "types.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <type_traits>

#include <iostream>

namespace vaspml::math
{

/** Fills the given output vector with factorial numbers.
 *
 * @param f Factorial numbers output vector (determines maximum integer, see note).
 *
 * @note The size of the output vector on entry determines the maximum integer
 * number to which factorials should be computed. For factorials up to
 * @f$N_\max@f$ reserve space for one more entry in the output vector, i.e.
 * f.size() = @f$N_\max + 1@f$. The first entry is 0! = 1.
 */
void factorial(Vec1Real& f);
/** Fills the given output vector with double factorials for odd numbers.
 *
 * @param f Double factorial numbers output vector (determines maximum integer, see note).
 *
 * @note The size of the output vector on entry determines the maximum integer number to
 * which double factorials should be computed. To map between the output
 * vector index i and N use the expressions @f$N = 2i + 1@f$ and @f$i = (N - 1)/2@f$.
 * Hence, the size of the output vector should be f.size() = @f$(N_\max + 1) / 2@f$.
 *
 * Example: get N!! up to N = 7 requires a vector with size 4:
 *
 * * i = 0 => N = 1 => N!! = 1
 * * i = 1 => N = 3 => N!! = 3
 * * i = 2 => N = 5 => N!! = 15
 * * i = 3 => N = 7 => N!! = 105
 */
void doubleFactorialOdd(Vec1Real& f);
/** Computes the Clebsch-Gordan coefficient @f$\left<j_1 m_1 j_2 m_2|J M\right>@f$ using Racah
 * recursion relations.
 *
 * @param j1 Index @f$j_1@f$.
 * @param j2 Index @f$j_2@f$.
 * @param J Index @f$J@f$.
 * @param m1 Index @f$m_1@f$.
 * @param m2 Index @f$m_2@f$.
 * @param M Index @f$M@f$.
 * @param factorials Precomputed factorial numbers.
 * @return Clebsch-Gordan coefficient @f$\left<j_1 m_1 j_2 m_2|J M\right>@f$.
 *
 * @note Precomputed factorial numbers need to be passed to this function, use factorial() for
 * this.
 *
 * @todo Find out how many factorials we actually need, default in Fortran code up to 40! makes
 * little sense. I think 4*jmax + 2 should be sufficient. Applies to clebschGordanM0() as well.
 * Alternatively or additionally add debugging checks.
 */
Real clebschGordan(Int j1, Int j2, Int J, Int m1, Int m2, Int M, const Vec1Real& factorials);
/** Computes the Clebsch-Gordan coefficient @f$\left<j_1 0 j_2 0|J 0\right>@f$ using Racah
 * recursion relations.
 *
 * @param j1 Index @f$j_1@f$.
 * @param j2 Index @f$j_2@f$.
 * @param J Index @f$J@f$.
 * @param factorials Precomputed factorial numbers.
 * @return Clebsch-Gordan coefficient @f$\left<j_1 0 j_2 0|J 0\right>@f$.
 *
 * @note Precomputed factorial numbers need to be passed to this function, use factorial() for
 * this.
 */
Real clebschGordanM0(Int j1, Int j2, Int J, const Vec1Real& factorials);
/** Computes Gaussian function for a vector of numbers.
 *
 * @param x x-values for which Gaussian should be computed.
 * @param f Gaussian values for given x-values.
 * @param width Gaussian width parameter.
 *
 * Gaussian is assumed to be centered around zero.
 *
 * @warning Size of output vector must be equal to input vector on entry.
 */
void gaussian(const Vec1Real& x, Vec1Real& f, Real width = 1.0);
/** Computes Gaussian function and its derivative for a vector of numbers.
 *
 * @param x x-values for which Gaussian should be computed.
 * @param f Gaussian values for given x-values.
 * @param df Gaussian derivative for given x-values.
 * @param width Gaussian width parameter.
 *
 * Gaussian is assumed to be centered around zero.
 *
 * @warning Size of output vectors must be equal to input vector on entry.
 */
void gaussianAndDerivative(const Vec1Real& x, Vec1Real& f, Vec1Real& df, Real width = 1.0);
/** Iteratively Computes spherical Bessel function of the first kind for one given value.
 *
 * @param x x-value for which the Bessel function should be computed.
 * @param nu Order of the Bessel function.
 * @return Bessel function values for given x-value.
 *
 * @note This function gets automatically called by #sphericalBessel if the argument is below 1.0.
 */
Real sphericalBesselIterative(const Real x, const UInt nu);
/** Computes spherical Bessel function of the first kind for one given value.
 *
 * @param x x-value for which the Bessel function should be computed.
 * @param nu Order of the Bessel function.
 * @return Bessel function values for given x-value.
 */
Real sphericalBessel(const Real x, const UInt nu);
/** Computes spherical Bessel function of the first kind and its derivative for one given value.
 *
 * @param x x-value for which the Bessel function should be computed.
 * @param jnu Output Bessel function value.
 * @param djnu Output Bessel function derivative.
 * @param nu Order of the Bessel function.
 * @return Bessel function values for given x-value.
 */
void sphericalBesselAndDerivative(const Real x, Real& jnu, Real& djnu, const UInt nu);
/** Computes spherical Bessel function of the first kind for a vector of numbers.
 *
 * @param x x-values for which the Bessel function should be computed.
 * @param jnu Bessel function values for given x-values.
 * @param nu Order of the Bessel function.
 *
 * @warning Size of output vector must be equal to input vector on entry.
 */
void sphericalBessel(const Vec1Real& x, Vec1Real& jnu, const UInt nu);
/** Computes spherical Bessel function of the first kind and its derivative for a vector of numbers.
 *
 * @param x x-values for which the Bessel function should be computed.
 * @param jnu Bessel function values for given x-values.
 * @param djnu Bessel function derivatives for given x-values.
 * @param nu Order of the Bessel function.
 *
 * @warning Size of output vector must be equal to input vector on entry.
 */
void sphericalBesselAndDerivative(const Vec1Real& x, Vec1Real& jnu, Vec1Real& djnu, const UInt nu);
/** Computes multiple roots of a spherical Bessel function of the first kind.
 *
 * @param roots x-values for which the Bessel function is zero, in ascending order (see note).
 * @param nu Order of the Bessel function.
 *
 * @note The size of the output vector on entry to this function determines
 * how many roots are calculated.
 */
void sphericalBesselRoots(Vec1Real& roots, const UInt nu);
/** Computes modified spherical Bessel function of the first kind and its derivative.
 *
 * @warning The output vectors' first dimensions determine the maximum order of the Bessel
 * functions. Hence, on entry the condition `inu.size()` = @f$l_\max + 1@f$ must be satisfied. The
 * second dimension must match with the size of the input vector. In summary,
 * * index 1: Bessel function orders @f$l = 0,\ldots,l_\max@f$.
 * * index 2: Corresponds to input vector x-values.
 *
 * @param x x-values for which the Bessel function should be computed.
 * @param inu Bessel function values for given x-values.
 * @param dinu Bessel function derivatives for given x-values.
 */
void modifiedSphericalBesselAndDerivative(const Vec1Real& x, Vec2Real& inu, Vec2Real& dinu);
/** Computes product of exponential and modified spherical Bessel function of the first kind.
 *
 * @param a2 Input scalar @f$a^2@f$ value.
 * @param ab Input vector @f$ab@f$, i.e., should be a product of scalar @f$a@f$ and vector @f$b@f$.
 * @param b2 Input vector @f$b^2@f$.
 * @param p Output product of the exponential and Bessel function expression.
 *
 * This function is used to compute parts of @f$h_{nl}@f$ found in expression
 * (19) in Phys. Rev. B 100, 014105 (see also https://arxiv.org/pdf/1904.12961.pdf).
 * Given the entire expression
 * @f[
 * h_{nl}\left(r\right) = \frac{4\pi}{\left(\sqrt{2\sigma_\mathrm{atom}^{2}\pi}\right)^{3}}
 * f_\mathrm{cut}\left(r\right)
 * \int  _{0}^{\infty} \chi_{nl}\left(r'\right)
 * \mathrm{exp}\left(-\frac{r^{2}+r'^{2}}{2\sigma_\mathrm{atom}^{2}}\right)
 * \iota_{l}\left(\frac{rr'}{\sigma_{\mathrm  {atom}}^{2}}\right)
 * r'^{2}dr'
 * @f]
 * this function computes the product of the exponential and modified spherical
 * Bessel function (first kind), i.e.,
 * @f[
 * \mathrm{exp}\left(-a^2-b^2\right) \iota_{l}\left(ab\right),
 * @f]
 * where @f$a^2@f$ is passed as a scalar value and @f$ab@f$ and @f$b^2@f$ are passed as vectors of
 * values on a grid. The function computes the product above for all points in the vectors @f$ab@f$
 * and @f$b^2@f$ and for all Bessel function orders up to a given maximum integer number.
 *
 * The reason for computing this product instead of the individual terms is that the modified
 * spherical Bessel function alone quickly diverges for large arguments. However, it can be
 * expressed in terms of hyperbolic sine and cosine, i.e., exponential functions. Multiplying with
 * the exponential function @f$\mathrm{exp}\left(-a^2-b^2\right)@f$ results in terms of
 * @f$\mathrm{exp}\left(-a^2-b^2 \pm ab\right)@f$ which allows the arguments to cancel before the
 * exponentiation is carried out.
 *
 * Expressions for small @f$ab@f$ values are derived from the limit of the modified spherical Bessel
 * function for small arguments:
 * @f[
 * \iota_{l}(x) \approx \frac{x^l}{(2l + 1)!!}
 * @f]
 *
 * @warning The output vector's first dimension determines the maximum order of the Bessel
 * functions. Hence, on entry the condition `p.size()` = @f$l_\max + 1@f$ must be satisfied. The
 * second dimension must match with the size of the input vectors. In summary,
 * * index 1: Bessel function orders @f$l = 0,\ldots,l_\max@f$.
 * * index 2: Corresponds to @f$ab@f$ and @f$b^2@f$ dimension.
 */
void expModSphBessel(const Real a2, const Vec1Real& ab, const Vec1Real& b2, Vec2Real& p);

/** wrapper compute modified spherical Bessel function for max order = 0
 * @param Int  n   maximal order of the bessel functions to compute
 * @param Int  n   number of grid points
 * @param Real r   array which contains the grid and is of size n
 * @param Real ri  scalar which stores an offset
 * @param Real rj  arrray containing a second grid, (Jonathan thinks this is redundant)
 * @param Real ril contains Bessel function on return
 */
void modifiedSphericalBesselInterface(const Int                n,
                                      const std::vector<Real>& r,
                                      const Real               ri,
                                      const std::vector<Real>& rj,
                                      std::vector<Real>&       ril);

/** wrapper compute modified spherical Bessel function for max order > 0
 * @param Int lmax maximal order of the bessel functions to compute
 * @param Int  n   number of grid points
 * @param Real r   array which contains the grid and is of size n
 * @param Real ri  scalar which stores an offset
 * @param Real rj  arrray containing a second grid, (Jonathan thinks this is redundant)
 * @param Real ril contains Bessel function on return
 */
void modifiedSphericalBesselInterface(const Int                       lmax,
                                      const Int                       n,
                                      const std::vector<Real>&        r,
                                      const Real                      ri,
                                      const std::vector<Real>&        rj,
                                      std::vector<std::vector<Real>>& ril);

/** compute spherical Bessel function and derivative for max order = 0
 * should be used from modifiedSphericalBesselInterface in math_c
 * @param Real derivative contains derivative of the function on output
 * @param Int n           number of grid points
 * @param Real r          array which contains the grid and is of size n
 * @param Real function   contains Bessel function on return
 */
void modifiedSphericalBesselDerivative(std::vector<Real>&       derivative,
                                       const Int                n,
                                       const std::vector<Real>& r,
                                       std::vector<Real>&       function);

/** compute integral over std::vector with trapez, and uniform
 * @param func function to integrate
 * @param dx   grid spacing
 */
template<typename T>
T trapezIntegral(const std::vector<Real>& func, const Real dx)
{
    Real integral = (Real)0.5 * func[0];
    integral += (Real)0.5 * func[func.size() - 1];
    integral += std::accumulate(func.begin() + 1, func.end() - 1, 0);
    return integral * dx;
}

/** compute integral of squared function on the spherical radial grid x
 * @param func function to integrate
 * @param x    radial grid hast to be uniform
 * @param dx   grid spacing
 */
template<typename T>
T trapezSphericalL1Integral(const std::vector<T>& func, const std::vector<T>& x, const Real dx)
{
    assert((func.size() == x.size())
           && "trapezSphericalL1Integral: dimension of function and xaxis do not match");
    Real tmp = x[0] * x[0];
    Real integral = (Real)0.5 * func[0] * tmp;
    tmp = x[func.size() - 1] * x[func.size() - 1];
    integral += (Real)0.5 * func[func.size() - 1] * tmp;
    for (auto i = 1; i < func.size() - 1; i++)
    {
        tmp = x[i] * x[i];
        integral += func[i] * tmp;
    }
    return integral * dx;
}

/** compute integral of squared function on the spherical radial grid x
 * @param func function to integrate
 * @param x    radial grid hast to be uniform
 * @param dx   grid spacing
 */
template<typename T>
T trapezSphericalL2Integral(const std::vector<T>& func, const std::vector<T>& x, const Real dx)
{
    assert((func.size() == x.size())
           && "trapezSphericalL2Integral: dimension of function and xaxis do not match");
    Real tmp = func[0] * x[0];
    Real integral = (Real)0.5 * tmp * tmp;
    tmp = func[func.size() - 1] * x[func.size() - 1];
    integral += (Real)0.5 * tmp * tmp;
    for (size_t i = 1; i < func.size() - 1; i++)
    {
        tmp = func[i] * x[i];
        integral += tmp * tmp;
    }
    return integral * dx;
}
/** multiply vector times scalar, note a copy of the vector is made
 *
 * @param Real vec    input vector to multiply
 * @param Real scalar scalar to multiply vec with
 * @param returned will be the multiplied vector
 *
 */
template<typename T>
std::vector<T> vectorTimesScalar(const std::vector<T>& vec, const T scalar)
{
    std::vector<T> result = vec;

    std::transform(result.begin(),
                   result.end(),
                   result.begin(),
                   std::bind(std::multiplies<T>(), std::placeholders::_1, scalar));
    return result;
}
/** multiply vector times scalar, no copy of the input vector will be made
 *@param vec    -> input vector to multiply element wise by scalar, vector will be altered
 *@param scalar -> sclar to multiply vector with
 */
template<typename T>
void vectorTimesScalarNoCopy(std::vector<T>& vec, const T scalar)
{
    std::transform(vec.begin(),
                   vec.end(),
                   vec.begin(),
                   std::bind(std::multiplies<T>(), std::placeholders::_1, scalar));
}

/*******************************************************************************************
 * multiply 2D vector by a scalar without making a copy of the data.
 *
 * @param vec vector which will be rescaled in the routine. Note the original vector will
 * be destroyed.
 * @param scalar by which the 2D vector will be rescaled
 *******************************************************************************************/
template<typename T>
void vectorTimesScalarNoCopy(std::vector<std::vector<T>>& vec, const T scalar)
{
    std::for_each(vec.begin(),
                  vec.end(),
                  [scalar](std::vector<T>& innerVector)
                  {
                      std::transform(innerVector.begin(),
                                     innerVector.end(),
                                     innerVector.begin(),
                                     [scalar](T value) { return value * scalar; });
                  });
}
/** square vector elementwise
 * @param Real vec input vector to square
 * @param output is a Real vector containing the input vector squared
 */
template<typename T>
std::vector<T> squareVector(const std::vector<T>& vec)
{
    std::vector<T> result = vec;

    std::transform(result.begin(), result.end(), result.begin(), [](T x) { return x * x; });
    return result;
}
/** compute dot product between two vectors
 *@param Real vecA input vector has to have same length as vecB
 *@param Real vecB input vector to compute dot product with
 *@param return is the dot product between the two vectors
 */
template<typename T>
T dotProduct(const std::vector<T>& vecA, const std::vector<T>& vecB)
{
    assert((vecA.size() == vecB.size())
           && "dotProduct: dimension of function and xaxis do not match");
    return std::inner_product(vecA.begin(), vecA.end(), vecB.begin(), (T)0);
}

/** do simpson integration on a integer grid
 * @param data array containing the the function to integrate
 */
template<typename T>
T simpsonIntegration(const std::vector<T>& data)
{

    T result = 0;
    T temp1 = 0;
    T temp2 = 0;

    for (auto i = 2; i < data.size(); i += 2)
    {
        temp1 += data[i];
        temp2 += data[i - 1];
    }
    temp1 *= 4.0;
    temp2 *= 2.0;
    result = temp1 + temp2 + data[data.size() - 1] + data[0];
    result /= 3.0;
    return result;
}

/** do simpson integration on a integer grid, 2-dimensional
 * @param data  2d vector containing the function values
 * @param returns the 2dimensional integral on integer grid
 */
template<typename T>
T simpsonIntegration(const std::vector<std::vector<T>>& data)
{
    std::vector<T> integral;
    integral.resize(data.size());
    Real result;

    for (auto i = 0; i < data.size(); i++) { integral[i] = D1_simps(data[i]); }
    result = D1_simps(integral);
    return result;
}

/** compute square of number
 * @param x input number to square
 */
template<typename T>
static inline Real computeSquare(T x)
{
    return x * x;
}

/** compute square of std::vector
 * @param vec_in input vector to square elementqise
 *        on output it contains the squared entries
 */
template<typename T>
void squareVector(std::vector<T>& vec_in)
{
    std::transform(vec_in.begin(), vec_in.end(), vec_in.begin(), computeSquare<T>);
}

/**
 * sum up all elements of the supplied vector
 *@param input vector of numeric type
 */
template<typename T>
T sumVector(const std::vector<T>& input)
{
    return std::accumulate(input.begin(), input.end(), (T)0);
}

/*******************************************************************************************
 * Multiplies two vectors elementwise. Computes Hadarmad product
 *
 * @param vectorA first input vector which has to match the size of vectorB
 * @param vectorB second vectr which is multiplied with vectorA
 *
 * @f[
 * r_{i} = a_{i} * b_{i}
 * @f]
 *******************************************************************************************/
template<typename T>
std::vector<T> elementwiseProduct(const std::vector<T>& vectorA, const std::vector<T>& vectorB)
{
    VASPML_DEBUG_L1(
        if (vectorA.size() != vectorB.size())
        {
            global_scope::tutor.bug("ERROR: math::elementwiseProduct( const std::vector<T>& "
                                    "vectorA, const std::vector<T>& vectorB ) \n"
                                    "length of input vectors does not match");
        }
    );
    std::vector<T> result(vectorA.size());
    std::transform(vectorA.begin(),
                   vectorA.end(),
                   vectorB.begin(),
                   result.begin(),
                   std::multiplies<T>());
    return result;
}

/*******************************************************************************************
 * computes an integer power of some supplied value by explicit for loop
 *
 * @param base value of which the power will be computed
 * @param exp determines the order of the power which is computed of base
 *
 * For small integer exponents an explicit loop to compute the nth power of a function
 * is often faster than the c++ std::pow function. For large integer powers the function
 * intPowSquareAlgo should be prefered.
 *******************************************************************************************/
template<typename T>
inline T intPowLoop(const T& base, const Int& exp)
{

    VASPML_DEBUG_L1(
        if (exp < 0)
        {
            global_scope::tutor.bug("ERROR: inline T intPowLoop( const T& base, const Int& exp )\n"
                                    "exp < 0. Function can only compute powers exp >= 0");
        }
    );

    T result = 1.0;
    for (Int i = 0; i < exp; i++) { result *= base; }
    return result;
}

/*******************************************************************************************
 * computes an integer power of some supplied value by use of explicit cases
 *
 * @param base value of which the power will be computed
 * @param exp determines the order of the power which is computed of base
 *
 * For small integer exponents cases with hard coding the power can be faster than
 * the c++ std::pow function. For large integer powers the function
 * intPowSquareAlgo is called. This function was in tests a little slower than the intPowLoop
 * but still faster than std::pow
 *******************************************************************************************/
template<typename T>
inline T intPowCase(const T& base, const Int& exp)
{
    VASPML_DEBUG_L1(
        if (exp < 0)
        {
            global_scope::tutor.bug("inline T intPowCase( const T& base, const Int& exp ) \n"
                                    "exp < 0. Function can only compute powers exp >= 0");
        }
    );
    switch (exp)
    {
    case 0:
        return 1;
    case 1:
        return base;
    case 2:
        return base * base;
    case 3:
        return base * base * base;
    case 4:
        return base * base * base * base;
    case 5:
        return base * base * base * base * base;
    case 6:
        return base * base * base * base * base * base;
    case 7:
        return base * base * base * base * base * base * base;
    case 8:
        return base * base * base * base * base * base * base * base;
    case 9:
        return base * base * base * base * base * base * base * base * base;
    case 10:
        return base * base * base * base * base * base * base * base * base * base;
    default:
        return intPowLoop(base, exp);
    }
}

/*******************************************************************************************
 * computes an integer power of some supplied value by use of a square algorithm
 *
 * @param baseIn value of which the power will be computed
 * @param expIn determines the order of the power which is computed of base
 *
 * For large integer exponents the square algorithm will be faster than intPowCase, intPowLoop
 * and the std::pow algorithm. But for small integer powers the intPowLoop or intPowCase
 * should be preferred. Testing with gnu compiler showed around exp=16 this function is
 * fastest.
 *******************************************************************************************/
template<typename T>
inline T intPowSquareAlgo(const T& baseIn, const Int& expIn)
{
    VASPML_DEBUG_L1(
        if (expIn < 0)
        {
            global_scope::tutor.bug(
                "inline T intPowSquareAlgo( const T& baseIn, const Int& expIn )\n"
                "exp < 0. Function can only compute powers exp >= 0");
        }
    );
    T   result = 1.0;
    Int exp = expIn;
    T   base = baseIn;
    while (exp)
    {
        if (exp % 2 == 1) { result *= base; }
        base *= base;
        exp /= 2;
    }
    return result;
}

/*******************************************************************************************
 * computes integer power of some supplied value by use of a square algorithm or explicit loop
 *
 * @param base value of which the power will be computed
 * @param exp determines the order of the power which is computed of base
 * For small integer exponents the explicit loop algorithm is used and for larger exponents
 * the intPowSquareAlgo is used.
 *******************************************************************************************/
template<typename T>
inline T intPowMixAlgo(const T& baseIn, const Int& expIn)
{
    VASPML_DEBUG_L1(
        if (expIn < 0)
        {
            global_scope::tutor.bug("inline T intPowMixAlgo( const T& baseIn, const Int& expIn )\n"
                                    "exp < 0. Function can only compute powers exp >= 0");
        }
    );
    if (expIn < 10) return intPowLoop(baseIn, expIn);
    else return intPowSquareAlgo(baseIn, expIn);
}

/*******************************************************************************************
 * compute average value of input vector
 *
 * @param dataIn vector of which elements the average is computed
 *
 * @note if the input vector is empty the returned value is assigned to zero
 *******************************************************************************************/
template<typename T>
T average(const std::vector<T>& dataIn)
{
    if (dataIn.empty()) return (T)0;
    return sumVector(dataIn) / (T)dataIn.size();
}

/*******************************************************************************************
 * compute minimum value of input vector
 *
 * @param dataIn vector of which elements the minimum is computed
 *
 * @note if the input vector is empty the returned value is assigned to zero
 *******************************************************************************************/
template<typename T>
T minimum(const std::vector<T>& dataIn)
{
    if (dataIn.empty()) return (T)0;
    return *std::min_element(dataIn.cbegin(), dataIn.cend());
}

/*******************************************************************************************
 * compute maximum value of input vector
 *
 * @param dataIn vector of which elements the maximum is computed
 *
 * @note if the input vector is empty the returned value is assigned to zero
 *******************************************************************************************/
template<typename T>
T maximum(const std::vector<T>& dataIn)
{
    if (dataIn.empty()) return (T)0;
    return *std::max_element(dataIn.cbegin(), dataIn.cend());
}

/*******************************************************************************************
 * compute prime numbers up to a given number n using Eratosthenes sieve
 *
 * the sieve of Eratosthenes is an ancient algorithm for
 * finding all prime numbers up to any given limit. Routine was taken from vasp and was
 * originally written by Merzuk Kaltak
 *
 * @param[in] n upper bound to which prime numbers are computed
 * @param[out] primSqrt stores square roots of prime numbers which when squared are smaller n
 * @param[out] stores prime numbers up to n
 *******************************************************************************************/
template<typename T>
void computePrimeNumbers(const T n, std::vector<T>& primSqrt, std::vector<T>& primes)
{
    VASPML_DEBUG_L1(
        if constexpr (!std::is_integral<T>())
        {
            global_scope::tutor.bug("computePrimeFactor( const T n, std::vector<T>& "
                                    "primSqrt, std::vector<T>& primes )\n"
                                    "template data type is not an integral type");
        }
    );

    T                 maxPrime = std::floor(std::sqrt((Real)n));
    std::vector<bool> a(n, true);
    for (T i = 2; i <= maxPrime; i++)
    {
        T j = i * i;
        while (true)
        {
            a[j - 1] = false;
            j = j + i;
            if (j > n) break;
        }
    }

    T size = std::count(a.begin(), a.end(), true);
    primes.resize(size - 1);
    T j = 1;
    // skipping 1
    for (T i = 2; i <= n; i++)
    {
        if (!a[i - 1]) continue;
        primes[j - 1] = i;
        j++;
        if (j > (T)primes.size()) break;
    }

    j = 0;
    for (T i = 0; i < (Int)primes.size(); i++)
    {
        if (primes[i] * primes[i] > n)
        {
            j++;
            break;
        }
        j++;
    }
    primSqrt.resize(j);
    for (T i = 0; i < j; i++) primSqrt[i] = primes[i];
}

/*******************************************************************************************
 * compute prime factorization of given input number
 *
 * routine was copied from vasp and originally written by Merzuk Kaltak
 *
 * @param param[in] number for which the prime factorization is obtained
 * @param param[in]  primes prime numbers smaller than n, can be obtained by
 * computePrimeNumbers( const T n, std::vector<T>& primSqrt, std::vector<T>& primes )
 * @param[ out ] nexp prime factors of n
 *******************************************************************************************/
template<typename T>
void primeFactorization(const T n, const std::vector<T>& primes, std::vector<T>& nexp)
{
    VASPML_DEBUG_L1(
        if constexpr (!std::is_integral<T>())
        {
            global_scope::tutor.bug("void primeFactorization( const T n, std::vector<T>& nexp, "
                                    "const std::vector<T>& primes )\n"
                                    "template data type is not an integral type");
        }
    );
    T nprimes = primes.size();
    nexp.resize(primes.size(), (T)0);
    // 1 is not considered as prime
    T m = n;

    for (T ip = 0; ip < nprimes; ip++)
    {
        while (true)
        {
            if (m % primes[ip] == 0)
            {
                m = m / primes[ip];
                nexp[ip]++;
            }
            else { break; }
        }
    }
}
/*******************************************************************************************
 * returns true if supplied value is perfect square ie if it can be expressed as n*n=value
 *
 * routine was copied from vasp and originally written by Merzuk Kaltak
 *
 * @param[ in ] value value to be checked if it is expressible as n * n
 *******************************************************************************************/
template<typename T>
bool isPerfectSquare(const T value)
{
    VASPML_DEBUG_L1(
        if constexpr (!std::is_integral<T>())
        {
            global_scope::tutor.bug("bool isPerfectSquare( const T value )\n"
                                    "template data type is not an integral type");
        }
    );
    bool result;

    T ilimit = std::ceil(std::sqrt((Real)value));

    T i;
    for (i = 1; i <= ilimit; i++)
    {
        if (i * i == value) break;
    }

    if (i > ilimit) result = false;
    else result = true;
    return result;
}
/*******************************************************************************************
 * if the supplied value can be expressed as perfect square this function returns n
 *
 * routine was copied from vasp and originally written by Merzuk Kaltak
 *
 * @param[in] value which will be expressed as the return value n * n
 *******************************************************************************************/
template<typename T>
T perfectSquare(const T value)
{
    VASPML_DEBUG_L1(
      if constexpr ( !std::is_integral<T>() ){
         global_scope::tutor.bug( "T perfectSquare( const T value )\n"
                                  "template data type is not an integral type"
               );
      }
   );

    T ilimit = std::ceil(std::sqrt((Real)value));
    T result;
    for (result = 1; result <= ilimit; result++)
    {
        if (result * result == value) break;
    }
    return result;
}
/*******************************************************************************************
 * check if a supplied number is a prime number; returns true if number is prime
 *
 * routine was copied from vasp and originally written by Merzuk Kaltak
 *
 * @param[in] value value which will be check if it is a prime number
 * @param[in] primes prime numbers wich are smaller equal the supplied value; can be obtained
 * with void computePrimeNumbers( const T n, std::vector<T>& primSqrt, std::vector<T>& primes )
 *******************************************************************************************/
template<typename T>
bool isPrime(const T value, const std::vector<T>& primes)
{
    VASPML_DEBUG_L1(
        if constexpr (!std::is_integral<T>())
        {
            global_scope::tutor.bug("bool isPrime( const T value, const std::vector<T>& primes )\n"
                                    "template data type is not an integral type");
        }
    );
    bool result;
    T    ilimit = (T)primes.size();
    T    i = 0;
    for (T i = 0; i < ilimit; i++)
    {
        if (value == primes[i]) break;
    }
    if (i >= ilimit) result = false;
    else result = true;
    return result;
}

/*******************************************************************************************
 * Factorize an integer number into a product of two numbers.
 *
 * The factorization is done under the constraint that the two numbers are as close as
 * possible to each other. This routine is very helpful when decomposing a parallel
 * environment of MPI into a scalapack grid. The number of processors will be decomposed
 * in a 2 dimensional process grid where the number of rows is as close as possible to
 * the number of columns.
 * routine was copied from vasp and originally written by Merzuk Kaltak
 *
 * @param[in] nprocs is the number which will be decomposed into a product of nprow * npcol
 * @param[out] nprow is the first number
 *
 *******************************************************************************************/
template<typename T>
void fermatRazor(const T nprocs, T& nprow, T& npcol)
{
    VASPML_DEBUG_L1(
        if constexpr (!std::is_integral<T>())
        {
            global_scope::tutor.bug("void fermatRazor( const T nprocs, T nprow, T npcol )\n"
                                    "template data type is not an integral type");
        }
    );
    if (nprocs == 1)
    {
        nprow = 1;
        npcol = 1;
        return;
    }

    std::vector<T> primes;
    std::vector<T> allPrimes;
    // obtain all primes up to sqrt( nprocs  )
    computePrimeNumbers(nprocs, primes, allPrimes);
    std::vector<T> nexp;
    // obtain prime factorization of NPROCS
    primeFactorization(nprocs, primes, nexp);

    // first check if NPROCS is a perfect square
    if (isPerfectSquare(nprocs))
    {
        nprow = perfectSquare(nprocs);
        npcol = nprow;
    }
    else
    {
        // the idea is to factorize nprocs = a*b
        // with restriction |A - B| -> 0
        // fermat sieve works well if N is odd
        T nodd = nprocs / (std::pow(primes[0], nexp[0]));
        // in case nprocs is even
        if (nodd == 1)
        {
            nprow = std::pow(primes[0], nexp[0] / 2);
            npcol = std::pow(primes[0], nexp[0] / 2 + nexp[0] % 2);
        }
        // in case nprocs is not even
        // use fermats sieve on nodd
        else
        {
            // perfect square nothing to do
            if (isPerfectSquare(nodd))
            {
                nprow = perfectSquare(nodd);
                npcol = perfectSquare(nodd);
            }
            else
            {
                if (isPrime(nodd, allPrimes))
                {
                    nprow = 1;
                    npcol = nodd;
                }
                else
                {
                    // Fermat's way to factorize an odd integer:
                    // nodd = a**2- b**2 = ( a - b ) * ( a + b )
                    // and write a**2 - n = b**2
                    // and increase A by 1 until A**2 - N is a perfect square
                    T a = std::ceil(std::sqrt((Real)nodd));
                    T b = a * a - nodd;
                    while (true)
                    {
                        if (isPerfectSquare(b)) break;
                        a++;
                        b = a * a - nodd;
                    }
                    nprow = a + perfectSquare(b);
                    npcol = a - perfectSquare(b);
                } // nodd is prime
                for (T i = 0; i < nexp[0]; i++)
                {
                    if (nprow < npcol) nprow = nprow * 2;
                    else npcol = npcol * 2;
                }
            } // perfect square nodd
        } // nprocs not even
    } // perfec square

    if (nprow * npcol != nprocs)
    {
        global_scope::tutor.error("ERROR: void fermatRazor( const T nprocs, T& nprow, T& npcol )\n"
                                  "failed "
                                  + std::to_string(nprocs) + "  " + std::to_string(nprow)
                                  + std::to_string(npcol) + "\n");
    }
}

} // namespace vaspml::math
#endif
