#ifndef CUBIC_SPLINE_HPP
#define CUBIC_SPLINE_HPP

#include "ParallelEnvironemt.hpp"

#include <iostream>
#include <type_traits>
#include <vector>

namespace vaspml
{

template<class T>
class CubicSpline
{

    /** this class computes a cubic spline interpolation
     * for a input function -> function
     * given on a grid -> axis
     * The routine is adapted from the Cubic spline
     * interpolation given in
     * Numerical Recipes in C++: The Art of Scientific Computing 3rd Edition
     * by  William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery
     * page 120
     *  ISBN-10  0521880688
     *  ISBN-13  978-0521880688
     */

    // excluding non numeric types. Somehow c++
    // considers char as a numeric type
    static_assert(std::is_arithmetic<T>::value, "Wrong argument type in CubicSpline\n");

  public:
    void computey2(void);

    /*
     * compute first derivative of cubic spline function on supplied grid
     *@param x vector on which the derivatives of the spline are computed
     *@return derivtiaves of interpolated spline
     */
    std::vector<T> interpolateVectorDerivative(const std::vector<T>& x) const;
    /*
     * compute derivative of interpolated spline fuction  for supplied
     * @param x value on which to interpolate
     */
    T interpolateValueDerivative(const T x) const;
    /** compute an interpolated function value at position x
     * @param x   value at which the function value should be computed
     */
    T interpolateValue(const T x) const;
    /** compute an interpolated function values for given x array
     * @param x  array containingn x-values on which the function has to
     *                be computed
     */
    std::vector<T> interpolateVector(const std::vector<T>& x) const;

    /*
     * compute spline interpolated function and it's derivative
     *@param x value on which to interpolate
     *@param func    interpolated function value at x, return value
     *@param funcDer derivative of interpolated function at x, return value
     */
    void computeSplineAndDerivative(const T x, T& func, T& funcDer) const;

    /*
     * compute spline interpolated function and it's derivative
     *
     *@param x vector on which to interpolate
     *@param func    interpolated function value at x-positions, return array
     *@param funcDer derivative of interpolated function at x-positions, return array
     */
    VASPML_NV_HOST_DEVICE
    void computeSplineAndDerivative(const std::vector<T>& x,
                                    std::vector<T>&       func,
                                    std::vector<T>&       funcDer) const;

    /// default constructor, initializing empty object
    CubicSpline(void) {};
    /** constructor which computes the first and second derivatives
     *  during construction
     * @param xaxis     support points on which function is given
     * @param function  function values at support points
     * @param yp1       first derivative of function at first support point
     * @param ypn       first derivative of function at last support point
     */
    CubicSpline(const std::vector<T>& xaxis,
                const std::vector<T>& function,
                const T               yp1 = 1.e99,
                const T               ypn = 1.e99);
    /// give second derivatives to the user
    std::vector<T> get_y2(void) const;
    /// give first derivatives to the user
    std::vector<T> get_u(void) const;

    /** set all needed parameters for the CubicSpline object at once
     * @param xaxis     support points on which function is given
     * @param function  function values at support points
     * @param yp1       first derivative of function at first support point
     * @param ypn       first derivative of function at last support point
     */
    void setInput(const std::vector<T>& xaxis_in,
                  const std::vector<T>& function_in,
                  const T               yp1_in = 1.e99,
                  const T               ypn_in = 1.e99);
    /// set xaxis on which interpolation quantities are
    /// computed
    void setXaxis(const std::vector<T>& axis_in);
    /// set function to interpolate
    void setFunction(const std::vector<T>& function_in);
    /// set first derivative on first point
    void setYp1(const T yp1_in);
    /// set first derivative on last point
    void setYpn(const T ypn_in);
    /// set all input values and compute the second derivative
    void update(const std::vector<T>& xaxis_in,
                const std::vector<T>& function_in,
                const T               yp1_in = 1.e99,
                const T               ypn_in = 1.e99);

  private:
    /// xaxis on which the function to interpolate is given
    std::vector<T> xaxis;
    /// function values corresponding to xaxis
    std::vector<T> function;
    /// second derivatives of function
    std::vector<T> y2;
    /// first derivative of the function
    std::vector<T> u;
    /// first derivative of the function
    std::vector<T> v;
    /// grid spacing
    T grid_spacing;
    /// derivative at first grid point
    T yp1;
    /// derivative at first grid point
    T ypn;
};

template<class T>
CubicSpline<T>::CubicSpline(const std::vector<T>& axis_in,
                            const std::vector<T>& function_in,
                            const T               yp1_in,
                            const T               ypn_in)
{
    xaxis = axis_in;
    function = function_in;
    yp1 = yp1_in;
    ypn = ypn_in;
    computey2();
}

template<class T>
std::vector<T> CubicSpline<T>::get_y2(void) const
{
    return y2;
}

template<class T>
std::vector<T> CubicSpline<T>::get_u(void) const
{
    return u;
}

template<class T>
void CubicSpline<T>::computey2(void)
{

    T            p;
    T            qn;
    T            sig;
    T            un;
    unsigned int n = xaxis.size();
    // set array sizes
    y2.resize(n);
    u.resize(n);
    // store gridspacing
    grid_spacing = xaxis[1] - xaxis[0];

    if (yp1 > 0.99e99)
    {
        y2[0] = (T)0.0;
        u[0] = (T)0.0;
    }
    else
    {
        y2[0] = (T)-0.5;
        u[0] = ((T)3.0 / (xaxis[1] - xaxis[0]))
             * ((function[1] - function[0]) / (xaxis[1] - xaxis[0]) - yp1);
    }

    for (size_t i = 1; i < n - 1; i++)
    {
        sig = (xaxis[i] - xaxis[i - 1]) / (xaxis[i + 1] - xaxis[i - 1]);
        p = sig * y2[i - 1] + (T)2;
        y2[i] = (sig - (T)1) / p;
        u[i] = (function[i + 1] - function[i]) / (xaxis[i + 1] - xaxis[i])
             - (function[i] - function[i - 1]) / (xaxis[i] - xaxis[i - 1]);
        u[i] = ((T)6.0 * u[i] / (xaxis[i + 1] - xaxis[i - 1]) - sig * u[i - 1]) / p;
    }

    if (ypn > 0.99e99)
    {
        qn = (T)0.0;
        un = (T)0.0;
    }
    else
    {
        qn = (T)0.5;
        un = ((T)3 / (xaxis[n - 1] - xaxis[n - 2]))
           * (ypn - (function[n - 1] - function[n - 2]) / (xaxis[n - 1] - xaxis[n - 2]));
    }

    y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + (T)1.0);
    for (int i = n - 2; i >= 0; i--) { y2[i] = y2[i] * y2[i + 1] + u[i]; }
}

template<class T>
T CubicSpline<T>::interpolateValue(const T x) const
{

    size_t klow = (size_t)(x / grid_spacing);
    size_t khigh = klow + 1;
    T      h = xaxis[khigh] - xaxis[klow];
    if (h == 0) { throw("Bad value in CubicSpline::interpolateValue"); }
    T a = (xaxis[khigh] - x) / h;
    T b = (x - xaxis[klow]) / h;
    return a * function[klow] + b * function[khigh]
         + ((a * a * a - a) * y2[klow] + (b * b * b - b) * y2[khigh]) * h * h / (T)6;
}

template<class T>
T CubicSpline<T>::interpolateValueDerivative(const T x) const
{

    size_t klow = (size_t)(x / grid_spacing);
    size_t khigh = klow + 1;
    T      h = xaxis[khigh] - xaxis[klow];
    if (h == 0) { throw("Bad value in CubicSpline::interpolateValue"); }
    T a = (xaxis[khigh] - x) / h;
    T b = (x - xaxis[klow]) / h;
    return (((T)3 * a * a - (T)1) * y2[klow] + ((T)3 * b * b - (T)1) * y2[khigh]) * h / (T)6;
}

template<class T>
std::vector<T> CubicSpline<T>::interpolateVector(const std::vector<T>& x) const
{
    std::vector<T> y(x.size());
    for (auto i = 0; i < x.size(); i++) { y[i] = interpolateValue(x[i]); }
    return y;
}

template<class T>
std::vector<T> CubicSpline<T>::interpolateVectorDerivative(const std::vector<T>& x) const
{
    std::vector<T> y(x.size());
    for (auto i = 0; i < x.size(); i++) { y[i] = interpolateValueDerivative(x[i]); }
    return y;
}
template<class T>
void CubicSpline<T>::computeSplineAndDerivative(const T x, T& f, T& fder) const
{

    size_t klow = (size_t)(x / grid_spacing);
    size_t khigh = klow + 1;

    T h = xaxis[khigh] - xaxis[klow];
#ifndef USE_NVIDIA_GPU
    if (h == 0)
    {
        std::cout << "Bad value in CubicSpline::interpolateValue" << std::endl;
        std::cout << "Upper bound of grid reached" << std::endl;
        throw;
    }
#endif

    T a = (xaxis[khigh] - x) / h;
    T b = (x - xaxis[klow]) / h;

    f = a * function[klow] + b * function[khigh]
      + ((a * a * a - a) * y2[klow] + (b * b * b - b) * y2[khigh]) * h * h / (T)6;

    fder = (function[khigh] - function[klow]) / h
         - (((T)3 * a * a - (T)1) * y2[klow] - ((T)3 * b * b - (T)1) * y2[khigh]) * h / (T)6;
}

template<class T>
void CubicSpline<T>::computeSplineAndDerivative(const std::vector<T>& x,
                                                std::vector<T>&       func,
                                                std::vector<T>&       funcDer) const
{
    func.resize(x.size());
    funcDer.resize(x.size());
    for (auto i = 0; i < x.size(); i++) { computeSplineAndDerivative(x[i], func[i], funcDer[i]); }
}

template<class T>
void CubicSpline<T>::update(const std::vector<T>& xaxis_in,
                            const std::vector<T>& function_in,
                            const T               yp1_in,
                            const T               ypn_in)
{
    setXaxis(xaxis_in);
    setFunction(function_in);
    setYp1(yp1_in);
    setYpn(ypn_in);
    computey2();
}

template<class T>
void CubicSpline<T>::setInput(const std::vector<T>& xaxis_in,
                              const std::vector<T>& function_in,
                              const T               yp1_in,
                              const T               ypn_in)
{
    setXaxis(xaxis_in);
    setFunction(function_in);
    setYp1(yp1_in);
    setYpn(ypn_in);
}

template<class T>
void CubicSpline<T>::setXaxis(const std::vector<T>& xaxis_in)
{
    xaxis = xaxis_in;
}

template<class T>
void CubicSpline<T>::setFunction(const std::vector<T>& function_in)
{
    function = function_in;
}

template<class T>
void CubicSpline<T>::setYp1(const T yp1_in)
{
    yp1 = yp1_in;
}

template<class T>
void CubicSpline<T>::setYpn(const T ypn_in)
{
    ypn = ypn_in;
}

} // namespace vaspml

#endif
