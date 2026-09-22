#ifndef BASISFUNCTIONSANGULAR_HPP
#define BASISFUNCTIONSANGULAR_HPP

#include "BasisFunctionsRadialSpline.hpp"

#include <cstddef>

namespace vaspml
{

namespace math
{
class SphericalHarmonics;
}

class BasisFunctionsAngular : public BasisFunctionsRadialSpline
{

  public:
    BasisFunctionsAngular();
    BasisFunctionsAngular(math::CutoffType  cutOffType_in,
                          Real              cutOff_in,
                          Real              widthBroadening_in,
                          std::size_t       nGrid_in,
                          std::size_t       maxOrderBessel_in,
                          std::size_t       nRootsBessel_in,
                          BasisFunctionType type,
                          std::size_t       angularFiltering_in = 9999999,
                          Real              filterScale_in = 0);

    //**********************************************************************************************
    // Regular member functions.
    //**********************************************************************************************
    void computeAngularBasis(const Vec1Real& rhat,
                             const Vec1Real& rnorm,
                             Vec1Real&       ylm,
                             Vec1Real&       ylmd) const;

    //**********************************************************************************************
    // Getters and setters.
    //**********************************************************************************************
    // Setters
    void set_lmax(std::size_t lmax_in);
    // Getters
    const Vec1Real& get_filteringFactors() const;
    std::size_t     get_ldim() const;

  private:
    /// Maximum angular quantum number.
    std::size_t lmax;
    /// Number of @f$l,m@f$ combinations for given @f$l_{\max}$, equals size of flattened
    /// SphericalHarmonics arrays.
    std::size_t                               ldim;
    std::shared_ptr<math::SphericalHarmonics> angularBasis;
    /** Angular filtering type.
     *
     * Angular filtering decreases importance of angular basis with increasing l.
     */
    std::size_t angularFiltering;
    /// Scaling parameter for filtering function.
    Real filterScale;
    /// Computes filteringFactors.
    void computeAngularFilteringCoefficient();
    /// Filtering factors for each angular quantum number.
    Vec1Real filteringFactors;
};

} // namespace vaspml

#endif
