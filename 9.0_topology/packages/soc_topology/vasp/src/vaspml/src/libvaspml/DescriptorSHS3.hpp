#ifndef DESCRIPTORSHS3_HPP
#define DESCRIPTORSHS3_HPP

#include "ArrayResizing.hpp"
#include "Descriptor.hpp"
#include "DescriptorSHS2.hpp"
#include "ParallelEnvironemt.hpp"
#include "TypeMap.hpp"
#include "nearest_neighbor.hpp"
#include "types.hpp"

#include <array>
#include <cstddef>
#include <map>
#include <memory>
#include <vector>

namespace vaspml
{

class DescriptorSHS3 : public Descriptor
{
  public:
    /***************************************************************************************
     * @class DescriptorSHS3
     * @brief Three-body descriptor and it's derivative relevant stuff.
     *
     * @descriptorList list containing the used descriptors. This list is obtained from the
     * vasl ML_FF file
     * @angularFilterOn switching on and off angular filtering
     * @angularFilterType determine type of angular filtering function
     * @angularFilterScaling scaling parameter of angular filtering function
     * @param nRoots number of roots of modified spherical Bessel functions for various angular
     *parameter
     * @maxOrder maximal order of spherical Bessel functions/ spherical harmonics
     * @param spectrum_vasp optional shared_ptr where the SHS3 in vasp format is stored
     * @param spectrum_pair pair wise SHS3 spectrum where atom resolution is included
     *
     * This class is used to create a three-body descriptor. This descriptor
     * is derived from the two-body descriptor.
     * This class can for example be created with as part of a descriptor collector
     * DescriptorMap which is a std::map.
     * @code
     * DescriptorMap descriptors;
     * @endcode
     * It can then be accessed with the key "SHS3-3-body", but compared to the
     * regular three-body descriptor the shared pointer needs to be called
     * with the keyword "LinearScalingDescriptor".
     * The filling would look like the following (in the Frame class):
     * @code
     * descriptors["SHS3-3-body"] = std::make_shared<DescriptorSHS3>(
     *    inputParameters["SHS3-3-body-descriptor-list"].dcget<ShVec2Int>(),
     *    inputParameters["SHS3-3-body-angular-filter-on"].cget<bool>(),
     *    inputParameters["SHS3-3-body-angular-filter-type"].cget<Int>(),
     *    inputParameters["SHS3-3-body-angular-filter-scale"].cget<Real>(),
     *    basis3Body.get_nValidRoots(),
     *    inputParameters["SHS3-3-body-max-angular-number"].cget<Int>(),
     *    inputParameters["SHS3-3-body-weight"].cget<Real>(),
     *    inputParameters["SHS3-3-body-is-normalized"].cget<bool>(),
     *    std::get<ShVec2Real>(arrayData.at("3-body-SHS3-spectrum")),
     *    std::get<ShVec2Real>(arrayData.at("3-body-SHS3-pair-spectrum")));
     * @endcode
     *
     ***************************************************************************************/
    DescriptorSHS3(const Vec2Int&        descriptorList,
                   const bool            angularFilterOn,
                   const Int             angularFilterType,
                   const Real            angularFilterScaling,
                   const Vec1Size_t&     nRoots,
                   const Int             maxOrder,
                   const Real            weight = 1.0,
                   const bool            isNormalized = false,
                   const ShVec2Real&     spectrum_vasp = nullptr,
                   const ShVec2Real&     spectrum_pair = nullptr,
                   const DescriptorType  descriptorType = DescriptorType::none,
                   const ExecutionPolicy algoExecution = ExecutionPolicy::cpuSingleCore);

    /****************************************************************************************
     * Compute SHS3 (3-body spectrum ) for a given neighbor list configuration.
     *
     * @param shec spherical harmonics expansion coefficients for a given neighbor list
     * @typeMap map class which associates machine-learning force field types with types in
     * current structure
     *
     ****************************************************************************************/
    void computeSHS3(const std::shared_ptr<DescriptorSHS2>& shec, const TypeMap& typeMap);

    /*******************************************************************************************
     * Getter for 3-body descriptor (spectrum) in vasp format for a given central atom.
     *
     * @centralAtom central atom index which mactches neighbor list and expansion
     * coefficients shec
     *******************************************************************************************/
    const Vec1Real& get_SHS3Atom_vasp(const Int centralAtom) const;

    /*******************************************************************************************
     * Function returns same thing as get_SHS3Atom_vasp, but implements the virtual function of
     * Descriptor used in the DescriptorCollector class.
     *
     * @param centralAtom central atom index for which the coefficients are returned
     *******************************************************************************************/
    const Vec1Real& get_descriptor(const std::size_t centralAtom) const;

    /*******************************************************************************************
     * Returns the size of 3-body descriptor.
     *
     * @param centralAtom central atom index for which the coefficients are returned
     *******************************************************************************************/
    Int get_sizeDescriptor(const std::size_t centralAtom) const;

    /*******************************************************************************************
     * Rescale 3-body descriptor (spectrum).
     *
     * @param centralAtom central atom index for which descriptor is rescaled
     * @param centralAtom central atom index for which descriptor is rescaled
     *******************************************************************************************/
    void rescale_descriptor(const std::size_t centralAtom, const Real scaleFactor);

    /*******************************************************************************************
     * Compute the derivative of a kernel matrix times the derivatives of the descriptors of this
     *class.
     *
     * @param derivativeMatrix kernel matrix derivative
     * @param forcePreContract on output contains derivativeMatrix times descriptor derivative
     * @param typeMap translates atom types from current structure to types stored in force field
     *
     * @f[ -\sum\limits_{i}\mathbf{L}_{i} \frac{d\mathbf{X}_{i}}{dp_{nn'l}^{iJJ'}}
     *      \frac{dp_{nn'l}^{iJJ'}}{dc_{n''lm}^{iJ''}}@f]
     *******************************************************************************************/
    void compute_forcePreContract(const Vec2Real& derivativeMatrix,
                                  Vec2Real&       forcePreContract,
                                  const TypeMap&  typeMap) const;

    /*******************************************************************************************
     * Return the size of the second dimension (number of Cnlm) of the array forcePreContract.
     *
     * @param atomIndex index of central atom for which memory allocation is done
     *
     * This routine is used in the Predictor class for allocation (resize) of the
     * forcePreContract array. The size of forcePreContract is
     * [central atom ][ neighborType x l x nRadial x m ].
     *******************************************************************************************/
    std::size_t get_forcePreContractSize(const std::size_t atomIndex) const;

    /*******************************************************************************************
     * Compute centralForces and pairForces which are needed to compute total force on ion.
     *
     * @param forcePreContract array which stores the derivative of the
     *        kernelMatrix times the descriptor derivative
     * @param pairForces on output storing pair forces in same format as neighbor list
     * @param centralForces force acting on central atom by making derivative wrt to atom itself
     *
     * Routine is called from Predictor class to be able to predict total force
     * in predictor class. Routine implements a wrapper to the DescriptorSHS2 class
     * which implements the functionality.
     * @note For more implementation details check
     * compute_forcePreContractSparse. Currently
     * only compute_forcePreContractSparse is used.
     *******************************************************************************************/
    void computeForceTerms(const Vec2Real& forcePreContract,
                           Vec2Real&       pairForces,
                           Vec1Real&       centralForces) const;

  private:
    /***************************************************************************************
     * Loop over atoms for computeSHS3SingleAtom.
     *
     * Single CPU version, using a for loop.
     *
     * @typeMap map class which associates machine-learning force field types with types in
     * current structure
     ****************************************************************************************/
    void computeSHS3CPU(const TypeMap& typeMap);
    /***************************************************************************************
     * Loop over atoms for computeSHS3SingleAtom.
     *
     * GPU version, using standard library algorithms.
     *
     * @shec shared pointer to two-body descriptor
     * @typeMap map class which associates machine-learning force field types with types in
     * current structure
     ****************************************************************************************/
    void computeSHS3GPU(const TypeMap& typeMap);
    /*******************************************************************************************
     * Compute the product of the derivative matrix times the 2-body descriptor times a scalar.
     *
     * Single CPU version, uses c++ standard algorithms for computations.
     *
     * @param derivativeMatrix stores derivative of kernel
     *matrix.@f[L^{iJ_{1}J_{2}}_{ln_{1}n_{2}}@f]
     * @param forcePreContract stores the result. The product of the derivativeMatrix and 2 body
     *descriptor
     * @param typeMap gives the ability to convert between structure and force field types
     *
     * @f[
     *f^{iJ_{1}}_{ln_{1}m}=\sum_{J_{2}}\sum_{n_{2}}a_{l}*L^{iJ_{1}J_{2}}_{ln_{1}n_{2}}c^{iJ_{2}}_{l
     *n_{2}m} @f]
     * @f[
     *\tilde{f}^{iJ_{2}}_{ln_{2}m}\sum_{J_{1}}\sum_{n_{1}}a_{l}*L^{iJ_{1}J_{2}}_{ln_{1}n_{2}}c^{iJ_{1}}_{l
     *n_{1}m} @f] these terms are then summed up.
     * @f[ F^{iJ}_{lnm}=f^{iJ}_{lnm}+f^{iJ}_{lnm}@f]
     *
     * @note The algorithm does not compute the two terms from the equations first and sums them
     *afterwards. The algorithm uses an index list representation to compute both terms in one loop.
     *******************************************************************************************/
    void compute_forcePreContractSparseSingleCPU(const Vec2Real& derivativeMatrix,
                                                 Vec2Real&       forcePreContract,
                                                 const TypeMap&  typeMap) const;
    /*******************************************************************************************
     * Compute the product of the derivative matrix times the 2-body descriptor times a scalar.
     *
     * GPU version, uses c++ standard algorithms for computations.
     *
     * @param derivativeMatrix stores derivative of kernel
     *matrix.@f[L^{iJ_{1}J_{2}}_{ln_{1}n_{2}}@f]
     * @param forcePreContract stores the result. The product of the derivativeMatrix and 2 body
     *descriptor
     * @param typeMap gives the ability to convert between structure and force field types
     *
     * @f[
     *f^{iJ_{1}}_{ln_{1}m}=\sum_{J_{2}}\sum_{n_{2}}a_{l}*L^{iJ_{1}J_{2}}_{ln_{1}n_{2}}c^{iJ_{2}}_{l
     *n_{2}m} @f]
     * @f[
     *\tilde{f}^{iJ_{2}}_{ln_{2}m}\sum_{J_{1}}\sum_{n_{1}}a_{l}*L^{iJ_{1}J_{2}}_{ln_{1}n_{2}}c^{iJ_{1}}_{l
     *n_{1}m} @f] these terms are then summed up.
     * @f[ F^{iJ}_{lnm}=f^{iJ}_{lnm}+f^{iJ}_{lnm}@f]
     *
     * @note The algorithm does not compute the two terms from the equations first and sums them
     *afterwards. The algorithm uses an index list representation to compute both terms in one loop.
     *******************************************************************************************/
    void compute_forcePreContractSparseGPU(const Vec2Real& derivativeMatrix,
                                           Vec2Real&       forcePreContract,
                                           const TypeMap&  typeMap) const;

    //      void compute_forcePreContractSparseGPUOpenACC( const Vec2Real& derivativeMatrix,
    //                                           Vec2Real& forcePreContract,
    //                                           const TypeMap& typeMap ) const;

    /*******************************************************************************************
     * Compute index lists for compute_forcePreContract array. This is a necessary prerequisite
     * for calling compute_forcePreContractSparse.
     *
     * computes: \n
     * sparseMap_derivativeMatrix \n
     * sparseMap_SHS2 \n
     * sparseMap_central \n
     * preFactorPreContract \n
     * To implement the product of the matrix @f$L^{iJ_{1}J_{2}}_{ln_{1}n_{2}}@f$ and
     * @f$\frac{d\ p^{iJ_{1}J_{2}}_{ln_{1}n_{2}}}{d\ c^{iJ_{3}}_{ln_{3}m}}@f$ the
     * following selection rules have to be obeyed:
     * @f[ \frac{d\ p^{iJ_{1}J_{2}}_{ln_{1}n_{2}}}{d\ c^{iJ_{3}}_{ln_{3}m}}=
     *    a_{l}\left[ \delta_{n_{1}n_{3}} \delta_{J_{1}J_{3}}c^{iJ_{2}}_{ln_    {2}m_{1}} +
     *    \delta_{n_{2}n_{3}} \delta_{J_{2}J_{3}}c^{iJ_{1}}_{ln_{1}m} \right] @f]
     *
     * The lists are implemented as sparse lists for
     * @f$\frac{d\ p^{iJ_{1}J_{2}}_{ln_{1}n_{2}}}{d\ c^{iJ_{3}}_{ln_{3}m}}@f$.
     *******************************************************************************************/
    void prepareSparseforcePreContract(const TypeMap& typeMap);

    /*******************************************************************************************
     * Make list of the active descriptors used in the current force field.
     *
     * @descriptorList list containing the used descriptors. This list is obtained from the
     * vasp ML_FF file
     * @nRoots the zeros of the modified spherical Bessel functions as computed
     * in desriptorRadialSpline
     * @maxOrder maximal order of the angular parameter of the spherical harmonics
     *******************************************************************************************/
    void make_sparseList(const Vec2Int&    descriptorList,
                         const Vec1Size_t& nRoots,
                         const Int         maxOrder);

    /*******************************************************************************************
     * Precompute anfgular filtering functions.
     *
     * @param maxOrder maximal order of angular parameter of spherical harmonics
     *
     * @note more information can be found https://www.vasp.at/wiki/index.php/ML_IAFILT2
     *******************************************************************************************/
    void make_angularFilter(const Int maxOrder);

    /*******************************************************************************************
     * Allocate storage arrays for spectrumVasp.
     *
     * @typeMap map class which associates machine-learning force field types with types in
     * current structure
     *
     * @shec spherical harmonics expansion coefficients precomputed for current neighbor list
     *instance
     *******************************************************************************************/
    void allocateArrays(const TypeMap& typeMap);
    void resizeArraysGPU(const TypeMap& typeMap);
    /*******************************************************************************************
     * Compute the 3-body spectrum (stored in spectrumVasp) for given central atom.
     *
     * @shec  spherical harmonics expansion coefficients from which the 3-body spectrum is computed
     * @map class which associates machine-learning force field types with types in
     * current structure
     * @atom central atom index ordered as in neighhbor list or in shec
     *******************************************************************************************/
    void computeSHS3SingleAtom(const DescriptorSHS2& shec, const std::size_t atom);
    /*****************************************************
     * Number of elements which are present in force field.
     *****************************************************/
    std::size_t numberElementsMLFF;
    /****************************************************
     * Number of elements which are present in structure.
     ****************************************************/
    std::size_t numberElementsStructure;
    /*******************************************
     * Contains the SHS3 (3-body descriptor).
     *
     * First index is the central atom index.
     * Second index is a combined index of
     * type1, type2, l, n0, n1.
     ******************************************/
    ShVec2Real      spectrumVasp;
    ArrayResizing2D spectrumVaspSize;
    /*******************************************
     * Pair wise storage of SHS3 (3-body descriptor).
     *
     * Not used at the moment.
     ******************************************/
    ShVec2Real spectrumPair;

    /****************************************************
     * Number of atoms for which power spectrum is stored.
     ****************************************************/
    Int numberAtoms;

    /********************************************
     * List of centram atom types.
     *
     * @Note The type lists here need a converter
     * from force field type to structure type.
     *********************************************/
    std::vector<std::vector<std::size_t>> type0List;

    /***********************************************
     * List of neighbor types 1.
     *
     * This list is sorted such that one can loop
     * over the whole power spectrum in a linear way.
     ***********************************************/
    std::vector<std::vector<std::size_t>> type1List;

    /************************************************
     * List of neighbor types 2.
     *
     * This list is sorted such that one can loop
     * over the whole power spectrum in a linear way.
     ************************************************/
    std::vector<std::vector<std::size_t>> type2List;

    /************************************************
     * List of angular number of spherical harmonic l.
     *
     * This list is sorted such that one can loop
     * over the whole power spectrum in a linear way.
     ************************************************/
    std::vector<std::vector<std::size_t>> angularList;

    /******************************************************
     * Tist of first radial number ( number of roots ) of spherical bessel function.
     *
     * This list is sorted such that one can loop
     * over the whole power spectrum in a linear way.
     ******************************************************/
    std::vector<std::vector<std::size_t>> n0List;

    std::vector<std::size_t> n0ListSize;

    /******************************************************************
     * Tist of second radial number ( number of roots ) of spherical bessel function.
     *
     * This list is sorted such that one can loop
     * over the whole power spectrum in a linear way.
     ******************************************************************/
    std::vector<std::vector<std::size_t>> n1List;

    /*********************************************************
     * Weight factors to compute n0List n1List terms.
     *
     * The weight factor will be 1 for n0List == n1List
     * and will be sqrt(2) otherwise. This alllows to compute
     * only the upper triangular matrix of n0 n1 terms.
     *
     *********************************************************/
    Vec2Real weightFactor;

    /****************************************************
     * Controls whether the angular filtering is switched on or not.
     *
     * For more information check the vasp wiki
     * https://www.vasp.at/wiki/index.php/ML_IAFILT2
     *
     ****************************************************/
    bool angularFilterOn;

    /****************************************************
     * Determine which type of angular filtering is used.
     ****************************************************/
    Int angularFilterType;

    /****************************************************
     * Drray stores the precomputed angular filter.
     *
     * Dimension is maxOrder + 1, filter value for every l is stored.
     ****************************************************/
    Vec1Real angularFilter;

    /****************************************************************
     * Scaling factor to determine the width
     * of the distribution with which the angular filtering is done.
     ****************************************************************/
    Real angularFilterScaling;

    /********************************************************************
     * Store SHS2 descriptor which was used to build the SHS3 descriptor.
     ********************************************************************/
    std::shared_ptr<const DescriptorSHS2> descriptorSHS2;

    /************************************************************
     * Check if the sparse maps for sparsePreContract were set up.
     ************************************************************/
    bool isSparseforcePreContractReady;

    /*******************************************************************************************
    * Sparse map which shows the entries in the derivativeMatrix which
    * are needed for the computation of the array forcePreContract.
    *
    * The derivativeMatrix is given by
    * @f[ L^{iJ_{1}J_{2}}_{ln_{1}n_{2}} =
          \frac{\partial K(\hat{\mathbf{X}}_{i},\hat{\mathbf{X}}_{B})}
          {\partial \hat{\mathbf{X}}_{i}} @f]
    *@f[ i @f] central atom index
    *@f[J_{1} @f] index denoting first neighbor type
    *@f[J_{2} @f] index denoting second neighbor type
    *@f[ l @f] indexing angular parameter
    *@f[ n_{1} @f] radial index belonging to neighbor type 1
    *@f[ n_{2} @f] radial index belonging to neighbor type 2
    * The order in which the indices are described denotes storage order.
    * The sparse list stores those elements of the derivativeMatrix for which the equation is
    non-zero.
    *
    * The actual derivative matrix is not stored in this class, but for example in Predictor.hpp.
    *******************************************************************************************/
    std::vector<std::vector<std::size_t>> sparseMap_derivativeMatrix;

    /**************************************************************
     * Sparse map which shows the entries in the SHS2-descriptor which
     * are needed for the computation of the array forcePreContract.
     *
     * forcePreContract is computed by the following equation. Einstein sum convention is used
     * over @f[J_{2}@f] and @f[ n_{2}@f] in first term. And in the second term the sums are taken
     over
     * @f[J_{1}@f] and @f[n_{2}@f]. The central atom index i is summed over in both terms
     * @f[
       a_{l}*L^{iJ_{1}J_{2}}_{ln_{1}n_{2}}c^{iJ_{2}}_{l
     n_{2}m}+a_{l}*L^{iJ_{1}J_{2}}_{ln_{1}n_{2}}c^{iJ_{1}}_{l n_{1}m}
       @f]
     *@f[ L^{iJ_{1}J_{2}}_{ln_{1}n_{2}} @f] derivativeMatrix
     *@f[ c^{iJ_{2}}_{l n_{2}m} @f] descriptorSHS2
     *@f[ a_{l} @f] preFactorPreContract
     **************************************************************/
    std::vector<std::vector<std::size_t>> sparseMap_SHS2;

    /*******************************************************************************************
    * Sparse map which shows the entries in the forcePreContract array (central atom index)
    * which are needed for the computation of the array forcePreContract
    *
    * forcePreContract is computed by the following equation. Einstein sum convention is used
    * over @f[J_{2}@f] and @f[ n_{2}@f] in first term. And in the second term the sums are taken
    over
    * @f[J_{1}@f] and @f[n_{2}@f]
    * @f[
      a_{l}*L^{iJ_{1}J_{2}}_{ln_{1}n_{2}}c^{iJ_{2}}_{l
    n_{2}m}+a_{l}*L^{iJ_{1}J_{2}}_{ln_{1}n_{2}}c^{iJ_{1}}_{l n_{1}m}
      @f]
    * The result will therefore be a matrix @f[ Y^{kJ}_{lnm}@f] this map stores to which entry of
    * Y the result is written.
    *******************************************************************************************/
    std::vector<std::vector<std::size_t>> sparseMap_central;

    /*******************************************************************************************
    * Pre-factors for the forcePreContract array computation
    * which are needed for the computation of the array forcePreContract.
    *
    * The prefactors can be expressed with Kronecker-Deltas:
    * @f[
      f^{iJ_{1}J_{2}}_{ln_{1}n_{2}} =
      \begin{cases}
      2 * a_{l} &\text{if } J_{1}=J_{2} \land n_{1}=n_{2} \\
      a_{l} &\text{if } J_{1}\neq J_{2} \land n_{1}=n_{2}\\
      \sqrt{2}a_{l} &\text{else}
      \end{cases}
    @f]
     the factor @f[a_{l}@f] is stored in angularFilter.
    *******************************************************************************************/
    std::vector<std::vector<Real>> preFactorPreContract;

    Vec1Size_t      centralAtomIndex;
    ArrayResizing1D centralAtomIndexSize;
    Vec1Size_t      descriptorSize;
    TypeMap         typeMapLoc;
};

} // namespace vaspml
 
#endif
