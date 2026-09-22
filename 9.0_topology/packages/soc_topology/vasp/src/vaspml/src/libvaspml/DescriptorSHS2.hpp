#ifndef DESCRIPTORSHS2_HPP
#define DESCRIPTORSHS2_HPP

#include "ArrayResizing.hpp"
#include "BasisFunctions.hpp"
#include "Descriptor.hpp"
#include "Linalg.hpp"
#include "ParallelEnvironemt.hpp"
#include "ThermodynamicIntegration.hpp"
#include "TypeMap.hpp"
#include "nearest_neighbor.hpp"
#include "types.hpp"

#include <cstddef>
#include <memory>

namespace vaspml
{

class DescriptorSHS2 : public Descriptor
{
  public:
    /*******************************************************************************************
     * @class DescriptorSHS3
     * @brief Two-body descriptor and it's derivative relevant stuff.
     *
     * @param weight descriptor weight when computing normalization constant
     * @param isNormalized determines if descriptor will partake in normalization
     * @param clnmPair supply if memory for clnmPair is supplied from outside
     * @param clnmPairDerivative supply if memory for clnmPair is supplied from outside
     * @param clnmVasp supply if memory for clnmPair is supplied from outside
     * @param clnm_clnmDerivativeCentralVasp supply if memory for clnmPair is supplied from outside
     *
     * This class is used to create a two-body (SHS2) descriptor.
     * This class can for example be created as part of a descriptor collector
     * DescriptorMap which is of std::map type.
     * @code
     * DescriptorMap descriptors;
     * @endcode
     * It can then be accessed with the key "SHS2-2-body" for the two-body descriptors
     * or by "SHS2-3-body" for the two-body part of the thre-body descriptors.
     * The filling for the "SHS2-2-body" descriptors would look like the following (in the Frame
     *class):
     * @code
     * descriptors["SHS2-2-body"] = std::make_shared<DescriptorSHS2>(
     * inputParameters["SHS2-2-body-weight"].cget<Real>(),
     * inputParameters["SHS2-2-body-is-normalized"].cget<bool>(),
     * std::get<ShVec2Real>(arrayData.at("2-body-SHS2-pair")),
     * std::get<ShVec2Real>(arrayData.at("2-body-SHS2-pair-derivative")),
     * std::get<ShVec2Real>(arrayData.at("2-body-SHS-vasp")),
     * std::get<ShVec2Real>(arrayData.at("2-body-derivative-SHS2-central-vasp")));
     * @endcode
     *
     *******************************************************************************************/
    DescriptorSHS2(const Real            weight = (Real)1.0,
                   const bool            isNormalized = false,
                   const ShVec2Real&     clnmPair = nullptr,
                   const ShVec2Real&     clnmPairDerivative = nullptr,
                   const ShVec2Real&     clnmVasp = nullptr,
                   const ShVec2Real&     clnmDerivativeCentralVasp = nullptr,
                   const DescriptorType  descriptorType = DescriptorType::none,
                   const ExecutionPolicy algoExecution = ExecutionPolicy::cpuSingleCore);

    /****************************************************************************************
     * Updating the spherical harmonic pair coefficients
     * with the current neighbor list and a supplied basis set.
     *
     * function is only wrapper to decide what kind of algorithm execution will be choosen.
     * For implementation check functions updatePairCoefficientsSingleCPU,
     *  updatePairCoefficientsGPU
     *
     * @param nn_list nearest neighbor list
     * @param basisFunctions basis functions for descriptor
     ****************************************************************************************/
    void updatePairCoefficients(const std::shared_ptr<NearestNeighborNSquare>& nn_list,
                                const std::shared_ptr<BasisFunctions>&         basisFunctions);
    void updateVaspCoefficients(const std::shared_ptr<NearestNeighborNSquare>& nn_list);
    /****************************************************************************************
     * Compute clnmVasp, this is a sum over the clnmPair over the nearest neighbors,
     * by keeping type resolution.
     * This routine selects which algorithm to use (CPU or openACC GPU).
     *******************************************************************************************/
    void computeVaspCoefficientsFromPairCoefficients(void);

    /****************************************************************************************
     * Compute clnmVasp, this is a sum over the clnmPair over the nearest neighbors,
     * by keeping type resolution.
     * CPU implementation using standard library algorithms.
     *******************************************************************************************/
    void computeVaspCoefficientsFromPairCoefficientsCPU(void);

    /****************************************************************************************
     * Compute clnmVasp, this is a sum over the clnmPair over the nearest neighbors,
     * by keeping type resolution.
     * GPU implementation using standard library algorithms.
     *******************************************************************************************/
    void computeVaspCoefficientsFromPairCoefficientsGPU(void);
    /****************************************************************************************
     * Compute thermodynamic integration coefficients.
     * @param thermodynVars contains storage arrays for thermodynamic integration
     *******************************************************************************************/
    void computeThermodynamicIntegration(ThermodynIntegration& thermodynVars);

    /****************************************************************************************
     * Get the clnmPair vector for a certain centralAtom and all its neighbors.
     * @param centralAtom central atom index for which cnlm pairs are obtained
     *
     * @note The order of the clnm pair is neighbor,l,n,m.
     *******************************************************************************************/
    const Vec1Real& get_clnmPair(const std::size_t centralAtom) const;

    /****************************************************************************************
     * Obtain the clnmPairDerivative vector for a certain centralAtom and all its neighbors.
     * @param centralAtom central atom index for which cnlm pairs are obtained
     *
     * @note The order of the clnm pair is neighbor,xyz,l,n,m.
     *******************************************************************************************/
    const Vec1Real& get_clnmPairDerivative(const std::size_t centralAtom) const;

    /****************************************************************************************
     * Get the intermediate expansion coefficients clnmVasp for a certain centralAtom.
     *
     * This is a version of cnlm_pair which was summed over neighbors with same types.
     *
     * @param centralAtom atom index for which the coefficients are returned
     *
     * @note format is atom type, l,n,m.
     *******************************************************************************************/
    const Vec1Real& get_clnmVasp(const std::size_t centralAtom) const;

    /****************************************************************************************
     * Same as get_clnmVasp, but implements the virtual function of
     * Descriptor used in the DescriptorCollector class.
     *
     * @param centralAtom central atom index for which the coefficients are returned
     *******************************************************************************************/
    const Vec1Real& get_descriptor(const std::size_t centralAtom) const;
    const Vec2Real& get_descriptor(void) const;

    /****************************************************************************************
     * Get size of spherical harmonics expansion coefficient for a certain central atom.
     *
     * @param centralAtom central atom index for which size is returned
     *******************************************************************************************/
    Int get_sizeDescriptor(const std::size_t centralAtom) const;

    /****************************************************************************************
     * Rescale spherical harmonic expansion coeffcients for a certain central atom.
     *
     * @param centralAtom central atom index for which the coefficients are returned
     * @param scaleFactor scaling factor by which the coefficients are rescaled
     *******************************************************************************************/
    void rescale_descriptor(const std::size_t centralAtom, const Real scaleFactor);

    /****************************************************************************************
     * Get the clnmDerivativeCentralVasp vector for a certain centralAtom.
     *******************************************************************************************/
    const Vec1Real& get_clnmDerivativeCentralVasp(const std::size_t centralAtom);

    /****************************************************************************************
     * Return clnmVasp for a certain central atom and a nearest neighbor.
     *
     * @param central_atom central index as nearest neighbor list index
     * @param neighbor     nearest neighbor index as in the NearestNeighborNSquare
     * @param l            angular quantum number of spherical haromnics
     * @param n            determines the nth root of the modified spherical Bessel function
     * @param m            magnetic quantum number of spherical harmonics
     *******************************************************************************************/
    const Real& get_clnmVasp(const std::size_t central_atom,
                             const std::size_t type_index,
                             const std::size_t l,
                             const std::size_t n,
                             const std::size_t m) const;

    /****************************************************************************************
     *
     * Get maximal order of spherical harmonics expansion coefficients.
     *
     * @return maxOrder-1
     *******************************************************************************************/
    std::size_t get_maxOrder(void) const;

    /****************************************************************************************
     *
     * Get number of roots/radial basis functions for given l.
     *
     * @param l order of spherical harmonic
     *
     * @return number of roots for order l
     *
     *******************************************************************************************/
    const std::size_t& get_nRootsOrder(const std::size_t l) const;

    /****************************************************************************************
     *
     * Get number of roots of radial part of basis functions.
     *
     * @return const referecce to number roots of basis function
     *
     *******************************************************************************************/
    const std::vector<std::size_t>& get_nRoots(void) const;

    /****************************************************************************************
     *
     *Return the index for a specifig @[fc^{i}_{lnm}@f].
     *
     *@atom atom index
     *@l angular index (index of angular basis function (spherical harmonic))
     *@n radial index (index of radial basis function)
     *@m second index to define the spherical harmonic
     *******************************************************************************************/
    Size_t get_Index(const std::size_t atom,
                     const std::size_t l,
                     const std::size_t n,
                     const std::size_t m) const;
    /****************************************************************************************
     *
     *Return the basis set size of the two-body descriptor.
     *
     *******************************************************************************************/
    std::size_t get_basisSetSize(void) const;

    /****************************************************************************************
     *
     * Get offset within clnm for a given l quantum number.
     *
     * Function gives back variable lOffset[l].
     *
     *@l angular quantum number
     *******************************************************************************************/
    std::size_t get_lOffset(const std::size_t l) const;

    /*******************************************************************************************
     * Computing the derivative of the kernel matrix times the derivatives of the descriptors of
     *this class. Overlaying function, will choose between compute_forcePreContractCPU and
     * compute_forcePreContractGPUOpenACC.
     *
     * This function can be used in combination with the classes Kernel and Predictor
     * to compute pairwise forces acting on atoms.
     * @param derivativeMatrix derivative of the Kernel matrix for specific atom type
     * @param forcePreContract is the derivative of the kernel matrix times the descriptor
     *derivative will be computed in this function and has to be zero on input
     * @param typeMap used to convert  force field types to structure types
     *******************************************************************************************/
    void compute_forcePreContract(const Vec2Real& derivativeMatrix,
                                  Vec2Real&       forcePreContract,
                                  const TypeMap&  typeMap) const;

    /*******************************************************************************************
     * Computing the derivative of the kernel matrix times the derivatives of the descriptors of
     *this class. CPU version which uses standard algorithms.
     *
     * This function can be used in combination with the classes Kernel and Predictor
     * to compute pairwise forces acting on atoms.
     * @param derivativeMatrix derivative of the Kernel matrix for specific atom type
     * @param forcePreContract is the derivative of the kernel matrix times the descriptor
     *derivative will be computed in this function and has to be zero on input
     * @param typeMap used to convert  force field types to structure types
     *******************************************************************************************/
    void compute_forcePreContractCPU(const Vec2Real& derivativeMatrix,
                                     Vec2Real&       forcePreContract,
                                     const TypeMap&  typeMap) const;

    /*******************************************************************************************
     * Computing the derivative of the kernel matrix times the derivatives of the descriptors of
     *this class. GPU version which uses standard algorithms.
     *
     * This function can be used in combination with the classes Kernel and Predictor
     * to compute pairwise forces acting on atoms.
     * @param derivativeMatrix derivative of the Kernel matrix for specific atom type
     * @param forcePreContract is the derivative of the kernel matrix times the descriptor
     *derivative will be computed in this function and has to be zero on input
     * @param typeMap used to convert  force field types to structure types
     *******************************************************************************************/
    void compute_forcePreContractGPU(const Vec2Real& derivativeMatrix,
                                     Vec2Real&       forcePreContract,
                                     const TypeMap&  typeMap) const;
    /*******************************************************************************************
     * Compute the size of the array where the pre contraction of the derivativeMatrix and the
     * descriptor is stored. Routine is needed in the class Predictor.
     *******************************************************************************************/
    std::size_t get_forcePreContractSize(const std::size_t atomIndex) const;
    /*******************************************************************************************
     * Returns the lOffsets used in index computation routines.
     * Takes into account, that only upper triangular part for nRadial1 and nRadial2 is stored.
     *
     * call with std::vector<std::size_t>& variableName = get_lOffset();
     *******************************************************************************************/
    const std::vector<std::size_t>& get_lOffset(void) const;
    /*******************************************************************************************
     * Compute pair forces and central forces derived from the kernel approach.
     * Routine decides whether CPU or GPUOpenACC version is called.
     *
     * @param forcePreContract stores the pre contraction of the derivativeMatrix times the 2-body
     *descriptor
     * @param pairForces stores the position derivative of the kernel with respect to
     * the neighbors of the considered central atom
     * @param centralForces stores the position derivative of the kernel with respect to the
     *considered central atom
     *******************************************************************************************/
    void computeForceTerms(const Vec2Real& forcePreContract,
                           Vec2Real&       pairForces,
                           Vec1Real&       centralForces) const;
    /*******************************************************************************************
     * Compute pair forces and central forces derived from the kernel approach.
     * CPU implementation using standard library algorithms.
     *
     * @param forcePreContract stores the pre contraction of the derivativeMatrix times the 2-body
     *descriptor
     * @param pairForces stores the position derivative of the kernel with respect to
     * the neighbors of the considered central atom
     * @param centralForces stores the position derivative of the kernel with respect to the
     *considered central atom
     *******************************************************************************************/
    void computeForceTermsCPU(const Vec2Real& forcePreContract,
                              Vec2Real&       pairForces,
                              Vec1Real&       centralForces) const;
    /*******************************************************************************************
     * Compute pair forces and central forces derived from the kernel approach.
     * GPU implementation using standard library algorithms.
     *
     * @param forcePreContract stores the pre contraction of the derivativeMatrix times the 2-body
     *descriptor
     * @param pairForces stores the position derivative of the kernel with respect to
     * the neighbors of the considered central atom
     * @param centralForces stores the position derivative of the kernel with respect to the
     *considered central atom
     *******************************************************************************************/
    void computeForceTermsGPU(const Vec2Real& forcePreContract,
                              Vec2Real&       pairForces,
                              Vec1Real&       centralForces) const;
    /*******************************************************************************************
     * setting the basis functions which are used for computation of local atom descriptors.
     *
     * @param basisFunctions sets the class basis functions to thois input
     *******************************************************************************************/
    void set_basisFunctions(const std::shared_ptr<BasisFunctions>& basisFunctions);
    void write_clnmVasp(const String& fname) const;

  private:
    /****************************************************************************************
     * Updating the spherical harmonic pair coefficients
     * with the current neighbor list and a supplied basis set.
     *
     * Algorithm will run on CPU. Function is called from updatePairCoefficients
     *
     * @param basisFunctions basis functions for descriptor
     ****************************************************************************************/
    void updatePairCoefficientsCPU(const std::shared_ptr<BasisFunctions>& basisFunctions);
    void updateVaspCoefficientsCPU(void);
    /****************************************************************************************
     * Updating the spherical harmonic pair coefficients
     * with the current neighbor list and a supplied basis set.
     *
     * Algorithm will run on GPU. Function is called from updatePairCoefficients
     *
     * @param nn_list nearest neighbor list
     * @param basisFunctions basis functions for descriptor
     ****************************************************************************************/
    void updatePairCoefficientsGPU(const std::shared_ptr<BasisFunctions>& basisFunctions);
    void updateVaspCoefficientsGPU(void);
    /*******************************************************************************************
     * Compute Offsets for indexing of the storage arrays.
     *
     * generates lOffset shifts between increments of angular quantum number l and stores
     * them to lOffset
     *******************************************************************************************/
    void computeOffsets(void);
    /*******************************************************************************************
     * Compute pair coefficients clnmPair and the derivatives clnmPairDerivative
     * for a certain central atom.
     *
     * @param index central atom index under consideration
     * @param basisFunctions class containing basis functions to get expansion coefficients
     * clnmPair and the derivatives clnmPairDerivative
     *******************************************************************************************/
    void computePairCoefficientSingleAtom(const std::size_t                      index,
                                          const std::shared_ptr<BasisFunctions>& basisFunctions,
                                          const NearestNeighborNSquare&          neighborList);
    void computeVaspCoefficientSingleAtom(const std::size_t index);
    /*******************************************************************************************
     * Compute basisSetSize, the total number of basis functions denoted by l,n,m.
     *
     * @param lmax maximum angular quantum number of spherical harmonic basis function
     *******************************************************************************************/
    std::size_t compute_BasisSetSize(const std::size_t               lmax,
                                     const std::vector<std::size_t>& nRoots);
    /****************************************************************************************
     * Allocating the arrays cnlm_pair and cnlm_pair_derivative according
     * to the nearest neighbor list topology.
     *
     * @param[ in ] ldim number of spherical harmonics
     *******************************************************************************************/
    void allocatePairCoefficientArrays(const Size_t& ldim);
    /****************************************************************************************
     * Allocate arrays clnmVasp and clnmDerivativeCentralVasp.
     *******************************************************************************************/
    void allocateArraysVaspFormat(void);
    /****************************************************************************************
     * Compute the index pairs to assign a specific descriptor derivative term to an entry in the
     * derivative of the kernel.
     *******************************************************************************************/
    std::size_t computeIndexDerivativeMatrix(std::size_t centralAtom,
                                             std::size_t typeIndexForceField,
                                             std::size_t radialIndex,
                                             std::size_t totalShift);
    void        resizePairArraysGPU(const Size_t& ldim);
    void        resizeArraysVaspFormatGPU(void);
    void        resize_centralAtomIndex(void);
    void        resize_radialBasisFunction(void);
    void        resize_ylm(const Size_t& ldim);
    /****************************************************************************************
     * clnmPair summed over nearest neighbors.
     *
     * type resolution is kept
     * first index is central atom
     * second index is combined index
     * type 1 all lnm
     * type 2 all lnm
     * .
     * .
     * .
     *******************************************************************************************/
    ShVec2Real      clnmVasp;
    ArrayResizing2D clnmVaspSize;
    /****************************************************************************************
     * clnmPairDerivative summed over nearest neighbors.
     *
     * type resolution is kept
     * first index is central atom
     * second index is combined index
     * x component type 1 all lnm
     * y component type 1 all lnm
     * z component type 1 all lnm
     * x component type 2 all lnm
     * y component type 2 all lnm
     * z component type 2 all lnm
     * .
     * .
     * .
     *******************************************************************************************/
    ShVec2Real      clnmDerivativeCentralVasp;
    ArrayResizing2D clnmDerivativeCentralVaspSize;
    /****************************************************************************************
     * Store expansion coefficients of the nearest neighbor environment.
     *
     *  index 1 is the central atom, index 2 is
     *  neighbor 1 all lnm
     *  neighbor 2 all lnm
     *  neighbor 3 all lnm...
     *******************************************************************************************/
    ShVec2Real      clnmPair;
    ArrayResizing2D clnmPairSize;
    /****************************************************************************************
     * Derivatives of the expansion coefficients of the nearest neighbor environment.
     *
     * index 1 central atom
     * index 2 neighbor atom
     * index 3 x all nlm, y all nlm, z all nlm
     *******************************************************************************************/
    ShVec2Real      clnmPairDerivative;
    ArrayResizing2D clnmPairDerivativeSize;
    /****************************************************************************************
     * Total number of basis functions (n,l,m).
     *******************************************************************************************/
    std::size_t basisSetSize;
    /****************************************************************************************
     * Number of roots of radial function.
     *
     * @todo Maybe this should become a shared pointer at some point, also in the descriptorRadial
     * spline routine, like this a copy constructor will be called.
     *******************************************************************************************/
    std::vector<std::size_t> nRoots;
    /****************************************************************************************
     * Total number of radial basis functions.
     *******************************************************************************************/
    std::size_t nRadialBasis;
    /****************************************************************************************
     * Number of atoms in the nearest neighbor list that was supplied to construct the object.
     *******************************************************************************************/
    std::size_t nAtoms;
    /****************************************************************************************
     * Max Order of spherical haromincs, stores max angular momentum + 1 for loop writing
     *convenience.
     *
     *@note When returned with getter the the plus 1 will be removed.
     *******************************************************************************************/
    std::size_t maxOrder;
    /****************************************************************************************
     * Array offset for indexing.
     *******************************************************************************************/
    std::vector<std::size_t> lOffset;
    /****************************************************************************************
     * To check if pair coefficients are computed.
     *******************************************************************************************/
    bool pairsComputed;
    /****************************************************************************************
     * To check if calculations in vasp format were carried out.
     *******************************************************************************************/
    bool vaspFormatComputed;
    /*******************************************************************************************
     * determines the execution policy for the computaionally demanding algorithms of class
     *******************************************************************************************/
    Vec1Real        radialBasisFunction;
    Vec1Real        radialBasisFunctionDerivative;
    ArrayResizing1D radialBasisFunctionSize;

    Vec2Real                        ylm;
    ArrayResizing2D                 ylmSize;
    Vec2Real                        ylmDerivative;
    ArrayResizing2D                 ylmDerivativeSize;
    Vec1Size_t                      centralAtomIndex;
    ArrayResizing1D                 centralAtomIndexSize;
    std::shared_ptr<BasisFunctions> basisFunctions;
    bool                            basisFunctionsSet;
    linalg::LinalgContext           linalgContext;
};

} // namespace vaspml

#endif
