#ifndef DESCRIPTOR_HPP
#define DESCRIPTOR_HPP

#include "ParallelEnvironemt.hpp"
#include "TypeMap.hpp"
#include "nearest_neighbor.hpp"
#include "types.hpp"

namespace vaspml
{

enum class DescriptorType
{
    none = 0,
    bodyOrder2 = 1,
    bodyOrder3 = 2,
    bodyOrder3LinearElement = 3
};

class DescriptorSHS2;
class DescriptorSHS3;
class DescriptorSHS3ReducedLinElem;
class Descriptor
{

  public:
    Descriptor(void);
    /*******************************************************************************************
     * @param weight weiggt factor of descriptor used during normalization
     * @param isNormalzied bool which defines if the descriptor is normalized in
     * DescriptorCollector
     *******************************************************************************************/
    Descriptor(const Real            weight,
               const bool            isNormalized,
               const DescriptorType  descriptorType,
               const ExecutionPolicy algoExecution = ExecutionPolicy::cpuSingleCore);
    /**
     *return the weight factor which determines the influence of the descriptor in a list of
     *descriptors
     */
    Real get_weight(void) const;
    /**
     * get controll parameter of the descriptor will be normalized in the DescriptorCollector for
     * example
     */
    bool get_isNormalized(void) const;
    /**
     *set the neighbor list from which the descriptor is computed
     */
    void set_neighborList(const std::shared_ptr<const NearestNeighborNSquare>& neighborList);
    /**
     * return neighbor list on wich the pair coefficients are based
     *
     * @return const reference to neighbor list
     */
    const NearestNeighborNSquare& get_neighborList(void) const;

    /**
     * return neighbor list on wich the pair coefficients are based
     *
     * @return syd::shared_ptr to neighbor list
     */
    std::shared_ptr<const NearestNeighborNSquare> get_neighborList_ptr(void) const;
    /** Get type index of central atom of some clnm
     *
     * @param atomIndex array entry of desired central atom
     */
    const Int&     get_typeIndexCentral(const std::size_t atomIndex) const;
    const Vec1Int& get_typeIndexCentral(void) const;
    /** Get type index of indxNeigh neighbor from central atom atomIndex
     *
     * @param atomIndex central atom
     * @param indxNeigh number of the neighbor in NearestNeighborNSquare
     */
    const Int& get_typeIndex(const std::size_t atomIndex, const std::size_t indxNeigh) const;
    /**
     * give number of atoms for which the expansion coefficients were computed
     *
     * @return number of atoms
     */
    Int get_nAtoms(void) const;
    /**
     *get number of atoms per type from nearest neighbor list
     */
    const Vec1Int& get_nAtomsType(void) const;
    /**
     * get descriptor array for certain central atom
     *@param centralAtom central atom index for which descriptor is retrieved
     */
    const Vec1Real& get_descriptor(const Size_t centralAtom) const;
    /**
     * get the size of the descriptor, get number of descriptor for central atom
     *@param centralAtom atom index for which the number of descriptors are returned
     *@NOTE
     *the number of descriptors is independent of the atom index. So all atoms
     *from for example the SHS2 descriptor for a given neighbor list will have
     *the same number of descriptors
     */
    Int get_sizeDescriptor(const std::size_t centralAtom) const;
    /**
     * rescale the descriptors for a certain central atom by a constant factor
     *
     *@param centralAtom atom index for which descriptors are rescaled
     *@param scaleFactor factor by which the elements of the descriptor are scaled
     */
    void rescale_descriptor(const std::size_t centralAtom, const Real scaleFactor);
    /*******************************************************************************************
     * Computing the product of the kernel derivative times the derivative of the descriptor
     * vector with respect to the SHS2-2-body descriptors.
     *
     * @param derivativeMatrix stores the derivative of the kernel matrix with respect
     * to descriptor vector
     * @param forcePreContract stores the result given by the product of the kernel
     * matrix times the descriptor
     * @param typeMap used to convert the structure types to the machine learning force field types
     *
     * This equation shows the what product is computed by the routine.
     * @f[\sum_{B}w_{B}\frac{\partial K( \hat{\mathbf{X}}_{i},\hat{\mathbf{X}}_{B} ) }
     * {\partial \hat{\mathbf{X}}_{i} } \frac{d\ \hat{\mathbf{X}}_{i}}{c^{iJ_{1}}_{lnm}}@f]
     * For a more information on computing this product check the child classes
     * DescriptorSHS2, DescriptorSHS3
     *******************************************************************************************/
    void compute_forcePreContract(const Vec2Real& derivativeMatrix,
                                  Vec2Real&       forcePreContract,
                                  const TypeMap&  typeMap) const;
    /*******************************************************************************************
     * compute size for the array where forcePreContract is saved to.
     *******************************************************************************************/
    Size_t get_forcePreContractSize(const std::size_t atomIndex) const;
    /*******************************************************************************************
     * computing central and pair force terms. pair force terms are derivatives of kernel matrix
     * with respect to neighbor atom positions and central forces are those where derivative is
     *taken w.r.t to central atom
     *
     * @param forcePreContract stores the derivative of the kernel matrix times the 2-body
     *descriptors. This array has to be multiplied by 2-body descriptor derivatives with respect to
     *positions to get force terms.
     * @param pairForces on output stores the derivative of the kernel matrix with respect to the
     *positions of the neighbor atoms for a certain central atom. which is the force pair between
     *two atoms
     * @param centralForces on output stores the derivative of the kernel matrix wrt to the
     *considered central atom
     *******************************************************************************************/
    void computeForceTerms(const Vec2Real& forcePreContract,
                           Vec2Real&       pairForces,
                           Vec1Real&       centralForces) const;

  private:
    /// weight of the descriptor when combining several to supervector in DescriptorCollector class
    Real weight;
    /// define if the descriptor will be normalized in the Descriptor collector
    bool isNormalized;
    /// shared ptr to neighbor list from which the descriptor was computed
    std::shared_ptr<const NearestNeighborNSquare> neighborList;
    /// check if derivative sparse map is already computed. Has to be done once only
    bool                  derivativeMapComputed;
    const DescriptorType  descriptorType;
    const ExecutionPolicy algoExecution;
    friend class DescriptorSHS2;
    friend class DescriptorSHS3;
    friend class DescriptorSHS3ReducedLinElem;
};
} // namespace vaspml
#endif
