#ifndef DESCRIPTORCOLLECTOR_HPP
#define DESCRIPTORCOLLECTOR_HPP

#include "ArrayResizing.hpp"
#include "Descriptor.hpp"
#include "Linalg.hpp"
#include "ParallelEnvironemt.hpp"
#include "SmartEnum.hpp"
#include "TypeMap.hpp"
#include "types.hpp"

#include <map>
#include <memory>

// debug
#include "utils.hpp"

namespace vaspml
{

enum class DescriptorStorage
{
    Type,
    CentralAtom,
};

class DescriptorCollector
{

  public:
    /**
     *@param storageOrder determines the storage order of the arrays.
     *Setting 0 makes type order and setting it 1
     *makes central atom ordering
     *@param normalizedDescriptors here a map array can be supplied
     *if the memory should not be managed by the class
     *itself
     @param normsAtom supply a storage space where the descriptor norms are stored per central atom
     */
    DescriptorCollector(
        Int storageOrder = 0,
        const std::map<String, ShVec2Real>& normalizedDescriptors = std::map<String, ShVec2Real>(),
        const ShVec1Real& normsAtom = nullptr,
        const ExecutionPolicy algoExecution = ExecutionPolicy::cpuSingleCore);
    /**
     * updating the descriptor collector class
     *
     *@param descriptors map which contains the descriptors that should be collected and normalized.
     *The objects in the map have to be inherited from the Descriptor class
     */
    void updateCollector(void);
    /**
     *get the normalized descriptors as array by const reference
     */
    const std::map<String, ShVec2Real>& get_descriptorsNormalized(void) const;
    /*******************************************************************************************
     * get norm over descriptor list for single atom
     *
     * @param atomIndx index of central atom
     *******************************************************************************************/
    Real get_normAtom(std::size_t atomIndx) const;
    /*******************************************************************************************
     * get normalization factors of central atoms
     *
     * normalization factors are computed over all descriptors of the specific atom
     *
     * @note call with const Vec1Real& x = get_normAtom(); to avoid making copy
     *******************************************************************************************/
    const Vec1Real& get_normAtom(void) const;
    /*******************************************************************************************
     * rearrange normalized descriptors in case that ML_FF file contains more types than
     * the actual structure.
     *
     * this routine selects CPU or GPU version
     *******************************************************************************************/
    void rearrangeSHS2Body(const TypeMap& typeMap);
    /*******************************************************************************************
     * rearrange normalized descriptors - cpu version using std::lib
     *
     * no parallel execution policy is used in routine
     *******************************************************************************************/
    void rearrangeSHS2BodyCPU(const TypeMap& typeMap);

    /*******************************************************************************************
     * rearrange normalized descriptors - gpu version using std::lib
     *******************************************************************************************/
    void rearrangeSHS2BodyGPU(const TypeMap& typeMap);
    /*******************************************************************************************
     * returns one of the collected Descriptor in non normalized form.
     *
     * Can for example be used to retrieve neighbor lists
     *******************************************************************************************/
    const Descriptor& getDescriptor(const std::string& key) const;

    void setDescriptorMap(const std::map<String, std::shared_ptr<Descriptor>>& descriptors);

    void writeDescriptorCollector(void) const;

  private:
    /**
     *Normalizing the supplied descriptors and storing them in the chosen format.
     *Choses whether CPU or GPU version is used.
     *
     *@param descriptors map storing the descriptors
     *@NOTE the descriptors are not normalized independently, but are normalized
     *as supervectors per central atom
     */
    void normalizeDescriptors(void);
    /**
     *Normalizing the supplied descriptors and storing them in the chosen format
     *
     *CPU implementation using for loops and standard library algorithms.
     *
     *@param descriptors map storing the descriptors
     *@NOTE the descriptors are not normalized independently, but are normalized
     *as supervectors per central atom
     */
    void normalizeDescriptorsCPU(void);
    /**
     *normalizing the supplied descriptors and storing them in the chosen format
     *GPU implementation using standard library algorithms.
     *
     *@param descriptors map storing the descriptors
     *@NOTE the descriptors are not normalized independently, but are normalized
     *as supervectors per central atom
     */
    void normalizeDescriptorsGPU(void);
    /**
     *allocate arrays descriptorsNormalized in type ordered or central atom format
     *@param descriptors map storing the descriptors to retrieve the needed sizes for the arrays
     */
    void allocateArrays(void);
    void allocateArraysCPU(void);
    void allocateArraysGPU(void);
    /**
     * compute norm of a single atom over the set of all descriptors
     *
     * @param atomIndex index of central atom for which the descriptor norm is computed
     */
    Real computeNormSingleAtom(std::size_t atomIndex);
    /**
     *rescale descriptors by the computed norms and multiply times the
     *descriptor weight store in type order format
     *
     *@param atomIndex index of central atom
     *@param norm computed norm of descriptor for certain central atom and summed over descriptor
     *types
     *@param descriptors map of classes to obtain the non-normalized descriptor elements
     */
    void normalizeDescriptorsSingleAtomTypeOrder(std::size_t atomIndex, const Real norm);
    /**
     *rescale descriptors by the computed norms and multiply times the
     *descriptor weight; store in Central atom format
     *
     *@param atomIndex index of central atom
     *@param norm computed norm of descriptor for certain central atom and summed over descriptor
     *types
     *@param descriptors map of classes to obtain the non-normalized descriptor elements
     */
    void normalizeDescriptorsSingleAtomCentralOrder(std::size_t atomIndex, const Real norm);
    /*******************************************************************************************
     * Resize descriptorsNormalized and fill with 0.
     *
     * This function will only be called if the number of force field types does not match
     * the number of types in the structure.
     *******************************************************************************************/
    void allocateDescriptorSHS2BodyNormalized(const TypeMap& typeMap);

    /**
     * vector of the descriptors needed by the descriptor collector class
     */
    std::map<String, std::shared_ptr<Descriptor>> descriptors;
    /**
     * std::map to store the normalized for the supplied descriptors.
     * @NOTE
     * the array can either be stored as [type][ atoms of type x number descriptors ] or
     * [ central atom ][ number descriptors ]
     */
    std::map<String, ShVec2Real>      descriptorsNormalized;
    std::map<String, ArrayResizing2D> descriptorsNormalizedSize;
    /**
     *smart enum to determine the storage order
     *
     *supported order types are type ordered. Type index determines first index
     *of storage arrays; Second index is a combination index of atoms x number of descriptors.
     *Central atom index ordering is second possibility. Then the first index is the central atom
     *index and the second index of the storage arrays is the number of descriptors
     */
    DescriptorStorage storageOrder;
    /**
     * storing the length of the array descriptorsNormalized
     *@NOTE
     * can be stored in a type or central atom fashion
     */
    std::map<String, ShVec1Int>       length_descriptorsNormalized;
    std::map<String, ArrayResizing1D> length_descriptorsNormalizedSize;
    /**
     * norm computed over descriptor list
     *
     *@Note order is always over central atom. This order will match those of the descriptors
     *so when doing
     *loop type:
     *    loop atom per type:
     *the order matches normsAtom when incrementing by one
     */
    ShVec1Real      normsAtom;
    Vec1Size_t      centralAtomIndex;
    ArrayResizing1D centralAtomIndexSize;
    // std::ofstream fileSHS2;
    // std::ofstream fileSHS3;
    /***********************************************************************
     * Number of atoms in the system
     ************************************************************************/
    Int                   numberAtoms;
    ExecutionPolicy       algoExecution;
    Vec1String            descriptorKeyList;
    linalg::LinalgContext linalgContext;
};

} // namespace vaspml

#endif
