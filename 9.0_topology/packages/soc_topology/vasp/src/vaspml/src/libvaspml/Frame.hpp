#ifndef FRAME_HPP
#define FRAME_HPP

#include "BasisFunctions.hpp"
#include "Descriptor.hpp"
#include "IoHandlerML_FF.hpp"
#include "Lattice.hpp"
#include "ParallelEnvironemt.hpp"
#include "SmartEnum.hpp"
#include "Structure.hpp"
#include "TypeMap.hpp"
#include "nearest_neighbor.hpp"
#include "types.hpp"

#include <map>
#include <memory>

namespace vaspml
{

using NeighborListMap = std::map<String, std::shared_ptr<NearestNeighborNSquare>>;
using BasisFunctionMap = std::map<String, std::shared_ptr<BasisFunctions>>;
using DescriptorMap = std::map<String, std::shared_ptr<Descriptor>>;

enum class Descriptor3BodyType
{
    StandardDescriptor,
    LinearScalingDescriptor
};

/***************************************************************************************************
 * This class contains one structure, the neighbor lists and its descriptors.
 **************************************************************************************************/
class Frame
{
  public:
    /***********************************************************************************************
     * Set up empty frame class instance.
     *
     * @warning Before use, the init() member function needs to be called.
     **********************************************************************************************/
    Frame() = default;
    /***********************************************************************************************
     * Initialize an empty frame class instance from given force field parameters.
     *
     * This method uses createFrameMemory(), createNeighborLists() and createDescriptors() to set up
     * a Frame instance ready to be supplied with an actual structure. It prepares the multi-type
     * map for all required dynamic memory, then creates NeighborList and Descriptor instances.
     *
     * @param inputParameters Instance of force field setup.
     * @param basis3Body Instance of BasisFunctionsAngular (?) class.
     **********************************************************************************************/
    void init(const IoHandlerML_FF& inputParameters,
              const BasisFunctions& basis3Body,
              const ExecutionPolicy algoExecution = ExecutionPolicy::cpuSingleCore);
    /***********************************************************************************************
     * Compute all descriptors, neighbor list is assumed to be present and correctly filled.
     *
     * @param basis3Body Instance of BasisFunctionsAngular (?) class.
     **********************************************************************************************/
    void update(const BasisFunctionMap& basisFunctions);
    void update(void);
    void set_basisFunctions(BasisFunctionMap& basisFunctions);
    /***********************************************************************************************
     * Read structure from file in POSCAR format and compute neighbor lists and descriptors.
     *
     * @param fileName Name of file containing structure in POSCAR format.
     * @param basis3Body Instance of BasisFunctionsAngular (?) class.
     **********************************************************************************************/
    void update(const String& fileName, const BasisFunctionMap& basisFunctions);
    /***********************************************************************************************
     * update descriptors from previosuly filled neighbor list array
     *
     * @param basisFunctions is a dictionary storing basis functions for various
     * types of descriptors as the SHS2-2-body and the SHS3-3-body.
     * @note SHS2-2-body and SHS2-3-body
     * use BasisFunctionsRadialSpline. SHS3-3-body uses BasisFunctionsAngular.
     * All three kinds are stored in the map basisFunctions.
     **********************************************************************************************/
    void update(const BasisFunctionMap& basisFunctions, const std::shared_ptr<TypeMap>& typeMap);
    //**********************************************************************************************
    // GETTERS
    //**********************************************************************************************
    const DescriptorMap&           get_descriptors() const;
    std::shared_ptr<const TypeMap> get_typeMap() const;

  protected:
    /***********************************************************************************************
     * Multi-type map collecting all memory required for this Frame and its member instances of
     * NeighborList and Descriptor.
     *
     * This will be set up in the init() function by createFrameMemory().
     **********************************************************************************************/
    MultiTypeMap data;
    /// A single structure for which neighbor lists and descriptors are calculated in this Frame.
    std::shared_ptr<Structure> structure;
    /***********************************************************************************************
     * Map of neighbor lists corresponding to different descriptors, cutoffs, etc.
     *
     * This will be set up in the init() function by createNeighborLists().
     **********************************************************************************************/
    NeighborListMap neighborLists;
    /***********************************************************************************************
     * Map of descriptors required for force field evaluation.
     *
     * This will be set up in the init() function by createDescriptors().
     **********************************************************************************************/
    DescriptorMap descriptors;
    /// Map from force field types to types in structure.
    std::shared_ptr<TypeMap> typeMap;
    /// List of types present in force field.
    ShCVec1String forceFieldTypes;
    /// List of types present in structure.
    ShVec1String structureTypes;
    /// Variable selecting the descriptor type for the 3-body descriptor
    Descriptor3BodyType descriptor3BodyType;
    /*******************************************************************************************
     * determines the execution policy for the computaionally demanding algorithms of class
     *******************************************************************************************/
    std::unique_ptr<const ExecutionPolicy> algoExecution;
};

/***************************************************************************************************
 * Generate dictionary entries for default VASP MLFF.
 *
 * @param inputParameters Instance of force field setup.
 *
 * @return A multi-type map with keys neccessary for default force field memory layout.
 **************************************************************************************************/
MultiTypeMap createFrameMemory(const IoHandlerML_FF& inputParameters);
/***************************************************************************************************
 * Generate dictionary with all neighbor lists required by a given force field setup.
 *
 * @param inputParameters Instance of force field setup.
 * @param arrayData Frame memory (multi-type map) with pre-defined keys for neighbor list memory.
 *
 * @return Map with all required instances of neighbor lists.
 *
 * @pre Currently createFrameMemory() needs to be called in advance to create the desired multi-type
 * map entries.
 **************************************************************************************************/
NeighborListMap createNeighborLists(const IoHandlerML_FF& inputParameters,
                                    const MultiTypeMap&   arrayData);
/***************************************************************************************************
 * Generate dictionary with all descriptors required by a given force field setup.
 *
 * @param inputParameters Instance of force field setup.
 * @param arrayData Frame memory (multi-type map) with pre-defined keys for descriptor memory.
 * @param basis3Body Instance of BasisFunctionsAngular (?) class.
 * @param descriptor3BodyType Enumeration of 3-body descriptor type.
 *
 * @return Map with all required instances of descriptors.
 *
 * @pre Currently createFrameMemory() needs to be called in advance to create the desired multi-type
 * map entries.
 *
 * @todo Get rid of basis3Body argument!
 **************************************************************************************************/
DescriptorMap createDescriptors(const IoHandlerML_FF&      inputParameters,
                                const MultiTypeMap&        arrayData,
                                const BasisFunctions&      basis3Body,
                                const Descriptor3BodyType& descriptor3BodyType,
                                const ExecutionPolicy      algoExecution);

} // namespace vaspml

#endif
