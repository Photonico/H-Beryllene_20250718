#include "Frame.hpp"
#include "DescriptorSHS2.hpp"
#include "DescriptorSHS3.hpp"
#include "DescriptorSHS3ReducedLinElem.hpp"
#include "constants.hpp"
#include "utils.hpp"

using namespace vaspml;

namespace vaspml
{

MultiTypeMap createFrameMemory(const IoHandlerML_FF& /*inputParameters*/)
{
    MultiTypeMap data;

    // clang-format off
    // poscar data
    data["lattice"]      = std::make_shared<Vec1Real>();
    data["positions"]    = std::make_shared<Vec1Real>();
    data["atom_types"]   = std::make_shared<Vec1String>();
    data["number_atoms"] = std::make_shared<Vec1Int>();
    data["type_index"]   = std::make_shared<Vec1Int>();

    // nearest neighbor data
    data["2-body-n_globalIndex"]                = std::make_shared<Vec2Int>();
    data["2-body-n_typeIndex"]                  = std::make_shared<Vec2Int>();
    data["2-body-n_typeIndexCentral"]           = std::make_shared<Vec1Int>();
    data["2-body-n_distances"]                  = std::make_shared<Vec2Real>();
    data["2-body-n_connectionVector"]           = std::make_shared<Vec2Real>();
    data["2-body-n_connectionVectorNormalized"] = std::make_shared<Vec2Real>();
    data["2-body-n_numberNeighbors"]            = std::make_shared<Vec1Int>();
    data["2-body-n_numberNeighborsType"]        = std::make_shared<Vec2Int>();
    data["2-body-n_nAtomsType"]                 = std::make_shared<Vec1Int>();
    data["2-body-n_centralAtomIndexPerType"]    = std::make_shared<Vec1Int>();

    data["3-body-n_globalIndex"]                = std::make_shared<Vec2Int>();
    data["3-body-n_typeIndex"]                  = std::make_shared<Vec2Int>();
    data["3-body-n_typeIndexCentral"]           = std::make_shared<Vec1Int>();
    data["3-body-n_distances"]                  = std::make_shared<Vec2Real>();
    data["3-body-n_connectionVector"]           = std::make_shared<Vec2Real>();
    data["3-body-n_connectionVectorNormalized"] = std::make_shared<Vec2Real>();
    data["3-body-n_numberNeighbors"]            = std::make_shared<Vec1Int>();
    data["3-body-n_numberNeighborsType"]        = std::make_shared<Vec2Int>();
    data["3-body-n_nAtomsType"]                 = std::make_shared<Vec1Int>();
    data["3-body-n_centralAtomIndexPerType"]    = std::make_shared<Vec1Int>();

    // 2-body terms
    data["2-body-SHS2-pair"]                    = std::make_shared<Vec2Real>();
    data["2-body-SHS2-pair-derivative"]         = std::make_shared<Vec2Real>();
    data["2-body-SHS-vasp"]                     = std::make_shared<Vec2Real>();
    data["2-body-derivative-SHS2-central-vasp"] = std::make_shared<Vec2Real>();

    // 2-body terms to create 3-body terms
    data["3-body-SHS2-pair"]                    = std::make_shared<Vec2Real>();
    data["3-body-SHS2-pair-derivative"]         = std::make_shared<Vec2Real>();
    data["3-body-SHS-vasp"]                     = std::make_shared<Vec2Real>();
    data["3-body-derivative-SHS2-central-vasp"] = std::make_shared<Vec2Real>();

    // 3-body terms
    data["3-body-SHS3-spectrum"]      = std::make_shared<Vec2Real>();
    data["3-body-SHS3-pair-spectrum"] = std::make_shared<Vec2Real>();
    // clang-format on

    return data;
}

NeighborListMap createNeighborLists(const IoHandlerML_FF& inputParameters,
                                    const MultiTypeMap&   arrayData)
{
    NeighborListMap neighborLists;

    neighborLists["2-body"] = std::make_shared<NearestNeighborNSquare>(
        inputParameters["SHS2-2-body-cutoff"].cget<Real>(),
        true,
        false,
        std::get<ShVec2Int>(arrayData.at("2-body-n_globalIndex")),
        std::get<ShVec2Int>(arrayData.at("2-body-n_typeIndex")),
        std::get<ShVec1Int>(arrayData.at("2-body-n_typeIndexCentral")),
        std::get<ShVec2Real>(arrayData.at("2-body-n_distances")),
        std::get<ShVec2Real>(arrayData.at("2-body-n_connectionVector")),
        std::get<ShVec2Real>(arrayData.at("2-body-n_connectionVectorNormalized")),
        std::get<ShVec1Int>(arrayData.at("2-body-n_numberNeighbors")),
        std::get<ShVec2Int>(arrayData.at("2-body-n_numberNeighborsType")),
        std::get<ShVec1Int>(arrayData.at("2-body-n_nAtomsType")),
        std::get<ShVec1Int>(arrayData.at("2-body-n_centralAtomIndexPerType")));

    neighborLists["3-body"] = std::make_shared<NearestNeighborNSquare>(
        inputParameters["SHS3-3-body-cutoff"].cget<Real>(),
        true,
        false,
        std::get<ShVec2Int>(arrayData.at("3-body-n_globalIndex")),
        std::get<ShVec2Int>(arrayData.at("3-body-n_typeIndex")),
        std::get<ShVec1Int>(arrayData.at("3-body-n_typeIndexCentral")),
        std::get<ShVec2Real>(arrayData.at("3-body-n_distances")),
        std::get<ShVec2Real>(arrayData.at("3-body-n_connectionVector")),
        std::get<ShVec2Real>(arrayData.at("3-body-n_connectionVectorNormalized")),
        std::get<ShVec1Int>(arrayData.at("3-body-n_numberNeighbors")),
        std::get<ShVec2Int>(arrayData.at("3-body-n_numberNeighborsType")),
        std::get<ShVec1Int>(arrayData.at("3-body-n_nAtomsType")),
        std::get<ShVec1Int>(arrayData.at("3-body-n_centralAtomIndexPerType")));

    return neighborLists;
}

DescriptorMap createDescriptors(const IoHandlerML_FF&      inputParameters,
                                const MultiTypeMap&        arrayData,
                                const BasisFunctions&      basis3Body,
                                const Descriptor3BodyType& descriptor3BodyType,
                                const ExecutionPolicy      algoExecution)
{
    DescriptorMap descriptors;

    descriptors["SHS2-2-body"] = std::make_shared<DescriptorSHS2>(
        inputParameters["SHS2-2-body-weight"].cget<Real>(),
        inputParameters["SHS2-2-body-is-normalized"].cget<bool>(),
        std::get<ShVec2Real>(arrayData.at("2-body-SHS2-pair")),
        std::get<ShVec2Real>(arrayData.at("2-body-SHS2-pair-derivative")),
        std::get<ShVec2Real>(arrayData.at("2-body-SHS-vasp")),
        std::get<ShVec2Real>(arrayData.at("2-body-derivative-SHS2-central-vasp")),
        DescriptorType::bodyOrder2,
        algoExecution);

    descriptors["SHS2-3-body"] = std::make_shared<DescriptorSHS2>(
        inputParameters["SHS3-3-body-weight"].cget<Real>(),
        inputParameters["SHS3-3-body-is-normalized"].cget<bool>(),
        std::get<ShVec2Real>(arrayData.at("3-body-SHS2-pair")),
        std::get<ShVec2Real>(arrayData.at("3-body-SHS2-pair-derivative")),
        std::get<ShVec2Real>(arrayData.at("3-body-SHS-vasp")),
        std::get<ShVec2Real>(arrayData.at("3-body-derivative-SHS2-central-vasp")),
        DescriptorType::bodyOrder2,
        algoExecution);

    if (descriptor3BodyType == Descriptor3BodyType::StandardDescriptor)
    {
        descriptors["SHS3-3-body"] = std::make_shared<DescriptorSHS3>(
            inputParameters["SHS3-3-body-descriptor-list"].dcget<ShVec2Int>(),
            inputParameters["SHS3-3-body-angular-filter-on"].cget<bool>(),
            inputParameters["SHS3-3-body-angular-filter-type"].cget<Int>(),
            inputParameters["SHS3-3-body-angular-filter-scale"].cget<Real>(),
            basis3Body.get_nValidRoots(),
            inputParameters["SHS3-3-body-max-angular-number"].cget<Int>(),
            inputParameters["SHS3-3-body-weight"].cget<Real>(),
            inputParameters["SHS3-3-body-is-normalized"].cget<bool>(),
            std::get<ShVec2Real>(arrayData.at("3-body-SHS3-spectrum")),
            std::get<ShVec2Real>(arrayData.at("3-body-SHS3-pair-spectrum")),
            DescriptorType::bodyOrder3,
            algoExecution);
    }
    else if (descriptor3BodyType == Descriptor3BodyType::LinearScalingDescriptor)
    {
        descriptors["SHS3-3-body"] = std::make_shared<DescriptorSHS3ReducedLinElem>(
            inputParameters["SHS3-3-body-descriptor-list"].dcget<ShVec2Int>(),
            inputParameters["SHS3-3-body-angular-filter-on"].cget<bool>(),
            inputParameters["SHS3-3-body-angular-filter-type"].cget<Int>(),
            inputParameters["SHS3-3-body-angular-filter-scale"].cget<Real>(),
            basis3Body.get_nValidRoots(),
            inputParameters["SHS3-3-body-max-angular-number"].cget<Int>(),
            inputParameters["SHS3-3-body-weight"].cget<Real>(),
            inputParameters["SHS3-3-body-is-normalized"].cget<bool>(),
            std::get<ShVec2Real>(arrayData.at("3-body-SHS3-spectrum")),
            std::get<ShVec2Real>(arrayData.at("3-body-SHS3-pair-spectrum")),
            DescriptorType::bodyOrder3LinearElement,
            algoExecution);
    }

    return descriptors;
}

} // namespace vaspml

void Frame::init(const IoHandlerML_FF& inputParameters,
                 const BasisFunctions& basis3Body,
                 const ExecutionPolicy algoExecution)
{
    data = createFrameMemory(inputParameters);
    neighborLists = createNeighborLists(inputParameters, data);
    this->algoExecution = std::make_unique<ExecutionPolicy>(algoExecution);

    structure = std::make_shared<Structure>(std::get<ShVec1Real>(data["lattice"]),
                                            std::get<ShVec1Real>(data["positions"]),
                                            std::get<ShVec1Int>(data["type_index"]),
                                            std::get<ShVec1String>(data["atom_types"]),
                                            std::get<ShVec1Int>(data["number_atoms"]));

    forceFieldTypes = inputParameters["types"].cget<ShVec1String>();

    switch (inputParameters["descriptor-type"].cget<Int>())
    {
    case 0:
        this->descriptor3BodyType = Descriptor3BodyType::StandardDescriptor;
        break;
    case 1:
        this->descriptor3BodyType = Descriptor3BodyType::LinearScalingDescriptor;
        break;
    }
    descriptors = createDescriptors(inputParameters,
                                    data,
                                    basis3Body,
                                    descriptor3BodyType,
                                    *this->algoExecution);
    return;
}

void Frame::set_basisFunctions(BasisFunctionMap& basisFunctions)
{

    for (const auto& x : neighborLists)
    {
        String          descKey = "SHS2-" + x.first;
        DescriptorSHS2* desc = static_cast<DescriptorSHS2*>(descriptors[descKey].get());
        desc->set_basisFunctions(basisFunctions.at(x.first));
    }
}

void Frame::update(const BasisFunctionMap& basisFunctions)
{
    // compute SHS2 descriptors
    VASPML_PROFILING_START("DescriptorSHS2");
    for (auto& [key, item] : neighborLists)
    {
        String descKey = "SHS2-" + key;
        item->compute_centralAtomIndexPerType();
        DescriptorSHS2* desc = static_cast<DescriptorSHS2*>(descriptors[descKey].get());
        desc->updatePairCoefficients(item, basisFunctions.at(key));
        desc->computeVaspCoefficientsFromPairCoefficients();
        // desc->write_clnmVasp(descKey + ".dat");
    }
    VASPML_PROFILING_STOP("DescriptorSHS2");

    // compute SHS3 and SHS3LinElem descriptor
    VASPML_PROFILING_START("DescriptorSHS3");
    std::shared_ptr<DescriptorSHS2> desc2 =
        std::static_pointer_cast<DescriptorSHS2>(descriptors["SHS2-3-body"]);
    if (descriptor3BodyType == Descriptor3BodyType::StandardDescriptor)
    {
        DescriptorSHS3* desc3 = static_cast<DescriptorSHS3*>(descriptors["SHS3-3-body"].get());
        desc3->computeSHS3(desc2, *typeMap);
    }
    else if (descriptor3BodyType == Descriptor3BodyType::LinearScalingDescriptor)
    {
        DescriptorSHS3ReducedLinElem* desc3 =
            static_cast<DescriptorSHS3ReducedLinElem*>(descriptors["SHS3-3-body"].get());
        desc3->computeSHS3(desc2, *typeMap);
    }
    VASPML_PROFILING_STOP("DescriptorSHS3");

    return;
}

void Frame::update(void)
{

    VASPML_PROFILING_START("DescriptorSHS2");
    for (auto& [key, item] : neighborLists)
    {
        String descKey = "SHS2-" + key;
        item->compute_centralAtomIndexPerType();
        DescriptorSHS2* desc = static_cast<DescriptorSHS2*>(descriptors[descKey].get());
        desc->updateVaspCoefficients(item);
        // desc->write_clnmVasp(descKey + ".dat");
    }
    VASPML_PROFILING_STOP("DescriptorSHS2");

    // compute SHS3 and SHS3LinElem descriptor
    VASPML_PROFILING_START("DescriptorSHS3");
    std::shared_ptr<DescriptorSHS2> desc2 =
        std::static_pointer_cast<DescriptorSHS2>(descriptors["SHS2-3-body"]);
    if (descriptor3BodyType == Descriptor3BodyType::StandardDescriptor)
    {
        DescriptorSHS3* desc3 = static_cast<DescriptorSHS3*>(descriptors["SHS3-3-body"].get());
        desc3->computeSHS3(desc2, *typeMap);
    }
    else if (descriptor3BodyType == Descriptor3BodyType::LinearScalingDescriptor)
    {
        DescriptorSHS3ReducedLinElem* desc3 =
            static_cast<DescriptorSHS3ReducedLinElem*>(descriptors["SHS3-3-body"].get());
        desc3->computeSHS3(desc2, *typeMap);
    }
    VASPML_PROFILING_STOP("DescriptorSHS3");

    return;
}

void Frame::update(const String& fileName, const BasisFunctionMap& basisFunctions)
{

    structure->readPoscar(fileName);
    structure->cartesianToDirect();
    // structure->position_direct_to_cartesian();
    structureTypes = structure->get_typeNames();

    typeMap = std::make_shared<TypeMap>(*forceFieldTypes, *structureTypes);
    // update neighbor lists
    for (auto& list : neighborLists)
    {
        list.second->computeNearestNeighborsDirectCoordinates(*structure);
        // list.second->computeNearestNeighborsCartesianCoordinates(*structure);
        // list.second->writeListToScreen();
    }

    update(basisFunctions);

    return;
}

void Frame::update(const BasisFunctionMap& basisFunctions, const std::shared_ptr<TypeMap>& typeMap)
{
    this->typeMap = typeMap;

    update(basisFunctions);

    return;
}

const DescriptorMap& Frame::get_descriptors(void) const
{
    return descriptors;
}

std::shared_ptr<const TypeMap> Frame::get_typeMap(void) const
{
    return typeMap;
}
