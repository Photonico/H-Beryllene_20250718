#include "InterfaceVASP.hpp"

#include "BasisFunctionsAngular.hpp"
#include "BasisFunctionsRadialSpline.hpp"
#include "ParallelEnvironemt.hpp"
#include "Timer.hpp"
#include "constants.hpp"
#include "utils.hpp"

#include <memory>

//debug
#include <iomanip>

using namespace vaspml;

extern "C"
{

void* createFrame()
{
    return new InterfaceVASP();
}

void destroyFrame(void** ptrFrame)
{
    //InterfaceVASP* const& frame = static_cast<InterfaceVASP*>(*ptrFrame);
    //if (frame->mpiMain->get_rank() == 0) { VASPML_PROFILING_WRITE(); }

    delete static_cast<InterfaceVASP*>(*ptrFrame);
    *ptrFrame = nullptr;
    return;
}

Real get_W1(void** ptrFrame)
{
    InterfaceVASP* const& frame = static_cast<InterfaceVASP*>(*ptrFrame);
    return frame->forceFieldData["SHS2-2-body-weight"].cget<Real>();
}

Real get_W2(void** ptrFrame)
{
    InterfaceVASP* const& frame = static_cast<InterfaceVASP*>(*ptrFrame);
    return frame->forceFieldData["SHS3-3-body-weight"].cget<Real>();
}

Real get_RCUT1(void** ptrFrame)
{
    InterfaceVASP* const& frame = static_cast<InterfaceVASP*>(*ptrFrame);
    return frame->forceFieldData["SHS2-2-body-cutoff"].cget<Real>();
}

Real get_RCUT2(void** ptrFrame)
{
    InterfaceVASP* const& frame = static_cast<InterfaceVASP*>(*ptrFrame);
    return frame->forceFieldData["SHS3-3-body-cutoff"].cget<Real>();
}

void resizeNeighborArrays(void** ptrFrame, const Int* nions, const char* keyIn)
{

    InterfaceVASP* const& frame = static_cast<InterfaceVASP*>(*ptrFrame);
    String                key = keyIn;
    frame->resizeNeighborArrays(*nions, key);
}

void set_nAtomsType(void** ptrFrame, const Int* ntypes, const Int* nAtomsType, const char* keyIn)
{
    InterfaceVASP* const& frame = static_cast<InterfaceVASP*>(*ptrFrame);
    String                key = keyIn;
    frame->set_nAtomsType(*ntypes, nAtomsType, key);
}

void fillNeighhborArrays(void**      ptrFrame,
                         const Int*  numberNeighbors,
                         const Int*  atomNumber,
                         const Int*  centralType,
                         const Int*  neighborIndex,
                         const Int*  neighborTypes,
                         const Real* neighborDist,
                         const Real* neighborConnect,
                         const char* keyIn)
{

    String                key = keyIn;
    InterfaceVASP* const& frame = static_cast<InterfaceVASP*>(*ptrFrame);
    frame->fillNeighhborArrays(*numberNeighbors,
                               *atomNumber,
                               *centralType,
                               neighborIndex,
                               neighborTypes,
                               neighborDist,
                               neighborConnect,
                               key);
}

void setupForceField(void** ptrFrame, MPI_Fint* mpiComm)
{
    InterfaceVASP* const& frame = static_cast<InterfaceVASP*>(*ptrFrame);

    using CT = vaspml::math::CutoffType;

    if (mpiComm != nullptr) frame->mpiMain = std::make_shared<MlMPI>(MPI_Comm_f2c(*mpiComm), false);

    IoHandlerML_FF& ffd = frame->forceFieldData;
    ffd.init("ML_FF", frame->mpiMain);
    ffd.mainReader();
    // convert units to atomic units. needed because of unit problem in old code
    ffd.convertUnitsToAtomicUnits();

    frame->basisFunctions["2-body"] = std::make_shared<BasisFunctionsRadialSpline>(
        static_cast<CT>(std::get<Int>(ffd["SHS2-2-body-cutoff-type"])),
        std::get<Real>(ffd["SHS2-2-body-cutoff"]),
        std::get<Real>(ffd["SHS2-2-body-smearing-param"]),
        std::get<Int>(ffd["SHS2-2-body-number-radial-grid-points"]),
        0,
        std::get<Int>(ffd["SHS2-2-body-max-radial-funcs"]),
        BasisFunctionType::bodyOrder2);

    frame->basisFunctions["3-body"] = std::make_shared<BasisFunctionsAngular>(
        static_cast<CT>(std::get<Int>(ffd["SHS3-3-body-cutoff-type"])),
        std::get<Real>(ffd["SHS3-3-body-cutoff"]),
        std::get<Real>(ffd["SHS3-3-body-smearing-param"]),
        std::get<Int>(ffd["SHS3-3-body-number-radial-grid-points"]),
        std::get<Int>(ffd["SHS3-3-body-max-angular-number"]),
        std::get<Int>(ffd["SHS3-3-body-max-radial-funcs"]),
        BasisFunctionType::bodyOrder3);
    // init from Frame
    //frame->init(ffd, *(frame->basisFunctions["3-body"]), ExecutionPolicy::gpuStdLib);
    frame->init(ffd, *(frame->basisFunctions["3-body"]), ExecutionPolicy::cpuSingleCore);
    frame->typeMapPtr = std::make_shared<TypeMap>();
    // call other init functions
    frame->polynomialKernel.init(ffd, frame->mpiMain);
    frame->predictor.init(ffd, frame->mpiMain);
    return;
}

void update(void** ptrFrame, const Real* volume)
{

    VASPML_PROFILING_START("UPDATE");
    InterfaceVASP* const& frame = static_cast<InterfaceVASP*>(*ptrFrame);
    frame->update(frame->basisFunctions, frame->typeMapPtr);
    frame->polynomialKernel.updatePolynomialKernel(*frame);
    frame->predictor.update(frame->polynomialKernel);
    frame->predictor.compute_stressTensor(frame->polynomialKernel.get_descriptorCollection(),
                                          *volume);
    VASPML_PROFILING_STOP("UPDATE");
}

void set_typeMap(void** ptrFrame, const char* typesIn)
{

    String     types = typesIn;
    Vec1String typeElements = string_tools::splitString(types, ";");
    if (string_tools::trim(typeElements.back()).empty()) typeElements.pop_back();

    InterfaceVASP* const& frame = static_cast<InterfaceVASP*>(*ptrFrame);
    frame->set_typeMap(typeElements);
}

void fillForceSingleAtom(void**     ptrFrame,
                         const Int* /* nions */,
                         const Int* centralVasp,
                         const Int* element,
                         Real*      force)
{

    InterfaceVASP* const& frame = static_cast<InterfaceVASP*>(*ptrFrame);
    const auto&           descriptorCollection = frame->polynomialKernel.get_descriptorCollection();
    Int                   centralAtom = *element;
    for (const String& key : constants::descriptorKeyList)
    {
        std::shared_ptr<const NearestNeighborNSquare> neighborList =
            descriptorCollection.getDescriptor(key).get_neighborList_ptr();
        force[3 * (*centralVasp)] += frame->predictor.get_centralForcesX(key, centralAtom);
        force[3 * (*centralVasp) + 1] += frame->predictor.get_centralForcesY(key, centralAtom);
        force[3 * (*centralVasp) + 2] += frame->predictor.get_centralForcesZ(key, centralAtom);
        for (std::size_t neighborAtom = 0; neighborAtom < neighborList->get_size(centralAtom);
             neighborAtom++)
        {
            const Int& neighborIndex = neighborList->get_globalIndex(centralAtom, neighborAtom);
            force[3 * neighborIndex] +=
                frame->predictor.get_pairForcesX(key, centralAtom, neighborAtom);
            force[3 * neighborIndex + 1] +=
                frame->predictor.get_pairForcesY(key, centralAtom, neighborAtom);
            force[3 * neighborIndex + 2] +=
                frame->predictor.get_pairForcesZ(key, centralAtom, neighborAtom);
        }
    }
}

void getStressTensor(void** ptrFrame, Real* stressTensor)
{

    InterfaceVASP* const& frame = static_cast<InterfaceVASP*>(*ptrFrame);
    stressTensor[0] = frame->predictor.get_totalStressTensor(0, 0);
    stressTensor[1] = frame->predictor.get_totalStressTensor(0, 1);
    stressTensor[2] = frame->predictor.get_totalStressTensor(0, 2);

    stressTensor[3] = frame->predictor.get_totalStressTensor(1, 0);
    stressTensor[4] = frame->predictor.get_totalStressTensor(1, 1);
    stressTensor[5] = frame->predictor.get_totalStressTensor(1, 2);

    stressTensor[6] = frame->predictor.get_totalStressTensor(2, 0);
    stressTensor[7] = frame->predictor.get_totalStressTensor(2, 1);
    stressTensor[8] = frame->predictor.get_totalStressTensor(2, 2);
}

void getPotentialEnergy(void** ptrFrame, Real* totalEnergy)
{

    InterfaceVASP* const& frame = static_cast<InterfaceVASP*>(*ptrFrame);
    (*totalEnergy) += frame->predictor.get_totalEnergy();
}

} // extern "C"

void InterfaceVASP::resizeNeighborArrays(const Int nions, const String& key)
{
    data[key + "-n_globalIndex"].dget<ShVec2Int>().resize(nions);
    data[key + "-n_typeIndex"].dget<ShVec2Int>().resize(nions);
    data[key + "-n_typeIndexCentral"].dget<ShVec1Int>().resize(nions);
    data[key + "-n_distances"].dget<ShVec2Real>().resize(nions);
    data[key + "-n_connectionVector"].dget<ShVec2Real>().resize(nions);
    data[key + "-n_connectionVectorNormalized"].dget<ShVec2Real>().resize(nions);
    data[key + "-n_numberNeighbors"].dget<ShVec1Int>().resize(nions);
    data[key + "-n_numberNeighborsType"].dget<ShVec2Int>().resize(nions);
    data[key + "-n_centralAtomIndexPerType"].dget<ShVec1Int>().resize(nions);
    neighborLists[key]->set_nAtoms(nions);
}

void InterfaceVASP::set_nAtomsType(const Int ntypes, const Int* nAtomsType, const String& key)
{

    neighborLists[key]->set_nTypes(ntypes);
    data[key + "-n_nAtomsType"].dget<ShVec1Int>().resize(ntypes);
    for (std::size_t type = 0; type < (std::size_t)ntypes; type++)
    {
        data[key + "-n_nAtomsType"].dget<ShVec1Int>()[type] = nAtomsType[type];
    }
}

void InterfaceVASP::fillNeighhborArrays(const Int     numberNeighbors,
                                        const Int     atomNumber,
                                        const Int     centralType,
                                        const Int*    neighborIndex,
                                        const Int*    neighborTypes,
                                        const Real*   neighborDist,
                                        const Real*   neighborConnect,
                                        const String& key)
{
    data[key + "-n_globalIndex"].dget<ShVec2Int>()[atomNumber].resize(numberNeighbors);
    data[key + "-n_typeIndex"].dget<ShVec2Int>()[atomNumber].resize(numberNeighbors);
    data[key + "-n_distances"].dget<ShVec2Real>()[atomNumber].resize(numberNeighbors);
    data[key + "-n_connectionVector"].dget<ShVec2Real>()[atomNumber].resize(3 * numberNeighbors);
    data[key + "-n_connectionVectorNormalized"].dget<ShVec2Real>()[atomNumber].resize(
        3 * numberNeighbors);

    data[key + "-n_numberNeighborsType"].dget<ShVec2Int>()[atomNumber].resize(
        data["2-body-n_nAtomsType"].dget<ShVec1Int>().size());
    std::fill(data[key + "-n_numberNeighborsType"].dget<ShVec2Int>()[atomNumber].begin(),
              data[key + "-n_numberNeighborsType"].dget<ShVec2Int>()[atomNumber].end(),
              (Real)0);

    data[key + "-n_typeIndexCentral"].dget<ShVec1Int>()[atomNumber] = centralType;
    data[key + "-n_numberNeighbors"].dget<ShVec1Int>()[atomNumber] = numberNeighbors;

    Vec1Int&  index = data[key + "-n_globalIndex"].dget<ShVec2Int>()[atomNumber];
    Vec1Int&  nTypes = data[key + "-n_typeIndex"].dget<ShVec2Int>()[atomNumber];
    Vec1Real& nDistance = data[key + "-n_distances"].dget<ShVec2Real>()[atomNumber];
    Vec1Real& nVectorNorm =
        data[key + "-n_connectionVectorNormalized"].dget<ShVec2Real>()[atomNumber];
    Vec1Real& nVector = data[key + "-n_connectionVector"].dget<ShVec2Real>()[atomNumber];
    Vec1Int&  nTypeCounter = data[key + "-n_numberNeighborsType"].dget<ShVec2Int>()[atomNumber];

    for (std::size_t nIndex = 0; nIndex < (std::size_t)numberNeighbors; nIndex++)
    {
        index[nIndex] = neighborIndex[nIndex];
        nTypes[nIndex] = neighborTypes[nIndex];
        nDistance[nIndex] = neighborDist[nIndex];
        nVectorNorm[3 * nIndex] = neighborConnect[3 * nIndex];
        nVectorNorm[3 * nIndex + 1] = neighborConnect[3 * nIndex + 1];
        nVectorNorm[3 * nIndex + 2] = neighborConnect[3 * nIndex + 2];
        nVector[3 * nIndex] = neighborConnect[3 * nIndex] * neighborDist[nIndex];
        nVector[3 * nIndex + 1] = neighborConnect[3 * nIndex + 1] * neighborDist[nIndex];
        nVector[3 * nIndex + 2] = neighborConnect[3 * nIndex + 2] * neighborDist[nIndex];
        // count number of neighors per type ( not really needed )
        nTypeCounter[neighborTypes[nIndex]]++;
    }
}

void InterfaceVASP::set_typeMap(const Vec1String& types)
{
    Vec1String ffTypes = forceFieldData["types"].dcget<ShVec1String>();
    typeMapPtr->update(ffTypes, types);
}
