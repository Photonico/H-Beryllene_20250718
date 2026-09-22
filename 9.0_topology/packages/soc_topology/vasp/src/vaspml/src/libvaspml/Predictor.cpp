#include "Predictor.hpp"
#include "ParallelEnvironemt.hpp"
#include "TagTranslator.hpp"
#include "Tutor.hpp"
#include "constants.hpp"
#include "debug.hpp"
#include "nearest_neighbor.hpp"
#include "utils.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

using namespace vaspml;

namespace vaspml
{
std::map<String, std::shared_ptr<ShmemArray2DVariableLen<Real>>> make_descriptorsRefConfsWeighted(
    const std::map<String, ShVec1Int>& featureSpaceSize,
    const std::map<String, Real>&      weights,
    const std::shared_ptr<MlMPI>&      mpiIn)
{
    std::map<String, std::shared_ptr<ShmemArray2DVariableLen<Real>>> arrays;
    for (const auto& [key, item] : featureSpaceSize)
    {
        arrays[key] = std::make_shared<ShmemArray2DVariableLen<Real>>(
            weights.at(key) > 0 ? *featureSpaceSize.at(key) : Vec1Int{0},
            mpiIn);
    }
    return arrays;
}

std::map<String, ShVec2Real> makeMaps2D(void)
{

    std::map<String, ShVec2Real> data;
    for (const String& key : constants::descriptorKeyList)
    {
        data[key] = std::make_shared<Vec2Real>();
    }
    return data;
}

std::map<String, ShVec1Real> makeMaps1D(void)
{
    std::map<String, ShVec1Real> data;
    for (const String& key : constants::descriptorKeyList)
    {
        data[key] = std::make_shared<Vec1Real>();
    }
    return data;
}

std::map<String, Vec1Size_t> makeMaps1DSize_t(void)
{
    std::map<String, Vec1Size_t> data;
    for (const String& key : constants::descriptorKeyList) data[key] = Vec1Size_t();
    return data;
}

std::map<String, ArrayResizing1D> makeResizeMap1D(void)
{
    std::map<String, ArrayResizing1D> data;
    for (const String& key : constants::descriptorKeyList) data[key] = ArrayResizing1D();
    return data;
}

std::map<String, ArrayResizing2D> makeResizeMap2D(void)
{
    std::map<String, ArrayResizing2D> data;
    for (const String& key : constants::descriptorKeyList) data[key] = ArrayResizing2D();
    return data;
}

std::shared_ptr<ShmemArray2DVariableLen<Real>> make_regressionCoefficients(
    const IoHandlerML_FF&         inputParameters,
    const std::shared_ptr<MlMPI>& mpiIn)
{

    const Vec2Real& regCoeff = inputParameters["regression-coeff"].dcget<ShVec2Real>();
    const Vec1Int&  sizes = inputParameters["number-local-reference-configs"].dcget<ShVec1Int>();
    std::shared_ptr<ShmemArray2DVariableLen<Real>> shmemPtr =
        std::make_shared<ShmemArray2DVariableLen<Real>>(sizes, mpiIn);
    shmemPtr->set_value(regCoeff);
    return shmemPtr;
}

} //namespace vaspml

template<>
SmartEnum<TotalEnergyType>::EnumMap SmartEnum<TotalEnergyType>::names = {
    {TotalEnergyType::IsolatedAtom,       "IsolatedAtom"         },
    {TotalEnergyType::AverageTrainEnergy, "AverageTrainingEnergy"},
};

Predictor::Predictor(const IoHandlerML_FF&         inputParameters,
                     const std::shared_ptr<MlMPI>& mpiIn,
                     const ShVec2Real&             energyArray,
                     const ShVec1Real&             atomicForces,
                     const ShVec1Real&             totalStressTensor,
                     const ExecutionPolicy&        algoExecution) :
    numberTypesForceField(inputParameters["number-types"].cget<Int>()),
    numberLocalRefConfs(std::make_shared<Vec1Int>(
        inputParameters["number-local-reference-configs"].dcget<ShVec1Int>())),
    numberDescriptors(inputParameters.get_numberDescriptorsMap()),
    featureSpaceSize(inputParameters.generate_featureSpaceMap()),
    weights(inputParameters.generate_weightsMap()),
    regressionCoefficients(make_regressionCoefficients(inputParameters, mpiIn)),
    referenceEnergyPerType(inputParameters["reference-energies"].cget<ShVec1Real>()),
    averageTrainEnergy(inputParameters["average-energy-per-atom"].cget<Real>()),
    descriptorsRefConfsWeighted(make_descriptorsRefConfsWeighted(featureSpaceSize, weights, mpiIn)),
    derivativeMatrix(makeMaps2D()),
    derivativeMatrixDim(makeMaps1DSize_t()),
    derivativeMatrixSize(makeResizeMap2D()),
    forcePreContract(makeMaps2D()),
    forcePreContractDim(makeMaps1DSize_t()),
    forcePreContractSize(makeResizeMap2D()),
    pairForces(makeMaps2D()),
    pairForcesDim(makeMaps1DSize_t()),
    pairForcesSize(makeResizeMap2D()),
    centralForces(makeMaps1D()),
    centralForcesSize(makeResizeMap1D()),
    stressTensor(makeMaps1D()),
    algoExecution(algoExecution)
{

    if (energyArray == nullptr) { this->energyArray = std::make_shared<Vec2Real>(); }
    else { this->energyArray = energyArray; }
    if (atomicForces == nullptr) { this->atomicForces = std::make_shared<Vec1Real>(); }
    else { this->atomicForces = atomicForces; }
    if (totalStressTensor == nullptr) { this->totalStressTensor = std::make_shared<Vec1Real>(); }
    else { this->totalStressTensor = totalStressTensor; }

    forceArraysComputed = false;

    compute_descriptorsRefConfsWeighted(inputParameters);

    switch (inputParameters["scale-type-total-energy"].cget<Int>())
    {
    case 1:
        this->energyType = TotalEnergyType::IsolatedAtom;
        break;
    case 2:
        this->energyType = TotalEnergyType::AverageTrainEnergy;
        break;
    }
}

void Predictor::init(const IoHandlerML_FF&         inputParameters,
                     const std::shared_ptr<MlMPI>& mpiIn,
                     const ShVec2Real&             energyArray,
                     const ShVec1Real&             atomicForces,
                     const ShVec1Real&             totalStressTensor,
                     const ExecutionPolicy&        algoExecution)
{

    numberTypesForceField = inputParameters["number-types"].cget<Int>();
    numberLocalRefConfs = std::make_shared<Vec1Int>(
        inputParameters["number-local-reference-configs"].dcget<ShVec1Int>());
    numberDescriptors = inputParameters.get_numberDescriptorsMap();
    featureSpaceSize = inputParameters.generate_featureSpaceMap();
    weights = inputParameters.generate_weightsMap();
    regressionCoefficients = make_regressionCoefficients(inputParameters, mpiIn);
    referenceEnergyPerType = inputParameters["reference-energies"].cget<ShVec1Real>();
    averageTrainEnergy = inputParameters["average-energy-per-atom"].cget<Real>();
    descriptorsRefConfsWeighted =
        make_descriptorsRefConfsWeighted(featureSpaceSize, weights, mpiIn);
    derivativeMatrix = makeMaps2D();
    derivativeMatrixDim = makeMaps1DSize_t();
    derivativeMatrixSize = makeResizeMap2D();
    forcePreContract = makeMaps2D();
    forcePreContractDim = makeMaps1DSize_t();
    forcePreContractSize = makeResizeMap2D();
    pairForces = makeMaps2D();
    pairForcesDim = makeMaps1DSize_t();
    pairForcesSize = makeResizeMap2D();
    centralForces = makeMaps1D();
    centralForcesSize = makeResizeMap1D();
    stressTensor = makeMaps1D();
    if (energyArray == nullptr) { this->energyArray = std::make_shared<Vec2Real>(); }
    else { this->energyArray = energyArray; }
    if (atomicForces == nullptr) { this->atomicForces = std::make_shared<Vec1Real>(); }
    else { this->atomicForces = atomicForces; }
    if (totalStressTensor == nullptr) { this->totalStressTensor = std::make_shared<Vec1Real>(); }
    else { this->totalStressTensor = totalStressTensor; }

    forceArraysComputed = false;
    this->algoExecution = algoExecution;

    compute_descriptorsRefConfsWeighted(inputParameters);

    switch (inputParameters["scale-type-total-energy"].cget<Int>())
    {
    case 1:
        this->energyType = TotalEnergyType::IsolatedAtom;
        break;
    case 2:
        this->energyType = TotalEnergyType::AverageTrainEnergy;
        break;
    }
}

void Predictor::compute_descriptorsRefConfsWeighted(const IoHandlerML_FF& inputParameters)
{

    for (std::string key : constants::descriptorKeyList)
    {
        std::string     tag = key + "-reference-configs";
        const Vec2Real& regCoeff = inputParameters["regression-coeff"].dcget<ShVec2Real>();
        if (weights[key] > 0)
        {
            const Vec2Real descriptor = inputParameters[tag].dcget<ShVec2Real>();
            for (std::size_t type = 0; type < descriptor.size(); type++)
            {
                std::size_t featureSpaceCounter = 0;
                for (std::size_t local_reference = 0;
                     local_reference < (std::size_t)(*numberLocalRefConfs)[type];
                     local_reference++)
                {
                    for (std::size_t desc_count = 0;
                         desc_count < (std::size_t)(*numberDescriptors[key])[type];
                         desc_count++)
                    {
#ifdef use_shmem
                        if ((*descriptorsRefConfsWeighted[key])["mpiShmem"].get_rank() == 0
                            and (*descriptorsRefConfsWeighted[key])["mpiInter"].get_rank() == 0)
#endif
                            descriptorsRefConfsWeighted[key]->set_value(
                                type,
                                featureSpaceCounter,
                                regCoeff[type][local_reference]
                                    * descriptor[type][featureSpaceCounter]);

                        featureSpaceCounter++;
                    }
                }
            }
        }
#ifdef use_shmem
        (*descriptorsRefConfsWeighted[key])["mpiShmem"].barrier();
        descriptorsRefConfsWeighted[key]->distributeInternode();
#endif
    }
}

void Predictor::compute_energyArrayDim(const Size_t numberTypesStruc)
{
    energyArrayDim.resize(numberTypesStruc);
    for (std::size_t typeStruc = 0; typeStruc < numberTypesStruc; typeStruc++)
    {
        std::size_t typeFF = typeMap->toType(typeStruc);
        energyArrayDim[typeStruc] = (*numberLocalRefConfs)[typeFF];
    }
}

void Predictor::allocate_energyArray(const Size_t numberTypesStruc)
{

    compute_energyArrayDim(numberTypesStruc);
    if (algoExecution == ExecutionPolicy::cpuSingleCore) allocate_energyArrayCPU(numberTypesStruc);
    else if (algoExecution == ExecutionPolicy::gpuStdLib) allocate_energyArrayGPU(numberTypesStruc);
}

void Predictor::allocate_energyArrayCPU(const std::size_t numberTypesStruc)
{
    energyArray->resize(numberTypesStruc);
    for (std::size_t typeStruc = 0; typeStruc < numberTypesStruc; typeStruc++)
    {
        (*energyArray)[typeStruc].resize(energyArrayDim[typeStruc]);
        std::fill((*energyArray)[typeStruc].begin(), (*energyArray)[typeStruc].end(), (Real)0);
    }
}

void Predictor::allocate_energyArrayGPU(const std::size_t numberTypesStruc)
{
    bool resize = energyArraySize.checkResize1Dim(numberTypesStruc);
    if (resize)
    {
        energyArraySize.resizeArray1Dim(*energyArray, energyArrayDim);
    } // second dimension can not change in execution mode. numberLocalRefConfs conserved
}

void Predictor::compute_derivativeMatrixDim(const ShVec1Int& numberAtomType)
{
    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        derivativeMatrixDim[key].resize(numberAtomType->size());
        for (Size_t typeStruc = 0; typeStruc < numberAtomType->size(); typeStruc++)
        {
            Size_t typeForceField = typeMap->toType(typeStruc);
            derivativeMatrixDim[key][typeStruc] =
                (*numberAtomType)[typeStruc] * (*numberDescriptors[key])[typeForceField];
        }
    }
}

void Predictor::allocate_derivativeMatrix(const ShVec1Int& numberAtomType)
{

    compute_derivativeMatrixDim(numberAtomType);
    if (algoExecution == ExecutionPolicy::cpuSingleCore)
        allocate_derivativeMatrixCPU(numberAtomType);
    else if (algoExecution == ExecutionPolicy::gpuStdLib)
        allocate_derivativeMatrixGPU(numberAtomType);
}

void Predictor::allocate_derivativeMatrixCPU(const ShVec1Int& numberAtomType)
{
    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        derivativeMatrix[key]->resize(numberAtomType->size());
        for (std::size_t typeStruc = 0; typeStruc < numberAtomType->size(); typeStruc++)
        {
            //std::size_t typeForceField  =  typeMap->toType( typeStruc );
            (*derivativeMatrix[key])[typeStruc].resize(derivativeMatrixDim[key][typeStruc]);
            std::fill((*derivativeMatrix[key])[typeStruc].begin(),
                      (*derivativeMatrix[key])[typeStruc].end(),
                      (Real)0);
        }
    }
}

void Predictor::allocate_derivativeMatrixGPU(const ShVec1Int& numberAtomType)
{
    bool resize;
    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        resize = derivativeMatrixSize[key].checkResize1Dim(numberAtomType->size());
        if (resize)
        {
            derivativeMatrixSize[key].resizeArray1Dim((*derivativeMatrix[key]),
                                                      derivativeMatrixDim[key]);
        }
        else
        {
            resize = derivativeMatrixSize[key].checkResize2Dim(derivativeMatrixDim[key]);
            if (resize)
            {
                derivativeMatrixSize[key].resizeArray2Dim((*derivativeMatrix[key]),
                                                          derivativeMatrixDim[key]);
            }
        }
        std::for_each(VASPML_PAR_UNSEQ derivativeMatrix[key]->begin(),
                      derivativeMatrix[key]->begin() + derivativeMatrixSize[key].act1Dim,
                      [&](Vec1Real& slice) { std::fill(slice.begin(), slice.end(), (Real)0); });
    }
}

void Predictor::compute_forcePreContractDim(const DescriptorCollector& descriptorCollection,
                                            const std::size_t          numberAtoms)
{
    for (const String& key : constants::descriptorKeyList)
    {
        forcePreContractDim[key].resize(numberAtoms);
        for (Size_t nAtom = 0; nAtom < numberAtoms; nAtom++)
        {
            const Descriptor& descriptor = descriptorCollection.getDescriptor(key);
            forcePreContractDim[key][nAtom] = descriptor.get_forcePreContractSize(nAtom);
        }
    }
}

void Predictor::allocate_forcePreContract(const DescriptorCollector& descriptorCollection,
                                          const std::size_t          numberAtoms)
{
    compute_forcePreContractDim(descriptorCollection, numberAtoms);
    if (algoExecution == ExecutionPolicy::cpuSingleCore) allocate_forcePreContractCPU(numberAtoms);
    else if (algoExecution == ExecutionPolicy::gpuStdLib) allocate_forcePreContractGPU(numberAtoms);
}

void Predictor::allocate_forcePreContractCPU(const std::size_t numberAtoms)
{
    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        forcePreContract[key]->resize(numberAtoms);
        for (Size_t nAtom = 0; nAtom < numberAtoms; nAtom++)
        {
            (*forcePreContract[key])[nAtom].resize(forcePreContractDim[key][nAtom]);
            std::fill((*forcePreContract[key])[nAtom].begin(),
                      (*forcePreContract[key])[nAtom].end(),
                      (Real)0);
        }
    }
}

void Predictor::allocate_forcePreContractGPU(const std::size_t numberAtoms)
{
    allocate_centralAtomIndex(numberAtoms);
    bool resize;
    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        // maps are not allowed in parallel region which is later needed in fill
        Vec2Real&        forcePreContract = *this->forcePreContract[key];
        ArrayResizing2D& forcePreContractSize = this->forcePreContractSize[key];
        resize = forcePreContractSize.checkResize1Dim(numberAtoms);
        if (resize)
        {
            forcePreContractSize.resizeArray1Dim(forcePreContract, forcePreContractDim[key]);
        }
        else
        {
            resize = forcePreContractSize.checkResize2Dim(forcePreContractDim[key]);
            if (resize)
            {
                forcePreContractSize.resizeArray2Dim(forcePreContract, forcePreContractDim[key]);
            }
        }
        std::for_each(VASPML_PAR_UNSEQ centralAtomIndex.cbegin(),
                      centralAtomIndex.cbegin() + centralAtomIndexSize.actDim,
                      [&](const Size_t& nAtom) mutable
                      {
                          std::fill(forcePreContract[nAtom].begin(),
                                    forcePreContract[nAtom].begin()
                                        + forcePreContractSize.actSize[nAtom],
                                    (Real)0);
                      });
    }
}

void Predictor::allocate_centralAtomIndex(const Size_t numberAtoms)
{
    bool resize = centralAtomIndexSize.checkResize(numberAtoms);
    if (resize)
    {
        centralAtomIndex.resize(numberAtoms);
        std::iota(centralAtomIndex.begin(), centralAtomIndex.end(), 0);
    }
}

void Predictor::compute_pairForcesDim(const DescriptorCollector& descriptorCollection)
{
    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        const Descriptor& descriptor = descriptorCollection.getDescriptor(key);
        std::shared_ptr<const NearestNeighborNSquare> neighborList =
            descriptor.get_neighborList_ptr();
        pairForcesDim[key].resize(neighborList->get_nAtoms());
        for (Size_t atom = 0; atom < pairForcesDim[key].size(); atom++)
        {
            // factor 3 for xyz
            pairForcesDim[key][atom] = 3 * neighborList->get_size(atom);
        }
    }
}

void Predictor::allocate_pairForces(const DescriptorCollector& descriptorCollection,
                                    const Size_t&              numberAtoms)
{
    compute_pairForcesDim(descriptorCollection);
    if (algoExecution == ExecutionPolicy::cpuSingleCore) allocate_pairForcesCPU(numberAtoms);
    else if (algoExecution == ExecutionPolicy::gpuStdLib) allocate_pairForcesGPU(numberAtoms);
}

void Predictor::allocate_pairForcesCPU(const Size_t& numberAtoms)
{
    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        pairForces[key]->resize(numberAtoms);
        for (std::size_t atom = 0; atom < pairForces[key]->size(); atom++)
        {
            // factor 3 for xyz
            (*pairForces[key])[atom].resize(pairForcesDim[key][atom]);
            std::fill((*pairForces[key])[atom].begin(), (*pairForces[key])[atom].end(), (Real)0);
        }
    }
}

void Predictor::allocate_pairForcesGPU(const Size_t& numberAtoms)
{

    bool resize;
    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        Vec2Real&        pairForces = *this->pairForces[key];
        ArrayResizing2D& pairForcesSize = this->pairForcesSize[key];
        resize = pairForcesSize.checkResize1Dim(numberAtoms);
        if (resize) { pairForcesSize.resizeArray1Dim(pairForces, pairForcesDim[key]); }
        else
        {
            resize = pairForcesSize.checkResize2Dim(pairForcesDim[key]);
            if (resize) { pairForcesSize.resizeArray2Dim(pairForces, pairForcesDim[key]); }
        }
        std::for_each(VASPML_PAR_UNSEQ centralAtomIndex.cbegin(),
                      centralAtomIndex.cbegin() + centralAtomIndexSize.actDim,
                      [&](const Size_t& indx)
                      {
                          std::fill(pairForces[indx].begin(),
                                    pairForces[indx].begin() + pairForcesSize.actSize[indx],
                                    (Real)0);
                      });
    }
}

void Predictor::allocate_centralForces(const Int numberIons)
{
    if (algoExecution == ExecutionPolicy::cpuSingleCore) allocate_centralForcesCPU(numberIons);
    else if (algoExecution == ExecutionPolicy::gpuStdLib) allocate_centralForcesGPU(numberIons);
}

void Predictor::allocate_centralForcesCPU(const Int numberIons)
{

    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        centralForces[key]->resize(3 * numberIons);
        std::fill(centralForces[key]->begin(), centralForces[key]->end(), (Real)0.0);
    }
}

void Predictor::allocate_centralForcesGPU(const Int numberIons)
{

    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        bool             resize = centralForcesSize[key].checkResize(3 * numberIons);
        Vec1Real&        centralForces = *this->centralForces[key];
        ArrayResizing1D& centralForcesSize = this->centralForcesSize[key];
        if (resize) { centralForces.resize(3 * numberIons); }
        std::fill(VASPML_PAR_UNSEQ centralForces.begin(),
                  centralForces.begin() + centralForcesSize.actDim,
                  (Real)0.0);
    }
}

void Predictor::allocate_stressTensor(void)
{
    if (algoExecution == ExecutionPolicy::cpuSingleCore) allocate_stressTensorCPU();
    else if (algoExecution == ExecutionPolicy::gpuStdLib) allocate_stressTensorGPU();
}

void Predictor::allocate_stressTensorCPU(void)
{
    for (const String& key : constants::descriptorKeyList)
    {
        stressTensor[key]->resize(9);
        std::fill(stressTensor[key]->begin(), stressTensor[key]->end(), (Real)0);
    }
    totalStressTensor->resize(9);
    std::fill(totalStressTensor->begin(), totalStressTensor->end(), (Real)0);
}

void Predictor::allocate_stressTensorGPU(void)
{

    bool resize = totalStressTensorSize.checkResize(9);
    if (resize)
    {
        totalStressTensor->resize(9);
        for (const String& key : constants::descriptorKeyList)
        {
            if (weights[key] <= 0) continue;
            stressTensor[key]->resize(9);
        }
    }
    // parallel to move/ keep the data on the device
    std::fill(VASPML_PAR_UNSEQ totalStressTensor->begin(), totalStressTensor->end(), (Real)0);
    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        Vec1Real& stressTensor = *this->stressTensor[key];
        std::fill(VASPML_PAR_UNSEQ stressTensor.begin(), stressTensor.end(), (Real)0);
    }
}

void Predictor::allocate_atomicForces(void)
{
    if (algoExecution == ExecutionPolicy::cpuSingleCore) allocate_atomicForcesCPU();
    else if (algoExecution == ExecutionPolicy::gpuStdLib) allocate_atomicForcesGPU();
}

void Predictor::allocate_atomicForcesCPU(void)
{

    // all descriptors have to have the same number of central atoms
    // choose one which is initialized
    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        atomicForces->resize(centralForces[key]->size());
        break;
    }
    std::fill(atomicForces->begin(), atomicForces->end(), (Real)0);
}

void Predictor::allocate_atomicForcesGPU(void)
{

    // all descriptors have to have the same number of central atoms
    // choose one which is initialized
    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        bool resize = atomicForcesSize.checkResize(centralForcesSize[key].actDim);
        if (resize) { atomicForces->resize(centralForcesSize[key].actDim); }
        break;
    }
    std::fill(VASPML_PAR_UNSEQ atomicForces->begin(),
              atomicForces->begin() + atomicForcesSize.actDim,
              (Real)0);
}

void Predictor::compute_tempForceVectorDim(const ShVec1Int& nAtomsType)
{
    tempForceVectorDim.resize(nAtomsType->size());
    for (std::size_t typeStruc = 0; typeStruc < nAtomsType->size(); typeStruc++)
        tempForceVectorDim[typeStruc] = (*nAtomsType)[typeStruc];
}

void Predictor::allocate_tempForceVector(const ShVec1Int& nAtomsType)
{
    compute_tempForceVectorDim(nAtomsType);
    if (algoExecution == ExecutionPolicy::cpuSingleCore) allocate_tempForceVectorCPU();
    else if (algoExecution == ExecutionPolicy::gpuStdLib) allocate_tempForceVectorGPU();
}

void Predictor::allocate_tempForceVectorCPU(void)
{
    tempForceVector.resize(tempForceVectorDim.size());
    for (Size_t typeStruc = 0; typeStruc < tempForceVectorDim.size(); typeStruc++)
    {
        tempForceVector[typeStruc].resize(tempForceVectorDim[typeStruc]);
    }
}

void Predictor::allocate_tempForceVectorGPU(void)
{
    bool resize = tempForceVectorSize.checkResize1Dim(tempForceVectorDim.size());
    if (resize) { tempForceVectorSize.resizeArray1Dim(tempForceVector, tempForceVectorDim); }
    else
    {
        resize = tempForceVectorSize.checkResize2Dim(tempForceVectorDim);
        if (resize) { tempForceVectorSize.resizeArray2Dim(tempForceVector, tempForceVectorDim); }
    }
}

void Predictor::update(const Kernel& kernel)
{
    typeMap = kernel.get_typeMap();
    updateEnergy(kernel);
    updateForces(kernel);
}

void Predictor::updateEnergy(const Kernel& kernel)
{

    const ShVec1Int& nAtomsType = kernel.get_nAtomsType();
    const Int        numberIons = vector_tools::sum(*nAtomsType);
    allocate_energyArray(nAtomsType->size());
    // TODO include the LTOTEN_SYSTEM tag for scaling
    Real scaleFactor = (Real)1.0 / (Real)numberIons;
    compute_referenceEnergyTotal(nAtomsType);
    totalEnergy = (Real)0;
    for (std::size_t typeStruc = 0; typeStruc < nAtomsType->size(); typeStruc++)
    {
        std::size_t typeForceField = typeMap->toType(typeStruc);
        compute_energyArray(kernel,
                            typeStruc,
                            typeForceField,
                            (*nAtomsType)[typeStruc],
                            scaleFactor);
        totalEnergy += ((Real)numberIons) * computeEnergyPerType(typeStruc, typeForceField);
    }
    finalizeEnergyComputation(nAtomsType);
}

void Predictor::updateForces(const Kernel& kernel)
{

    const ShVec1Int&           nAtomsType = kernel.get_nAtomsType();
    const Int                  numberIons = vector_tools::sum(*nAtomsType);
    const DescriptorCollector& descriptorCollection = kernel.get_descriptorCollection();
    Size_t                     numberAtoms = vector_tools::sum(*nAtomsType);

    allocate_derivativeMatrix(nAtomsType);
    allocate_forcePreContract(kernel.get_descriptorCollection(), numberAtoms);
    allocate_pairForces(descriptorCollection, numberAtoms);
    allocate_centralForces(numberIons);
    allocate_tempForceVector(nAtomsType);

    for (std::size_t typeStruc = 0; typeStruc < nAtomsType->size(); typeStruc++)
    {
        std::size_t typeForceField = typeMap->toType(typeStruc);
        compute_derivativeMatrix(kernel, tempForceVector[typeStruc], typeStruc, typeForceField);
    }

    compute_forcePreContract(descriptorCollection, *typeMap);
    computeForceArrays(descriptorCollection);
    forceArraysComputed = true;
}

void Predictor::compute_energyArray(const Kernel&     kernel,
                                    const std::size_t typeStruc,
                                    const std::size_t typeForceField,
                                    const std::size_t atomsPerType,
                                    const Real        scaleFactor)
{

    for (std::size_t atom = 0; atom < atomsPerType; atom++)
    {
        linalg::scaleVectorPlusVector(scaleFactor,
                                      kernel.get_kernelMatrixAtom(typeStruc, atom),
                                      (*energyArray)[typeStruc],
                                      (*numberLocalRefConfs)[typeForceField],
                                      linalgContext);
    }
}

Real Predictor::computeEnergyPerType(const std::size_t typeStruc, const std::size_t typeForceField)
{
    return linalg::dotProduct((*energyArray)[typeStruc],
                              regressionCoefficients->get_slice(typeForceField),
                              (*numberLocalRefConfs)[typeForceField],
                              linalgContext);
}

void Predictor::compute_referenceEnergyTotal(const ShVec1Int& atomsPerType)
{
    referenceEnergyTotal = (Real)0;
    for (std::size_t typeStruc = 0; typeStruc < atomsPerType->size(); typeStruc++)
    {
        std::size_t typeFF = typeMap->toType(typeStruc);
        //      for ( std::size_t atom = 0; atom < ( std::size_t )(*atomsPerType)[ typeStruc ];
        //      atom++ ){
        //         referenceEnergyTotal += (*referenceEnergyPerType)[ typeFF ];
        //      }
        referenceEnergyTotal += (*atomsPerType)[typeStruc] * (*referenceEnergyPerType)[typeFF];
    }
}

void Predictor::finalizeEnergyComputation(const ShVec1Int& atomsPerType)
{

    if (energyType == TotalEnergyType::IsolatedAtom) { totalEnergy += referenceEnergyTotal; }
    else if (energyType == TotalEnergyType::AverageTrainEnergy)
    {
        for (std::size_t typeStruc = 0; typeStruc < atomsPerType->size(); typeStruc++)
        {
            std::size_t typeFF = typeMap->toType(typeStruc);
            totalEnergy += (*atomsPerType)[typeStruc]
                         * ((*referenceEnergyPerType)[typeFF] + averageTrainEnergy);
        }
    }
}

void Predictor::compute_derivativeMatrix(const Kernel&     kernel,
                                         Vec1Real&         tempForceVector,
                                         const std::size_t typeStruc,
                                         const std::size_t typeForceField)
{

    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        kernel.computeDerivativeKernel((*derivativeMatrix[key])[typeStruc],
                                       tempForceVector,
                                       descriptorsRefConfsWeighted[key]->get_slice(typeForceField),
                                       regressionCoefficients->get_slice(typeForceField),
                                       typeStruc,
                                       typeForceField,
                                       key);
    }
}

void Predictor::compute_forcePreContract(const DescriptorCollector& descriptorCollection,
                                         const TypeMap&             typeMap)
{

    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        const Descriptor& descriptor = descriptorCollection.getDescriptor(key);
        descriptor.compute_forcePreContract(*derivativeMatrix[key],
                                            *forcePreContract[key],
                                            typeMap);
        //for ( int i = 0; i < (*forcePreContract[ key ]).size(); i++ )
        //{
        //    for ( int j = 0; j < (*forcePreContract[ key ])[i].size(); j++ )
        //    {
        //       std::cout << (*forcePreContract[ key ])[i][j] << std::endl;
        //    }
        //}
    }
}

void Predictor::computeForceArrays(const DescriptorCollector& descriptorCollection)
{

    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        const Descriptor& descriptor = descriptorCollection.getDescriptor(key);
        descriptor.computeForceTerms(*forcePreContract[key], *pairForces[key], *centralForces[key]);
    }
}

void Predictor::compute_atomicForces(const DescriptorCollector& descriptorCollection)
{

    if (forceArraysComputed)
    {
        allocate_atomicForces();
        switch (algoExecution)
        {
        case ExecutionPolicy::cpuSingleCore:
            compute_atomicForcesCPU(descriptorCollection);
            break;
        case ExecutionPolicy::gpuStdLib:
            VASPML_PARALLEL(
                compute_atomicForcesGPU( descriptorCollection );
            );
            break;
        }
    }
    else
    {
        throw std::runtime_error(
            "ERROR: void Predictor::compute_atomicForces( const DescriptorCollector& "
            "descriptorCollection )\n"
            "force arrays were not computed. Please call void Predictor::updateForcesAndStress( "
            "const Kernel& kernel ) before.");
    }
}

void Predictor::compute_atomicForcesCPU(const DescriptorCollector& descriptorCollection)
{
    Vec1Real& atomicForces = (*this->atomicForces);
    for (const String& key : constants::descriptorKeyList)
    {
        const NearestNeighborNSquare neighborList =
            *descriptorCollection.getDescriptor(key).get_neighborList_ptr();
        Vec2Real&      pairForces = (*this->pairForces[key]);
        Vec1Real&      centralForces = (*this->centralForces[key]);
        const Vec2Int& neighborIndex = neighborList.get_globalIndex();
        for (std::size_t centralAtom = 0; centralAtom < neighborList.get_nAtoms(); centralAtom++)
        {
            atomicForces[3 * centralAtom] += centralForces[3 * centralAtom];
            atomicForces[3 * centralAtom + 1] += centralForces[3 * centralAtom + 1];
            atomicForces[3 * centralAtom + 2] += centralForces[3 * centralAtom + 2];
            SumFunctorForce functor(atomicForces, pairForces[centralAtom]);
            std::for_each(neighborIndex[centralAtom].cbegin(),
                          neighborIndex[centralAtom].cend(),
                          functor);
        }
    }
}

void Predictor::compute_atomicForcesGPU(const DescriptorCollector& descriptorCollection)
{

    Vec1Real&                atomicForces = (*this->atomicForces);
    std::size_t              nAtoms = atomicForces.size();
    std::vector<std::size_t> centralAtomIndex(nAtoms);
    std::iota(centralAtomIndex.begin(), centralAtomIndex.end(), 0);
    for (const String& key : constants::descriptorKeyList)
    {
        const NearestNeighborNSquare neighborList =
            *descriptorCollection.getDescriptor(key).get_neighborList_ptr();
        Vec2Real&      pairForces = (*this->pairForces[key]);
        Vec1Real&      centralForces = (*this->centralForces[key]);
        const Vec2Int& neighborIndex = neighborList.get_globalIndex();
        std::for_each(VASPML_PAR_UNSEQ centralAtomIndex.cbegin(),
                      centralAtomIndex.cbegin() + centralAtomIndexSize.actDim,
                      [&](const std::size_t& centralAtom) mutable
                      {
                          atomicForces[3 * centralAtom] += centralForces[3 * centralAtom];
                          atomicForces[3 * centralAtom + 1] += centralForces[3 * centralAtom + 1];
                          atomicForces[3 * centralAtom + 2] += centralForces[3 * centralAtom + 2];

                          //SumFunctorForce functor( atomicForces, pairForces[ centralAtom ] );
                          //std::for_each( neighborIndex[ centralAtom ].cbegin(), neighborIndex[
                          //centralAtom ].cend(), functor );
                          std::size_t counter = 0;
                          std::for_each(neighborIndex[centralAtom].cbegin(),
                                        neighborIndex[centralAtom].cend(),
                                        [&](const Int& neighborIndex) mutable
                                        {
                                            atomicForces[3 * neighborIndex] +=
                                                pairForces[centralAtom][3 * counter];
                                            atomicForces[3 * neighborIndex + 1] +=
                                                pairForces[centralAtom][3 * counter + 1];
                                            atomicForces[3 * neighborIndex + 2] +=
                                                pairForces[centralAtom][3 * counter + 2];
                                            counter++;
                                        });
                      });
    }
}

void Predictor::compute_stressTensor(const DescriptorCollector& descriptorCollection,
                                     const Real&                volume)
{
    if (forceArraysComputed)
    {
        allocate_stressTensor();
        switch (algoExecution)
        {
        case ExecutionPolicy::cpuSingleCore:
            compute_stressTensorCPU(descriptorCollection, volume);
            break;
        case ExecutionPolicy::gpuStdLib:
            VASPML_PARALLEL(
                 compute_stressTensorGPU( descriptorCollection, volume );
             );
            break;
        }
    }
    else
    {
        global_scope::tutor.bug("ERROR: void Predictor::compute_atomicForces( const "
                                "DescriptorCollector& descriptorCollection )\n"
                                "force arrays were not computed. Please call void "
                                "Predictor::updateForcesAndStress( const Kernel& kernel ) before.");
    }
}

void Predictor::compute_stressTensorCPU(const DescriptorCollector& descriptorCollection,
                                        const Real&                volume)
{

    Vec1Real&   totalStressTensor = (*this->totalStressTensor);
    const Real& inverseVolume = (Real)-1.0 / volume;
    for (const String& key : constants::descriptorKeyList)
    {
        const NearestNeighborNSquare& neighborList =
            *descriptorCollection.getDescriptor(key).get_neighborList_ptr();
        Vec1Real&         stressTensor = (*this->stressTensor[key]);
        Vec2Real&         pairForces = (*this->pairForces[key]);
        const Vec2Int&    neighborIndex = neighborList.get_globalIndex();
        const Vec2Real    connectionVector = neighborList.get_connectionVector();
        const std::size_t nAtoms = neighborList.get_nAtoms();

        for (std::size_t centralAtom = 0; centralAtom < nAtoms; centralAtom++)
        {
            SumFunctorStress functor(stressTensor,
                                     pairForces[centralAtom],
                                     connectionVector[centralAtom],
                                     inverseVolume);
            std::for_each(neighborIndex[centralAtom].cbegin(),
                          neighborIndex[centralAtom].cend(),
                          functor);
        }
        // symmetrize
        stressTensor[3] = stressTensor[1];
        stressTensor[6] = stressTensor[2];
        stressTensor[7] = stressTensor[5];

        totalStressTensor[0] += stressTensor[0];
        totalStressTensor[1] += stressTensor[1];
        totalStressTensor[2] += stressTensor[2];
        totalStressTensor[4] += stressTensor[4];
        totalStressTensor[5] += stressTensor[5];
        totalStressTensor[8] += stressTensor[8];
    }
    // symmetrize total stress
    totalStressTensor[3] = totalStressTensor[1];
    totalStressTensor[6] = totalStressTensor[2];
    totalStressTensor[7] = totalStressTensor[5];
}

void Predictor::compute_stressTensorGPU(const DescriptorCollector& descriptorCollection,
                                        const Real&                volume)
{

    Vec1Real&   totalStressTensor = (*this->totalStressTensor);
    const Real& inverseVolume = (Real)-1.0 / volume;
    for (const String& key : constants::descriptorKeyList)
    {
        const NearestNeighborNSquare& neighborList =
            *descriptorCollection.getDescriptor(key).get_neighborList_ptr();
        Vec1Real&                stressTensor = (*this->stressTensor[key]);
        Vec2Real&                pairForces = (*this->pairForces[key]);
        const Vec2Int&           neighborIndex = neighborList.get_globalIndex();
        const Vec2Real           connectionVector = neighborList.get_connectionVector();
        const std::size_t        nAtoms = neighborList.get_nAtoms();
        std::vector<std::size_t> centralAtomIndex(nAtoms);
        std::iota(centralAtomIndex.begin(), centralAtomIndex.end(), 0);
        std::for_each( // seq,
            centralAtomIndex.cbegin(),
            centralAtomIndex.cbegin() + nAtoms,
            [&](const std::size_t& centralAtom) mutable
            {
                // functor can not be used on GPU because of the reference passing to contsructor
                // and copy passing in std::for_each
                //SumFunctorStress functor( stressTensor, pairForces[ centralAtom ],
                //         connectionVector[ centralAtom ], inverseVolume );
                Int counter = 0;
                std::for_each(
                    neighborIndex[centralAtom].cbegin(),
                    neighborIndex[centralAtom].cend(),
                    // NOTE (aS): Commented out unused parameter, is this needed here for
                    // some reason? Can this be expressed differently?
                    [&](const Int& /* neighborIndex */) mutable
                    {
                        // compute xx term
                        stressTensor[0] -= inverseVolume * pairForces[centralAtom][3 * counter]
                                         * connectionVector[centralAtom][3 * counter];
                        // compute xy term
                        stressTensor[1] -= inverseVolume * pairForces[centralAtom][3 * counter]
                                         * connectionVector[centralAtom][3 * counter + 1];
                        // compute xz term
                        stressTensor[2] -= inverseVolume * pairForces[centralAtom][3 * counter]
                                         * connectionVector[centralAtom][3 * counter + 2];
                        // compute yy term
                        stressTensor[4] -= inverseVolume * pairForces[centralAtom][3 * counter + 1]
                                         * connectionVector[centralAtom][3 * counter + 1];
                        // compute yz term
                        stressTensor[5] -= inverseVolume * pairForces[centralAtom][3 * counter + 1]
                                         * connectionVector[centralAtom][3 * counter + 2];
                        // compute zz term
                        stressTensor[8] -= inverseVolume * pairForces[centralAtom][3 * counter + 2]
                                         * connectionVector[centralAtom][3 * counter + 2];
                        counter++;
                    });
            });
        // symmetrize
        stressTensor[3] = stressTensor[1];
        stressTensor[6] = stressTensor[2];
        stressTensor[7] = stressTensor[5];

        totalStressTensor[0] += stressTensor[0];
        totalStressTensor[1] += stressTensor[1];
        totalStressTensor[2] += stressTensor[2];
        totalStressTensor[4] += stressTensor[4];
        totalStressTensor[5] += stressTensor[5];
        totalStressTensor[8] += stressTensor[8];
    }
    // symmetrize total stress
    totalStressTensor[3] = totalStressTensor[1];
    totalStressTensor[6] = totalStressTensor[2];
    totalStressTensor[7] = totalStressTensor[5];
}

const Real& Predictor::get_pairForcesX(const String&     key,
                                       const std::size_t centralAtom,
                                       const std::size_t neighborAtom) const
{
    VASPML_DEBUG_L1(
        if (centralAtom >= pairForces.at(key)->size())
        {
            throw std::runtime_error(
                "ERROR: const Real& Predictor::get_pairForcesX( const String& key, const "
                "std::size_t centralAtom, const std::size_t neighborAtom )const\n"
                "Index centralAtom is larger than number of atoms\n");
        }
        if (neighborAtom >= (*pairForces.at(key))[centralAtom].size() / 3)
        {
            throw std::runtime_error(
                "ERROR: const Real& Predictor::get_pairForcesX( const String& key, const "
                "std::size_t centralAtom, const std::size_t neighborAtom )const\n"
                "Index neighborAtom is larger than number of neighbors\n");
        }
    );
    return (*pairForces.at(key))[centralAtom][3 * neighborAtom];
}

const Real& Predictor::get_pairForcesY(const String&     key,
                                       const std::size_t centralAtom,
                                       const std::size_t neighborAtom) const
{
    VASPML_DEBUG_L1(
        if (centralAtom >= pairForces.at(key)->size())
        {
            throw std::runtime_error(
                "ERROR: const Real& Predictor::get_pairForcesY( const String& key, const "
                "std::size_t centralAtom, const std::size_t neighborAtom )const\n"
                "Index centralAtom is larger than number of atoms\n");
        }
        if (neighborAtom >= (*pairForces.at(key))[centralAtom].size() / 3)
        {
            throw std::runtime_error(
                "ERROR: const Real& Predictor::get_pairForcesY( const String& key, const "
                "std::size_t centralAtom, const std::size_t neighborAtom )const\n"
                "Index neighborAtom is larger than number of neighbors\n");
        }
    );
    return (*pairForces.at(key))[centralAtom][3 * neighborAtom + 1];
}

const Real& Predictor::get_pairForcesZ(const String&     key,
                                       const std::size_t centralAtom,
                                       const std::size_t neighborAtom) const
{
    VASPML_DEBUG_L1(
        if (centralAtom >= pairForces.at(key)->size())
        {
            throw std::runtime_error(
                "ERROR: const Real& Predictor::get_pairForcesZ( const String& key, const "
                "std::size_t centralAtom, const std::size_t neighborAtom )const\n"
                "Index centralAtom is larger than number of atoms\n");
        }
        if (neighborAtom >= (*pairForces.at(key))[centralAtom].size() / 3)
        {
            throw std::runtime_error(
                "ERROR: const Real& Predictor::get_pairForcesZ( const String& key, const "
                "std::size_t centralAtom, const std::size_t neighborAtom )const\n"
                "Index neighborAtom is larger than number of neighbors\n");
        }
    );
    return (*pairForces.at(key))[centralAtom][3 * neighborAtom + 2];
}

const std::tuple<const Real&, const Real&, const Real&> Predictor::get_pairForces(
    const String&     key,
    const std::size_t centralAtom,
    const std::size_t neighborAtom) const
{
    VASPML_DEBUG_L1(
        if (centralAtom >= pairForces.at(key)->size())
        {
            throw std::runtime_error(
                "ERROR: const Real& Predictor::get_pairForces( const String& key, const "
                "std::size_t centralAtom, const std::size_t neighborAtom )const\n"
                "Index centralAtom is larger than number of atoms\n");
        }
        if (neighborAtom >= (*pairForces.at(key))[centralAtom].size() / 3)
        {
            throw std::runtime_error(
                "ERROR: const Real& Predictor::get_pairForces( const String& key, const "
                "std::size_t centralAtom, const std::size_t neighborAtom )const\n"
                "Index neighborAtom is larger than number of neighbors\n");
        }
    );
    return std::tie(get_pairForcesX(key, centralAtom, neighborAtom),
                    get_pairForcesY(key, centralAtom, neighborAtom),
                    get_pairForcesZ(key, centralAtom, neighborAtom));
}

const Real& Predictor::get_centralForcesX(const String& key, const std::size_t centralAtom) const
{
    VASPML_DEBUG_L1(
        if (centralAtom >= centralForces.at(key)->size() / 3)
        {
            throw std::runtime_error("ERROR: const Real& Predictor::get_centralForcesX( const "
                                     "String& key, const std::size_t centralAtom )const\n"
                                     "Index centralAtom is larger than number of atoms\n");
        }
    );
    return (*centralForces.at(key))[3 * centralAtom];
}

const Real& Predictor::get_centralForcesY(const String& key, const std::size_t centralAtom) const
{
    VASPML_DEBUG_L1(
        if (centralAtom >= centralForces.at(key)->size() / 3)
        {
            throw std::runtime_error("ERROR: const Real& Predictor::get_centralForcesY( const "
                                     "String& key, const std::size_t centralAtom )const\n"
                                     "Index centralAtom is larger than number of atoms\n");
        }
    );
    return (*centralForces.at(key))[3 * centralAtom + 1];
}

const Real& Predictor::get_centralForcesZ(const String& key, const std::size_t centralAtom) const
{
    VASPML_DEBUG_L1(
        if (centralAtom >= centralForces.at(key)->size() / 3)
        {
            throw std::runtime_error("ERROR: const Real& Predictor::get_centralForcesZ( const "
                                     "String& key, const std::size_t centralAtom )const\n"
                                     "Index centralAtom is larger than number of atoms\n");
        }
    );
    return (*centralForces.at(key))[3 * centralAtom + 2];
}

const std::tuple<const Real&, const Real&, const Real&> Predictor::get_centralForces(
    const String&     key,
    const std::size_t centralAtom) const
{
    VASPML_DEBUG_L1(
        if (centralAtom >= centralForces.at(key)->size() / 3)
        {
            throw std::runtime_error("ERROR: const Real& Predictor::get_centralForces( const "
                                     "String& key, const std::size_t centralAtom )const\n"
                                     "Index centralAtom is larger than number of atoms\n");
        }
    );
    return std::tie(get_centralForcesX(key, centralAtom),
                    get_centralForcesY(key, centralAtom),
                    get_centralForcesZ(key, centralAtom));
}

const Real& Predictor::get_atomicForcesX(const std::size_t centralAtom) const
{
    VASPML_DEBUG_L1(
        if (centralAtom >= atomicForces->size() / 3)
        {
            throw std::runtime_error("ERROR: const Real& Predictor::get_atomicForcesX( const "
                                     "std::size_t centralAtom )const\n"
                                     "Index centralAtom is larger than number of atoms\n");
        }
    );
    return (*atomicForces)[3 * centralAtom];
}

const Real& Predictor::get_atomicForcesY(const std::size_t centralAtom) const
{
    VASPML_DEBUG_L1(
        if (centralAtom >= atomicForces->size() / 3)
        {
            throw std::runtime_error("ERROR: const Real& Predictor::get_atomicForcesY( const "
                                     "std::size_t centralAtom )const\n"
                                     "Index centralAtom is larger than number of atoms\n");
        }
    );
    return (*atomicForces)[3 * centralAtom + 1];
}

const Real& Predictor::get_atomicForcesZ(const std::size_t centralAtom) const
{
    VASPML_DEBUG_L1(
        if (centralAtom >= atomicForces->size() / 3)
        {
            throw std::runtime_error("ERROR: const Real& Predictor::get_atomicForcesZ( const "
                                     "std::size_t centralAtom )const\n"
                                     "Index centralAtom is larger than number of atoms\n");
        }
    );
    return (*atomicForces)[3 * centralAtom + 2];
}

const std::tuple<const Real&, const Real&, const Real&> Predictor::get_atomicForces(
    const std::size_t centralAtom) const
{

    VASPML_DEBUG_L1(
        if (centralAtom >= atomicForces->size() / 3)
        {
            throw std::runtime_error("ERROR: const Real& Predictor::get_atomicForces( const "
                                     "std::size_t centralAtom )const\n"
                                     "Index centralAtom is larger than number of atoms\n");
        }
    );
    return std::tie(get_atomicForcesX(centralAtom),
                    get_atomicForcesY(centralAtom),
                    get_atomicForcesZ(centralAtom));
}

const Real& Predictor::get_stressTensor(const String&     key,
                                        const std::size_t indx0,
                                        const std::size_t indx1) const
{

    VASPML_DEBUG_L1(
        if (indx0 >= 3)
        {
            throw std::runtime_error(
                "ERROR: Predictor::get_stressTensor( const String& key, const std::size_t indx0, "
                "const std::size_t indx1 )const\n"
                "Index indx0 is larger than 3 ( only 3 space dimensions available )\n");
        }
        if (indx1 >= 3)
        {
            throw std::runtime_error(
                "ERROR: Predictor::get_stressTensor( const String& key, const std::size_t indx0, "
                "const std::size_t indx1 )const\n"
                "Index indx1 is larger than 3 ( only 3 space dimensions available )\n");
        }
    );
    return (*stressTensor.at(key))[indx0 * 3 + indx1];
}

const Real& Predictor::get_totalStressTensor(const std::size_t indx0, const std::size_t indx1) const
{
    VASPML_DEBUG_L1(
        if (indx0 >= 3)
        {
            throw std::runtime_error(
                "ERROR: Predictor::get_totalStressTensor( const std::size_t indx0, const "
                "std::size_t indx1 )const\n"
                "Index indx0 is larger than 3 ( only 3 space dimensions available )\n");
        }
        if (indx1 >= 3)
        {
            throw std::runtime_error(
                "ERROR: Predictor::get_totalStressTensor( const std::size_t indx0, const "
                "std::size_t indx1 )const\n"
                "Index indx1 is larger than 3 ( only 3 space dimensions available )\n");
        }
    );
    return (*totalStressTensor)[indx0 * 3 + indx1];
}

const Real& Predictor::get_totalEnergy(void) const
{
    return totalEnergy;
}

void Predictor::write_atomicForceToScreen(void) const
{
    for (std::size_t atom = 0; atom < atomicForces->size(); atom += 3)
    {
        std::cout << str("%-24.16E %-24.16E %-24.16E\n",
                         (*atomicForces)[atom] * constants::FUNIT,
                         (*atomicForces)[atom + 1] * constants::FUNIT,
                         (*atomicForces)[atom + 2] * constants::FUNIT);
    }

    return;
}
