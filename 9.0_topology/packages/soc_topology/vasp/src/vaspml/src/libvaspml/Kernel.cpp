#include "Kernel.hpp"
#include "TagTranslator.hpp"
#include "constants.hpp"
#include "debug.hpp"
#include "math.hpp"

#include <algorithm>
#include <string>

#include "utils.hpp"

using namespace vaspml;

namespace vaspml
{

std::map<String, Real> generate_weightsMap(const IoHandlerML_FF& inputParameters)
{

    std::map<String, Real> weights;
    for (const String& key : constants::descriptorKeyList)
    {
        std::string tag = key + "-weight";
        weights[key] = inputParameters[tag].cget<Real>();
    }
    return weights;
}

std::map<String, ShVec1Int> generate_featureSpaceMap(const IoHandlerML_FF& inputParameters)
{
    std::map<String, ShVec1Int> featureSpaceSize;
    const Vec1Int& locRef = inputParameters["number-local-reference-configs"].dcget<ShVec1Int>();
    for (const String& key : constants::descriptorKeyList)
    {
        std::string tag = key + "-number-descriptors-per-type";
        if (std::get_if<Int>(&inputParameters[tag]))
        {
            featureSpaceSize[key] = std::make_shared<Vec1Int>(
                math::vectorTimesScalar(locRef, inputParameters[tag].cget<Int>()));
        }
        else if (std::get_if<ShVec1Int>(&inputParameters[tag]))
        {
            featureSpaceSize[key] = std::make_shared<Vec1Int>(
                math::elementwiseProduct(locRef, inputParameters[tag].dcget<ShVec1Int>()));
        }
    }
    return featureSpaceSize;
}

std::map<String, std::shared_ptr<ShmemArray2DVariableLen<Real>>> generate_ShmemArrayMap(
    const std::map<String, Real>&      weights,
    const std::map<String, ShVec1Int>& featureSpaceSize,
    const std::shared_ptr<MlMPI>&      mpiIn)
{
    std::map<String, std::shared_ptr<ShmemArray2DVariableLen<Real>>> arrays;
    for (auto& [key, item] : featureSpaceSize)
    {
        arrays[key] = std::make_shared<ShmemArray2DVariableLen<Real>>(
            weights.at(key) > 0 ? *featureSpaceSize.at(key) : Vec1Int{0},
            mpiIn);
    }
    return arrays;
}

} //namespace vaspml

Kernel::Kernel(const IoHandlerML_FF&         inputParameters,
               const std::shared_ptr<MlMPI>& mpiIn,
               const ShVec2Real&             kernelMatrix,
               const ExecutionPolicy         algoExecution) :
    weights(inputParameters.generate_weightsMap()),
    numberLocalRefConfs(std::make_shared<Vec1Int>(
        inputParameters["number-local-reference-configs"].dcget<ShVec1Int>())),
    numberDescriptors(inputParameters.get_numberDescriptorsMap()),
    featureSpaceSize(inputParameters.generate_featureSpaceMap()),
    descriptorsRefConfs(generate_ShmemArrayMap(weights, featureSpaceSize, mpiIn)),
    descriptorCollection(std::make_shared<DescriptorCollector>(0,
                                                               std::map<String, ShVec2Real>(),
                                                               nullptr,
                                                               algoExecution)),
    algoExecution(algoExecution)
{
    if (kernelMatrix == nullptr) { this->kernelMatrix = std::make_shared<Vec2Real>(); }
    else { this->kernelMatrix = kernelMatrix; }
    // filling descriptors into shared memory arrays
    fillDescriptor(inputParameters);
}

void Kernel::init(const IoHandlerML_FF&         inputParameters,
                  const std::shared_ptr<MlMPI>& mpiIn,
                  const ShVec2Real&             kernelMatrix,
                  const ExecutionPolicy         algoExecution)
{
    weights = inputParameters.generate_weightsMap();
    numberLocalRefConfs = std::make_shared<Vec1Int>(
        inputParameters["number-local-reference-configs"].dcget<ShVec1Int>());
    numberDescriptors = inputParameters.get_numberDescriptorsMap();
    featureSpaceSize = inputParameters.generate_featureSpaceMap();
    descriptorsRefConfs = generate_ShmemArrayMap(weights, featureSpaceSize, mpiIn);
    descriptorCollection = std::make_shared<DescriptorCollector>(0,
                                                                 std::map<String, ShVec2Real>(),
                                                                 nullptr,
                                                                 algoExecution);
    this->algoExecution = algoExecution;
    if (kernelMatrix == nullptr) { this->kernelMatrix = std::make_shared<Vec2Real>(); }
    else { this->kernelMatrix = kernelMatrix; }
    fillDescriptor(inputParameters);
}

void Kernel::fillDescriptor(const IoHandlerML_FF& inputParameters)
{
    // filling descriptors into shared memory arrays
    for (const String& key : constants::descriptorKeyList)
    {
        String tag = key + "-reference-configs";
        if (weights[key] > 0)
            descriptorsRefConfs[key]->set_value(inputParameters[tag].dcget<ShVec2Real>());
    }
}

void Kernel::updateKernel(const Frame& frame)
{
    const std::map<String, std::shared_ptr<Descriptor>>& structureDescriptors =
        frame.get_descriptors();
    VASPML_DEBUG_L1(
        for (std::size_t i = 0; i < constants::descriptorKeyList.size() - 1; i++)
        {
            for (std::size_t j = i + 1; j < constants::descriptorKeyList.size(); j++)
            {
                if (structureDescriptors.at(constants::descriptorKeyList[i])->get_nAtoms()
                    != structureDescriptors.at(constants::descriptorKeyList[j])->get_nAtoms())
                {
                    throw std::runtime_error(
                        "ERROR: void Kernel::update( const Frame& frame ) \n"
                        "the number of atoms does not agree in your descriptors");
                }
            }
        }
    );

    typeMap = frame.get_typeMap();
    nAtomsPerType = std::make_shared<Vec1Int>(
        structureDescriptors.at(constants::descriptorKeyList[0])->get_nAtomsType());
    allocate_kernelMatrix(*nAtomsPerType);
    setUpDescriptorCollector(structureDescriptors);
    compute_kernelMatrix();
}

void Kernel::compute_kernelMatrixDim(const Vec1Int& nAtomPerType)
{
    kernelMatrixDim.resize(typeMap->countStructureTypes());
    for (Size_t strucType = 0; strucType < typeMap->countStructureTypes(); strucType++)
    {
        Int forceFieldType = typeMap->toType(strucType);
        kernelMatrixDim[strucType] =
            nAtomPerType[strucType] * (*numberLocalRefConfs)[forceFieldType];
    }
}

void Kernel::allocate_kernelMatrix(const Vec1Int& nAtomPerType)
{
    compute_kernelMatrixDim(nAtomPerType);
    if (algoExecution == ExecutionPolicy::cpuSingleCore)
    {
        kernelMatrix->resize(typeMap->countStructureTypes());
        for (std::size_t strucType = 0; strucType < kernelMatrix->size(); strucType++)
        {
            (*kernelMatrix)[strucType].resize(kernelMatrixDim[strucType], (Real)0);
            std::fill((*kernelMatrix)[strucType].begin(),
                      (*kernelMatrix)[strucType].end(),
                      (Real)0);
        }
    }
    else if (algoExecution == ExecutionPolicy::gpuStdLib) { resizeKernelMatrixGPU(nAtomPerType); }
}

void Kernel::resizeKernelMatrixGPU(const Vec1Int& nAtomPerType)
{

    bool resize = kernelMatrixSize.checkResize1Dim(nAtomPerType.size());
    if (resize) { kernelMatrixSize.resizeArray1Dim(*kernelMatrix, kernelMatrixDim); }
    else
    {
        resize = kernelMatrixSize.checkResize2Dim(kernelMatrixDim);
        if (resize) kernelMatrixSize.resizeArray2Dim(*kernelMatrix, kernelMatrixDim);
    }
    // fill with zeros
    std::for_each(VASPML_PAR_UNSEQ
                  kernelMatrix->begin(),
                  kernelMatrix->begin() + kernelMatrixSize.act1Dim,
                  [](Vec1Real& slice) { std::fill(slice.begin(), slice.end(), (Real)0.0); });
}

void Kernel::setUpDescriptorCollector(
    const std::map<String, std::shared_ptr<Descriptor>>& structureDescriptors)
{

    // TODO
    descriptorCollection->setDescriptorMap(structureDescriptors);
    descriptorCollection->updateCollector();
    descriptorCollection->rearrangeSHS2Body(*typeMap);
    // descriptorCollection->writeDescriptorCollector();
}

void Kernel::compute_kernelMatrix(void)
{

    // unpack variables
    const std::map<String, ShVec2Real>& descriptorStructure =
        descriptorCollection->get_descriptorsNormalized();
    Vec1Int&  nAtomsPerType = (*this->nAtomsPerType);
    Vec1Int&  numberLocalRefConfs = (*this->numberLocalRefConfs);
    Vec2Real& kernelMatrix = (*this->kernelMatrix);

    for (const String& key : constants::descriptorKeyList)
    {
        if (weights[key] <= 0) continue;
        Vec1Int& numberDescriptors = (*this->numberDescriptors[key]);
        // extract descriptor of structure
        const Vec2Real& descriptor = (*descriptorStructure.at(key));
        for (std::size_t strucType = 0; strucType < nAtomsPerType.size(); strucType++)
        {
            Int            forceFieldType = typeMap->toType(strucType);
            constexpr Real one = (Real)1;
            linalg::matMul(linalg::NoTrans,
                           linalg::Trans,
                           nAtomsPerType[strucType],            // m
                           numberLocalRefConfs[forceFieldType], // n
                           numberDescriptors[forceFieldType],   // k
                           one,
                           descriptor[strucType],             // A\in\mathbb{R}^{Natom/times/Ndesc}
                           numberDescriptors[forceFieldType], // k
                           descriptorsRefConfs[key]->get_slice(
                               forceFieldType), //B^{T}\in\mathbb{R}^{Nref/times/Ndesc}
                           numberDescriptors[forceFieldType], // k
                           one,
                           kernelMatrix[strucType],
                           numberLocalRefConfs[forceFieldType],
                           linalgContext); // n
        }
    }
}

const std::map<String, std::shared_ptr<ShmemArray2DVariableLen<Real>>>
Kernel::get_descriptorsRefConfs(void) const
{
    return descriptorsRefConfs;
}

ShVec1Int Kernel::get_numberLocalRefConfs(void) const
{
    return numberLocalRefConfs;
}

const std::map<String, ShVec1Int>& Kernel::get_numberDescriptors(void) const
{
    return numberDescriptors;
}

const std::map<String, ShVec1Int>& Kernel::get_featureSpaceSize(void) const
{
    return featureSpaceSize;
}

const std::map<String, Real>& Kernel::get_weights(void) const
{
    return weights;
}

const ShVec2Real& Kernel::get_kernelMatrix(void) const
{
    return kernelMatrix;
}

const ShVec1Int& Kernel::get_nAtomsType(void) const
{
    return nAtomsPerType;
}

std::shared_ptr<const TypeMap> Kernel::get_typeMap(void) const
{
    return typeMap;
}

const Real* Kernel::get_kernelMatrixAtom(std::size_t type, std::size_t atomInType) const
{
    VASPML_DEBUG_L1(
        if (type >= nAtomsPerType->size())
        {
            throw std::runtime_error(
                "ERROR: const Real* Kernel::get_kernelMatrixAtom( std::size_t type,"
                "std::size_t atomInType )\n       type index "
                + std::to_string(type) + " out of bounds\n");
        }
    );
    VASPML_DEBUG_L1(
        if (atomInType >= (std::size_t)(*nAtomsPerType)[type])
        {
            throw std::runtime_error(
                "ERROR: const Real* Kernel::get_kernelMatrixAtom( std::size_t type,"
                "std::size_t atomInType )\n       atomInType index "
                + std::to_string(atomInType) + " out of bounds\n");
        }
    );

    Real*       ptr;
    std::size_t index = atomInType * (*numberLocalRefConfs)[typeMap->toType(type)];
    ptr = &(*kernelMatrix)[type][index];
    return ptr;
}

const DescriptorCollector& Kernel::get_descriptorCollection(void) const
{
    return *descriptorCollection;
}
