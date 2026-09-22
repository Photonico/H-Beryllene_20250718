#include "KernelPolynomial.hpp"
#include "Timer.hpp"
#include "constants.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <tuple>

using namespace vaspml;

KernelPolynomial::KernelPolynomial(const IoHandlerML_FF&         inputParameters,
                                   const std::shared_ptr<MlMPI>& mpiIn,
                                   const ShVec2Real&             polynomialDerivative,
                                   const ShVec2Real&             polynomialKernelNorm,
                                   const ExecutionPolicy         algoExecution) :
    Kernel(inputParameters, mpiIn, nullptr, algoExecution)
{

    kernelPower = inputParameters["SHS2-2-body-exponent"].cget<Int>();
    if (polynomialDerivative == nullptr)
    {
        this->polynomialDerivative = std::make_shared<Vec2Real>();
    }
    else { this->polynomialDerivative = polynomialDerivative; }
    if (polynomialKernelNorm == nullptr)
    {
        this->polynomialKernelNorm = std::make_shared<Vec2Real>();
    }
    else { this->polynomialKernelNorm = polynomialKernelNorm; }
}

void KernelPolynomial::init(const IoHandlerML_FF&         inputParameters,
                            const std::shared_ptr<MlMPI>& mpiIn,
                            const ShVec2Real&             polynomialDerivative,
                            const ShVec2Real&             polynomialKernelNorm,
                            const ExecutionPolicy         algoExecution)
{
    Kernel::init(inputParameters, mpiIn, nullptr, algoExecution);
    kernelPower = inputParameters["SHS2-2-body-exponent"].cget<Int>();
    this->algoExecution = algoExecution;
    if (polynomialDerivative == nullptr)
    {
        this->polynomialDerivative = std::make_shared<Vec2Real>();
    }
    else { this->polynomialDerivative = polynomialDerivative; }
    if (polynomialKernelNorm == nullptr)
    {
        this->polynomialKernelNorm = std::make_shared<Vec2Real>();
    }
    else { this->polynomialKernelNorm = polynomialKernelNorm; }
}

void KernelPolynomial::allocateDerivativeArrays(void)
{

    if (algoExecution == ExecutionPolicy::cpuSingleCore)
    {
        polynomialDerivative->resize(this->kernelMatrix->size());
        polynomialKernelNorm->resize(this->kernelMatrix->size());
        for (std::size_t type = 0; type < polynomialDerivative->size(); type++)
        {
            (*polynomialDerivative)[type].resize((*kernelMatrix)[type].size());
            (*polynomialKernelNorm)[type].resize((*kernelMatrix)[type].size());
        }
        polynomialKernel = kernelMatrix;
    }
    else if (algoExecution == ExecutionPolicy::gpuStdLib) { allocateDerivativeArraysGPU(); }
}

void KernelPolynomial::allocateDerivativeArraysGPU(void)
{

    const Size_t numberAtoms =
        descriptorCollection->getDescriptor("SHS2-2-body").get_neighborList().get_numberAtoms();

    bool resize = centralAtomIndexSize.checkResize(numberAtoms);
    // resize central atom index
    if (resize)
    {
        centralAtomIndex.resize(numberAtoms);
        std::iota(centralAtomIndex.begin(), centralAtomIndex.end(), 0);
    }

    resize = polynomialDerivativeSize.checkResize1Dim(nAtomsPerType->size());
    if (resize)
    {
        polynomialDerivativeSize.resizeArray1Dim(*polynomialDerivative, kernelMatrixDim);
        polynomialDerivativeSize.resizeArray1Dim(*polynomialKernelNorm, kernelMatrixDim);
        polynomialKernel = kernelMatrix;
    }
    else
    {
        resize = polynomialDerivativeSize.checkResize2Dim(kernelMatrixDim);
        if (resize)
        {
            polynomialDerivativeSize.resizeArray2Dim(*polynomialDerivative, kernelMatrixDim);
            polynomialDerivativeSize.resizeArray2Dim(*polynomialKernelNorm, kernelMatrixDim);
            polynomialKernel = kernelMatrix;
        }
    }
}

void KernelPolynomial::updatePolynomialKernel(const Frame& frame)
{
    updateKernel(frame);
    allocateDerivativeArrays();
    computePolynomialKernelTerms(*frame.get_typeMap());
}

void KernelPolynomial::computePolynomialKernelTerms(const TypeMap& typeMap)
{
    // unboxing variables
    Vec2Real& polynomialDerivative = (*this->polynomialDerivative);
    Vec2Real& polynomialKernel = (*this->polynomialKernel);
    Vec2Real& polynomialKernelNorm = (*this->polynomialKernelNorm);

    // transfer data in a std algorithm form
    std::vector<std::vector<std::tuple<Real&, Real&, Real&>>> data(polynomialDerivative.size());
    for (std::size_t i = 0; i < polynomialKernel.size(); i++)
    {
        for (std::size_t j = 0; j < polynomialKernel[i].size(); j++)
        {
            data[i].push_back(
                {polynomialKernel[i][j], polynomialDerivative[i][j], polynomialKernelNorm[i][j]});
        }
    }

    // actual calculations
    switch (algoExecution)
    {
    case ExecutionPolicy::cpuSingleCore:
        computePolynomialKernelTermsCPU(typeMap, data);
        break;
    case ExecutionPolicy::gpuStdLib:
        VASPML_PARALLEL(
            computePolynomialKernelTermsGPU(typeMap, data);
        );
        break;
    }
}

void KernelPolynomial::computePolynomialKernelTermsCPU(
    const TypeMap&                                            typeMap,
    std::vector<std::vector<std::tuple<Real&, Real&, Real&>>> data)
{
    const Vec1Real& normAtoms = descriptorCollection->get_normAtom();
    const Vec1Int&  nAtomsPerType = (*this->nAtomsPerType);
    const Vec1Int&  numberLocalRefConfs = (*this->numberLocalRefConfs);
    const Int&      kernelPower = (this->kernelPower);
    std::size_t     atomCounter = 0;
    for (std::size_t type = 0; type < typeMap.countStructureTypes(); type++)
    {
        std::size_t typeFF = typeMap.toType(type);
        std::size_t numberLocalRefConfsType = numberLocalRefConfs[typeFF];
        std::size_t increaseCounter = 0;
        for (std::size_t centralAtom = 0; centralAtom < (std::size_t)nAtomsPerType[type];
             centralAtom++)
        {
            Real norm = 0;
            if (normAtoms[atomCounter] > constants::EPS_TOL)
            {
                norm = (Real)kernelPower / normAtoms[atomCounter];
            }
            std::for_each( // par,
                data[type].begin() + increaseCounter,
                data[type].begin() + increaseCounter + (std::size_t)numberLocalRefConfsType,
                [&norm, &kernelPower](std::tuple<Real&, Real&, Real&>& kernelValue)
                {
                    Real factor = math::intPowMixAlgo(std::get<0>(kernelValue), kernelPower - 1);
                    std::get<1>(kernelValue) = factor * norm;
                    std::get<0>(kernelValue) = factor * std::get<0>(kernelValue);
                    std::get<2>(kernelValue) = std::get<0>(kernelValue) * norm;
                });
            atomCounter++;
            increaseCounter += numberLocalRefConfsType;
        }
    }
}

void KernelPolynomial::computePolynomialKernelTermsGPU(
    const TypeMap&                                            typeMap,
    std::vector<std::vector<std::tuple<Real&, Real&, Real&>>> data)
{
    const Vec1Real& normAtoms = descriptorCollection->get_normAtom();
    //const Vec1Int& nAtomsPerType       = (*this -> nAtomsPerType);
    const Vec1Int& numberLocalRefConfs = (*this->numberLocalRefConfs);
    const Int&     kernelPower = (this->kernelPower);

    // central atom indices of SHS2 and SHS3 are the same, so it is hopefully enough to take SHS2
    const Vec1Int& typeIndexCentral =
        descriptorCollection->getDescriptor("SHS2-2-body").get_typeIndexCentral();
    //const std::size_t numberAtoms    =  descriptorCollection -> getDescriptor("SHS2-2-body").
    //                                                            get_neighborList().get_numberAtoms();
    const Vec1Int& centralAtomIndexPerType = descriptorCollection->getDescriptor("SHS2-2-body")
                                                 .get_neighborList()
                                                 .get_centralAtomIndexPerType();
    std::for_each( // seq,
        centralAtomIndex.cbegin(),
        centralAtomIndex.cbegin() + centralAtomIndexSize.actDim,
        [&](const std::size_t& atomCounter) mutable
        {
            std::size_t type = typeIndexCentral[atomCounter];
            std::size_t typeFF = typeMap.toType(type);
            std::size_t numberLocalRefConfsType = numberLocalRefConfs[typeFF];
            std::size_t increaseCounter =
                centralAtomIndexPerType[atomCounter] * numberLocalRefConfsType;
            Real norm = 0;
            if (normAtoms[atomCounter] > constants::EPS_TOL)
            {
                norm = (Real)kernelPower / normAtoms[atomCounter];
            }
            std::for_each( // seq,
                data[type].begin() + increaseCounter,
                data[type].begin() + numberLocalRefConfsType + increaseCounter,
                [&](std::tuple<Real&, Real&, Real&>& kernelValue)
                {
                    Real factor = math::intPowMixAlgo(std::get<0>(kernelValue), kernelPower - 1);
                    std::get<1>(kernelValue) = factor * norm;
                    std::get<0>(kernelValue) = factor * std::get<0>(kernelValue);
                    std::get<2>(kernelValue) = std::get<0>(kernelValue) * norm;
                });
        });
}

void KernelPolynomial::computeDerivativeKernel(Vec1Real&          derivativeMatrix,
                                               Vec1Real&          forceVector,
                                               const Real*        descriptorsRefConfsWeighted,
                                               const Real*        regressionCoefficients,
                                               const std::size_t& typeStruc,
                                               const std::size_t& typeForceField,
                                               const String&      key) const
{

    linalg::matMul(linalg::NoTrans,
                   linalg::NoTrans,
                   (*nAtomsPerType)[typeStruc],
                   (*numberDescriptors.at(key))[typeForceField],
                   (*numberLocalRefConfs)[typeForceField],
                   std::sqrt(weights.at(key)),
                   (*polynomialDerivative)[typeStruc],
                   (*numberLocalRefConfs)[typeForceField],
                   descriptorsRefConfsWeighted,
                   (*numberDescriptors.at(key))[typeForceField],
                   (Real)0.0,
                   derivativeMatrix,
                   (*numberDescriptors.at(key))[typeForceField],
                   linalgContext);

    linalg::matVec(linalg::NoTrans,
                   (*nAtomsPerType)[typeStruc],
                   (*numberLocalRefConfs)[typeForceField],
                   -std::sqrt(weights.at(key)),
                   (*polynomialKernelNorm)[typeStruc],
                   (*numberLocalRefConfs)[typeForceField],
                   regressionCoefficients,
                   1,
                   0.0,
                   forceVector,
                   1,
                   linalgContext);

    const std::map<String, ShVec2Real>& descriptors =
        descriptorCollection->get_descriptorsNormalized();
    for (std::size_t atom = 0; atom < (std::size_t)(*nAtomsPerType)[typeStruc]; atom++)
    {
        Real* derivativeMatrixPtr =
            &derivativeMatrix[atom * (*numberDescriptors.at(key))[typeForceField]];
        const Real* descriptorPtr =
            &(*descriptors.at(key))[typeStruc][atom * (*numberDescriptors.at(key))[typeForceField]];
        linalg::scaleVectorPlusVector(forceVector[atom],
                                      descriptorPtr,
                                      derivativeMatrixPtr,
                                      (*numberDescriptors.at(key))[typeForceField],
                                      linalgContext);
    }
}
