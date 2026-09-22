#include "SpillingFactor.hpp"

#include "debug.hpp"
#include "math.hpp"
#include "utils.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

using namespace vaspml;

namespace vaspml
{
std::shared_ptr<ShmemArray2DVariableLen<Real>> makeInverseCovMatrixArray(
    const IoHandlerML_FF&         inputParameters,
    const std::shared_ptr<MlMPI>& mpiIn)
{

    const Vec1Int& numberlocRef =
        inputParameters["number-local-reference-configs"].dcget<ShVec1Int>();
    Vec1Int locRefSquared(numberlocRef.size());
    std::transform(numberlocRef.begin(),
                   numberlocRef.end(),
                   locRefSquared.begin(),
                   [](const Int& value) { return value * value; });
    std::shared_ptr<ShmemArray2DVariableLen<Real>> array;
    array = std::make_shared<ShmemArray2DVariableLen<Real>>(locRefSquared, mpiIn);
    for (std::size_t type = 0; type < numberlocRef.size(); type++)
    {
        for (std::size_t locRef = 0; locRef < (std::size_t)locRefSquared[type]; locRef++)
        {
#ifdef use_shmem
            if ((*array)["mpiShmem"].get_rank() == 0 and (*array)["mpiInter"].get_rank() == 0)
#endif
                array->set_value(
                    type,
                    locRef,
                    inputParameters["inverse-cov-matrix"].dcget<ShVec2Real>()[type][locRef]);
        }
    }
#ifdef use_shmem
    (*array)["mpiShmem"].barrier();
    array->distributeInternode();
#endif
    return array;
}
} //namespace vaspml

SpillingFactor::SpillingFactor(const IoHandlerML_FF&         inputParameters,
                               const std::shared_ptr<MlMPI>& mpiIn,
                               const ShVec2Real&             spillingFactor) :
    numberLocalRefConfs(inputParameters["number-local-reference-configs"].cget<ShVec1Int>()),
    inverseKernelMatrix(makeInverseCovMatrixArray(inputParameters, mpiIn)),
    isComputed(false),
    statisticsComputed(false),
    averageSpillingFactor(0.0),
    maxSpillingFactor(0.0),
    minSpillingFactor(0.0),
    averageSpillingFactorType(0),
    maxSpillingFactorType(0),
    minSpillingFactorType(0)
{
    if (this->spillingFactor == nullptr) { this->spillingFactor = std::make_shared<Vec2Real>(); }
    else { this->spillingFactor = spillingFactor; }
}

void SpillingFactor::computeInverse(void)
{

    const Vec1Int& numberLocalRefConfs = (*this->numberLocalRefConfs);
    for (std::size_t type = 0; type < numberLocalRefConfs.size(); type++)
    {
        linalg::computeInverseLU(inverseKernelMatrix->get_sliceReference(type),
                                 numberLocalRefConfs[type],
                                 numberLocalRefConfs[type]);
    }
}

void SpillingFactor::computeSpillingFactor(const Vec2Real& kernelMatrix,
                                           const Vec1Int&  nAtomsType,
                                           const TypeMap&  typeMap)
{

    ShmemArray2DVariableLen<Real>& inverseKernelMatrix = (*this->inverseKernelMatrix);
    const Vec1Int&                 numberLocalRefConfs = (*this->numberLocalRefConfs);
    Vec2Real&                      spillingFactor = (*this->spillingFactor);
    // compute the working array sizes for every atom type of MLFF
    const Vec1Int arraySizes = vector_tools::elementWiseProduct(nAtomsType, numberLocalRefConfs);
    Int           maxSize = vector_tools::get_max(arraySizes);
    Vec1Real      workingArray;
    workingArray.reserve(maxSize);
    vector_tools::allocate_vector(spillingFactor, nAtomsType);
    for (std::size_t type = 0; type < nAtomsType.size(); type++)
    {
        std::size_t typeFF = typeMap.toType(type);
        workingArray.resize(arraySizes[typeFF]);
        linalg::matMul(linalg::NoTrans,
                       linalg::NoTrans,
                       nAtomsType[type],
                       numberLocalRefConfs[typeFF],
                       numberLocalRefConfs[typeFF],
                       linalg::one,
                       kernelMatrix[type],
                       numberLocalRefConfs[typeFF],
                       inverseKernelMatrix.get_slice(typeFF),
                       numberLocalRefConfs[typeFF],
                       (Real)0.0,
                       workingArray,
                       numberLocalRefConfs[typeFF],
                       linalgContext);
        for (std::size_t atomInType = 0; atomInType < (std::size_t)nAtomsType[type]; atomInType++)
        {
            const Real* workingArrayPtr = &workingArray[atomInType * numberLocalRefConfs[typeFF]];
            const Real* kernelMatrixPtr =
                &kernelMatrix[type][atomInType * numberLocalRefConfs[typeFF]];
            spillingFactor[type][atomInType] =
                std::min(1.0,
                         std::abs(1.0
                                  - linalg::dotProduct(workingArrayPtr,
                                                       kernelMatrixPtr,
                                                       numberLocalRefConfs[typeFF],
                                                       linalgContext)));
        }
    }
    isComputed = true;
}

void SpillingFactor::computeStatistics(void)
{

    VASPML_DEBUG_L1(
        if (!isComputed)
        {
            throw std::runtime_error("ERROR: SpillingFactor::computeStatistics( void )\n"
                                     "Trying to compute statistics of spilling factor, but "
                                     "spilling factor was not computed yet");
        }
    );
    Vec2Real& spillingFactor = (*this->spillingFactor);
    averageSpillingFactorType.resize(spillingFactor.size());
    maxSpillingFactorType.resize(spillingFactor.size());
    minSpillingFactorType.resize(spillingFactor.size());
    for (std::size_t type = 0; type < averageSpillingFactorType.size(); type++)
    {
        averageSpillingFactorType[type] = math::average(spillingFactor[type]);
        minSpillingFactorType[type] = math::minimum(spillingFactor[type]);
        maxSpillingFactorType[type] = math::maximum(spillingFactor[type]);
    }

    averageSpillingFactor = math::average(averageSpillingFactorType);
    minSpillingFactor = math::minimum(averageSpillingFactorType);
    maxSpillingFactor = math::maximum(averageSpillingFactorType);
    statisticsComputed = true;
}

const Real& SpillingFactor::get_spillingFactor(const std::size_t type, const std::size_t atom) const
{
    return (*spillingFactor)[type][atom];
}

const Vec1Real& SpillingFactor::get_spillingFactor(const std::size_t type) const
{
    return (*spillingFactor)[type];
}

const Vec2Real& SpillingFactor::get_spillingFactor(void) const
{
    return *spillingFactor;
}

void SpillingFactor::writeToScreen(void) const
{
    VASPML_DEBUG_L1(
        if (!isComputed)
        {
            throw std::runtime_error(
                "ERROR: void SpillingFactor::writeToScreen( void ) const.\n"
                "Trying to access spillingFactor array which was not computed yet.\n"
                "Call void SpillingFactor::computeSpillingFactor( const Vec2Real& kernelMatrix, "
                "const Vec1Int& nAtomsType, const TypeMap& typeMap ) first");
        }
    );
    const Vec2Real& spillingFactor = (*this->spillingFactor);
    std::size_t     atom = 0;
    for (std::size_t type = 0; type < spillingFactor.size(); type++)
    {
        for (std::size_t atomInType = 0; atomInType < spillingFactor[type].size(); atomInType++)
        {
            String output =
                str("%-6zu -%6zu %24.16E ", atom, type, spillingFactor[type][atomInType]);
            std::cout << output << std::endl;
            atom++;
        }
    }
}

const Real& SpillingFactor::get_averageSpillingFactor(void) const
{
    VASPML_DEBUG_L1(
        if (!isComputed)
        {
            throw std::runtime_error(
                "ERROR: const Real& SpillingFactor::get_averageSpillingFactor( void ) const\n"
                "Trying to access statistics of spilling factor, but statistics were not computed "
                "yet."
                "Call void SpillingFactor::computeStatistics( void ) first");
        }
    );
    return averageSpillingFactor;
}

const Real& SpillingFactor::get_maxSpillingFactor(void) const
{
    VASPML_DEBUG_L1(
        if (!isComputed)
        {
            throw std::runtime_error(
                "ERROR: const Real& SpillingFactor::get_maxSpillingFactor( void ) const.\n"
                "Trying to access statistics of spilling factor, but statistics were not computed "
                "yet.\n"
                "Call void SpillingFactor::computeStatistics( void ) first");
        }
    );
    return maxSpillingFactor;
}

const Real& SpillingFactor::get_minSpillingFactor(void) const
{
    VASPML_DEBUG_L1(
        if (!isComputed)
        {
            throw std::runtime_error(
                "ERROR: const Real& SpillingFactor::get_minSpillingFactor( void ) const.\n"
                "Trying to access statistics of spilling factor, but statistics were not computed "
                "yet.\n"
                "Call void SpillingFactor::computeStatistics( void ) first");
        }
    );
    return minSpillingFactor;
}

const Vec1Real& SpillingFactor::get_averageSpillingFactorType(void) const
{
    VASPML_DEBUG_L1(
        if (!isComputed)
        {
            throw std::runtime_error(
                "ERROR: Vec1Real& SpllingFactor::get_averageSpillingFactorType( void ) const.\n"
                "Trying to access statistics of spilling factor, but statistics were not computed "
                "yet.\n"
                "Call void SpillingFactor::computeStatistics( void ) first");
        }
    );
    return averageSpillingFactorType;
}

const Vec1Real& SpillingFactor::get_maxSpillingFactorType(void) const
{
    VASPML_DEBUG_L1(
        if (!isComputed)
        {
            throw std::runtime_error(
                "ERROR: Vec1Real& SpllingFactor::get_maxSpillingFactorType( void ) const.\n"
                "Trying to access statistics of spilling factor, but statistics were not computed "
                "yet.\n"
                "Call void SpillingFactor::computeStatistics( void ) first");
        }
    );
    return maxSpillingFactorType;
}

const Vec1Real& SpillingFactor::get_minSpillingFactorType(void) const
{
    return minSpillingFactorType;
}

const Real& SpillingFactor::get_averageSpillingFactorType(const std::size_t type) const
{
    return averageSpillingFactorType[type];
}

const Real& SpillingFactor::get_maxSpillingFactorType(const std::size_t type) const
{
    return maxSpillingFactorType[type];
}

const Real& SpillingFactor::get_minSpillingFactorType(const std::size_t type) const
{
    return minSpillingFactorType[type];
}
