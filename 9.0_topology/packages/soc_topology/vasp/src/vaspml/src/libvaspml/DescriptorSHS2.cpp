#include "DescriptorSHS2.hpp"
#include "ParallelEnvironemt.hpp"
#include "Timer.hpp"
#include "Tutor.hpp"
#include "constants.hpp"
#include "debug.hpp"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <functional>
#include <iostream>
#include <stdexcept>

// debug
#include "utils.hpp"
#include <iomanip>

using namespace vaspml;

DescriptorSHS2::DescriptorSHS2(const Real            weight,
                               const bool            isNormalized,
                               const ShVec2Real&     clnmPair,
                               const ShVec2Real&     clnmPairDerivative,
                               const ShVec2Real&     clnmVasp,
                               const ShVec2Real&     clnmDerivativeCentralVasp,
                               const DescriptorType  descriptorType,
                               const ExecutionPolicy algoExecution) :
    Descriptor(weight, isNormalized, descriptorType, algoExecution)
{
    if (clnmPair == nullptr) { this->clnmPair = std::make_shared<Vec2Real>(); }
    else { this->clnmPair = clnmPair; }
    if (clnmPairDerivative == nullptr) { this->clnmPairDerivative = std::make_shared<Vec2Real>(); }
    else { this->clnmPairDerivative = clnmPairDerivative; }
    if (clnmVasp == nullptr) { this->clnmVasp = std::make_shared<Vec2Real>(); }
    else { this->clnmVasp = clnmVasp; }
    if (clnmDerivativeCentralVasp == nullptr)
    {
        this->clnmDerivativeCentralVasp = std::make_shared<Vec2Real>();
    }
    else { this->clnmDerivativeCentralVasp = clnmDerivativeCentralVasp; }

    pairsComputed = false;
    vaspFormatComputed = false;
    basisFunctionsSet = false;
}

std::size_t DescriptorSHS2::compute_BasisSetSize(const std::size_t               lmax,
                                                 const std::vector<std::size_t>& nRoots)
{

    std::size_t basisSetSize = 0;
    for (std::size_t l = 0; l < lmax; l++)
    {
        for (std::size_t m = 0; m < 2 * l + 1; m++)
        {
            for (std::size_t root = 0; root < nRoots[l]; root++) { basisSetSize++; }
        }
    }
    return basisSetSize;
}

std::size_t DescriptorSHS2::get_basisSetSize(void) const
{
    return basisSetSize;
}

void DescriptorSHS2::updatePairCoefficients(const std::shared_ptr<NearestNeighborNSquare>& nn_list,
                                            const std::shared_ptr<BasisFunctions>& basisFunctions)
{
    if (!nn_list->is_typeSorted())
    {
        String fname = __func__;
        global_scope::tutor.bug(
            "Neighbor list is not type sorted. Error in ComputeCLMs::updateCLMs");
    }
    set_neighborList(nn_list);
    // this should be removed
    nRoots = basisFunctions->get_nValidRoots();
    // plus 1 to be able to write loops with < smbol
    maxOrder = basisFunctions->get_maxOrder() + 1;
    basisSetSize = compute_BasisSetSize(maxOrder, nRoots);
    // allocate pair coefficient arrays
    nRadialBasis = basisFunctions->get_totalNumberBasisFunctions();
    computeOffsets();

    // allocate arrays for radial Basis function
    nAtoms = neighborList->get_numberAtoms();
    allocatePairCoefficientArrays(basisFunctions->get_ldim());
    resize_radialBasisFunction();
    VASPML_PROFILING_START("DescriptorSHS2::updatePairCoefficients");
    switch (algoExecution)
    {
    case ExecutionPolicy::cpuSingleCore:
        updatePairCoefficientsCPU(basisFunctions);
        break;
    case ExecutionPolicy::gpuStdLib:
        VASPML_PARALLEL(
            updatePairCoefficientsGPU(basisFunctions);
        );
        break;
    }
    pairsComputed = true;
    VASPML_PROFILING_STOP("DescriptorSHS2::updatePairCoefficients");
}

void DescriptorSHS2::updateVaspCoefficients(const std::shared_ptr<NearestNeighborNSquare>& nn_list)
{
    if (!nn_list->is_typeSorted())
    {
        String fname = __func__;
        global_scope::tutor.bug(
            "In function " + fname
            + " Neighbor list is not type sorted. Error in ComputeCLMs::updateCLMs");
    }
    if (!basisFunctionsSet)
    {
        String fname = __func__;
        global_scope::tutor.bug("In function " + fname
                                + " basis functions have to be set. Please call  ");
    }

    set_neighborList(nn_list);
    nAtoms = neighborList->get_numberAtoms();
    allocateArraysVaspFormat();
    resize_centralAtomIndex();
    resize_radialBasisFunction();
    resize_ylm(basisFunctions->get_ldim());
    VASPML_PROFILING_START("DescriptorSHS2::updateVaspCoefficients");
    switch (algoExecution)
    {
    case ExecutionPolicy::cpuSingleCore:
        updateVaspCoefficientsCPU();
        break;
    case ExecutionPolicy::gpuStdLib:
        VASPML_PARALLEL(
            updateVaspCoefficientsGPU();
        );
        break;
    }
    VASPML_PROFILING_STOP("DescriptorSHS2::updateVaspCoefficients");
    vaspFormatComputed = true;
}

void DescriptorSHS2::set_basisFunctions(const std::shared_ptr<BasisFunctions>& basisFunctions)
{
    this->basisFunctions = basisFunctions;
    // this should be removed
    nRoots = basisFunctions->get_nValidRoots();
    // plus 1 to be able to write loops with < smbol
    maxOrder = basisFunctions->get_maxOrder() + 1;
    basisSetSize = compute_BasisSetSize(maxOrder, nRoots);
    // allocate pair coefficient arrays
    nRadialBasis = basisFunctions->get_totalNumberBasisFunctions();
    computeOffsets();
    basisFunctionsSet = true;
}

void DescriptorSHS2::updatePairCoefficientsCPU(
    const std::shared_ptr<BasisFunctions>& basisFunctions)
{
    for (std::size_t atom = 0; atom < nAtoms; atom++)
    {
        computePairCoefficientSingleAtom(atom, basisFunctions, *neighborList);
    }
}

void DescriptorSHS2::updatePairCoefficientsGPU(
    const std::shared_ptr<BasisFunctions>& basisFunctions)
{

    const NearestNeighborNSquare& NNList = (*this->neighborList);
    std::for_each(VASPML_SEQ
                  centralAtomIndex.cbegin(),
                  centralAtomIndex.cend(),
                  [&](const std::size_t& atom)
                  { computePairCoefficientSingleAtom(atom, basisFunctions, NNList); });
}

void DescriptorSHS2::updateVaspCoefficientsCPU(void)
{
    for (std::size_t atom = 0; atom < nAtoms; atom++) { computeVaspCoefficientSingleAtom(atom); }
}

void DescriptorSHS2::updateVaspCoefficientsGPU(void)
{

    std::for_each(VASPML_PAR
                  centralAtomIndex.cbegin(),
                  centralAtomIndex.cbegin() + centralAtomIndexSize.actDim,
                  [&](const std::size_t& atom) { computeVaspCoefficientSingleAtom(atom); });
}

void DescriptorSHS2::computeVaspCoefficientSingleAtom(const std::size_t atom)
{

    auto& clnmVasp = (*this->clnmVasp)[atom];
    auto& clnmDerivativeCentralVasp = (*this->clnmDerivativeCentralVasp)[atom];
    // allocate arrays for aradial Basis function

    // compute angular descriptor
    const Vec1Real& distances = neighborList->get_distances(atom);
    const Vec1Real& connectionVectorNormalized = neighborList->get_connectionVectorNormalized(atom);
    const Vec1Int&  neighborTypeIndex = neighborList->get_typeIndex(atom);
    basisFunctions->computeAngularBasis(connectionVectorNormalized,
                                        distances,
                                        ylm[atom],
                                        ylmDerivative[atom]);
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    std::size_t nbasis = 0;
    std::size_t ylm_index = 0;
    std::size_t ylm_index_xyz = 0;
    std::size_t ylm_index_offset = 0;
    std::size_t ylm_index_xyz_offset = 0;

    std::size_t neighbor = 0;
    // compute radial basis functions
    neighbor = 0;
    std::for_each(
        VASPML_SEQ
        distances.cbegin(),
        distances.cend(),
        [&](const Real& r) mutable
        {
            // compute radial basis for atom pair
            basisFunctions->interpolate(r, radialBasisFunction, radialBasisFunctionDerivative);
            Real normed_x = connectionVectorNormalized[3 * neighbor];
            Real normed_y = connectionVectorNormalized[3 * neighbor + 1];
            Real normed_z = connectionVectorNormalized[3 * neighbor + 2];

            Size_t besselOffset = 0;
            Size_t nderivativeBasis = 0;
            Size_t typeOffset = neighborTypeIndex[neighbor] * basisSetSize;
            Size_t typeOffset3x = 3 * neighborTypeIndex[neighbor] * basisSetSize;
            for (std::size_t ll = 0; ll < maxOrder; ll++)
            {
                for (std::size_t nRoot = 0; nRoot < nRoots[ll]; nRoot++)
                {
                    ylm_index = ylm_index_offset;
                    ylm_index_xyz = ylm_index_xyz_offset;
                    for (std::size_t mm = 0; mm < 2 * ll + 1; mm++)
                    {
                        // compute expansion coefficient
                        clnmVasp[typeOffset + besselOffset + nRoot] +=
                            radialBasisFunction[besselOffset + nRoot] * ylm[atom][ylm_index];
                        // compute derivatives
                        Real x_component = ylm[atom][ylm_index]
                                             * radialBasisFunctionDerivative[besselOffset + nRoot]
                                             * normed_x
                                         + ylmDerivative[atom][ylm_index_xyz]
                                               * radialBasisFunction[besselOffset + nRoot];

                        Real y_component = ylm[atom][ylm_index]
                                             * radialBasisFunctionDerivative[besselOffset + nRoot]
                                             * normed_y
                                         + ylmDerivative[atom][ylm_index_xyz + 1]
                                               * radialBasisFunction[besselOffset + nRoot];

                        Real z_component = ylm[atom][ylm_index]
                                             * radialBasisFunctionDerivative[besselOffset + nRoot]
                                             * normed_z
                                         + ylmDerivative[atom][ylm_index_xyz + 2]
                                               * radialBasisFunction[besselOffset + nRoot];
                        clnmDerivativeCentralVasp[typeOffset3x + nderivativeBasis] += x_component;
                        clnmDerivativeCentralVasp[typeOffset3x + basisSetSize + nderivativeBasis] +=
                            y_component;
                        clnmDerivativeCentralVasp[typeOffset3x + 2 * basisSetSize
                                                  + nderivativeBasis] += z_component;
                        nbasis++;
                        nderivativeBasis++;
                        ylm_index++;
                        ylm_index_xyz += 3;
                    }
                }
                ylm_index_offset += 2 * ll + 1;
                ylm_index_xyz_offset += 6 * ll + 3;
                besselOffset += nRoots[ll];
            }
            neighbor++;
        });
}

void DescriptorSHS2::computePairCoefficientSingleAtom(
    const std::size_t                      atom,
    const std::shared_ptr<BasisFunctions>& basisFunctions,
    const NearestNeighborNSquare&          NNL)
{

    auto& clnmPair = *(this->clnmPair);
    auto& clnmPairDerivative = *(this->clnmPairDerivative);

    // allocate arrays for aradial Basis function

    // compute angular descriptor
    const Vec1Real& distances = NNL.get_distances(atom);
    const Vec1Real& connectionVectorNormalized = NNL.get_connectionVectorNormalized(atom);
    basisFunctions->computeAngularBasis(connectionVectorNormalized,
                                        distances,
                                        ylm[atom],
                                        ylmDerivative[atom]);
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    std::size_t nbasis = 0;
    std::size_t ylm_index = 0;
    std::size_t ylm_index_xyz = 0;
    std::size_t ylm_index_offset = 0;
    std::size_t ylm_index_xyz_offset = 0;

    std::size_t neighbor = 0;
    // compute radial basis functions
    neighbor = 0;
    std::for_each(
        VASPML_SEQ
        distances.cbegin(),
        distances.cend(),
        [&](const Real& r) mutable
        {
            // compute radial basis for atom pair
            basisFunctions->interpolate(r, radialBasisFunction, radialBasisFunctionDerivative);
            Real normed_x = connectionVectorNormalized[3 * neighbor];
            Real normed_y = connectionVectorNormalized[3 * neighbor + 1];
            Real normed_z = connectionVectorNormalized[3 * neighbor + 2];

            std::size_t besselOffset = 0;
            std::size_t nderivativeBasis = 0;
            for (std::size_t ll = 0; ll < maxOrder; ll++)
            {
                for (std::size_t nRoot = 0; nRoot < nRoots[ll]; nRoot++)
                {
                    ylm_index = ylm_index_offset;
                    ylm_index_xyz = ylm_index_xyz_offset;
                    for (std::size_t mm = 0; mm < 2 * ll + 1; mm++)
                    {
                        // compute expansion coefficient
                        clnmPair[atom][nbasis] =
                            radialBasisFunction[besselOffset + nRoot] * ylm[atom][ylm_index];
                        // compute derivatives
                        Real x_component = ylm[atom][ylm_index]
                                             * radialBasisFunctionDerivative[besselOffset + nRoot]
                                             * normed_x
                                         + ylmDerivative[atom][ylm_index_xyz]
                                               * radialBasisFunction[besselOffset + nRoot];

                        Real y_component = ylm[atom][ylm_index]
                                             * radialBasisFunctionDerivative[besselOffset + nRoot]
                                             * normed_y
                                         + ylmDerivative[atom][ylm_index_xyz + 1]
                                               * radialBasisFunction[besselOffset + nRoot];

                        Real z_component = ylm[atom][ylm_index]
                                             * radialBasisFunctionDerivative[besselOffset + nRoot]
                                             * normed_z
                                         + ylmDerivative[atom][ylm_index_xyz + 2]
                                               * radialBasisFunction[besselOffset + nRoot];

                        clnmPairDerivative[atom][3 * neighbor * basisSetSize + nderivativeBasis] =
                            x_component;
                        clnmPairDerivative[atom][3 * neighbor * basisSetSize + basisSetSize
                                                 + nderivativeBasis] = y_component;
                        clnmPairDerivative[atom][3 * neighbor * basisSetSize + 2 * basisSetSize
                                                 + nderivativeBasis] = z_component;

                        nbasis++;
                        nderivativeBasis++;
                        ylm_index++;
                        ylm_index_xyz += 3;
                    }
                }
                ylm_index_offset += 2 * ll + 1;
                ylm_index_xyz_offset += 6 * ll + 3;
                besselOffset += nRoots[ll];
            }
            neighbor++;
        });
}

void DescriptorSHS2::computeVaspCoefficientsFromPairCoefficients(void)
{
    if (!pairsComputed)
    {
        throw std::runtime_error(
            "ERROR in routine DescriptorSHS2::computeVaspCoefficientsFromPairCoefficients \n"
            "Pairs have to be computed first. Adapt your code such that \n"
            "void DescriptorSHS2::updatePairCoefficients( const "
            "std::shared_ptr<NearestNeighborNSquare>& nn_list,const "
            "std::unique_ptr<DescriptorBase>& descriptor )\n"
            "is called before \n");
    }

    VASPML_PROFILING_START("DescriptorSHS2::computeVaspCoefficientsFromPairCoefficients");
    allocateArraysVaspFormat();
    switch (algoExecution)
    {
    case ExecutionPolicy::cpuSingleCore:
        computeVaspCoefficientsFromPairCoefficientsCPU();
        break;
    case ExecutionPolicy::gpuStdLib:
        VASPML_PARALLEL(
            computeVaspCoefficientsFromPairCoefficientsGPU();
        );
        break;
    }
    VASPML_PROFILING_STOP("DescriptorSHS2::computeVaspCoefficientsFromPairCoefficients");
}

void DescriptorSHS2::computeVaspCoefficientsFromPairCoefficientsCPU(void)
{
    auto& clnmPair = *(this->clnmPair);
    auto& clnmPairDerivative = *(this->clnmPairDerivative);
    auto& clnmVasp = *(this->clnmVasp);
    auto& clnmDerivativeCentralVasp = *(this->clnmDerivativeCentralVasp);

    for (std::size_t index_atom = 0; index_atom < neighborList->get_nAtoms(); index_atom++)
    {
        for (std::size_t neighbor = 0; neighbor < neighborList->get_size(index_atom); neighbor++)
        {
            const auto& clnmPair_begin = clnmPair[index_atom].begin() + neighbor * basisSetSize;
            const auto& clnm_derivative_pair_begin =
                clnmPairDerivative[index_atom].begin() + 3 * neighbor * basisSetSize;
            const std::size_t& neighborType = neighborList->get_typeIndex(index_atom, neighbor);
            std::size_t        typeOffset = neighborType * basisSetSize;
            std::size_t        typeOffset3x = 3 * neighborType * basisSetSize;

            std::transform( // par_unseq,
                clnmPair_begin,
                clnmPair_begin + basisSetSize,
                clnmVasp[index_atom].begin() + typeOffset,
                clnmVasp[index_atom].begin() + typeOffset,
                std::plus<Real>());
            std::transform( // par_unseq,
                clnmDerivativeCentralVasp[index_atom].begin() + typeOffset3x,
                clnmDerivativeCentralVasp[index_atom].begin() + typeOffset3x + basisSetSize,
                clnm_derivative_pair_begin,
                clnmDerivativeCentralVasp[index_atom].begin() + typeOffset3x,
                std::minus<Real>());

            std::transform( // par_unseq,
                clnmDerivativeCentralVasp[index_atom].begin() + typeOffset3x + basisSetSize,
                clnmDerivativeCentralVasp[index_atom].begin() + typeOffset3x + 2 * basisSetSize,
                clnm_derivative_pair_begin + basisSetSize,
                clnmDerivativeCentralVasp[index_atom].begin() + typeOffset3x + basisSetSize,
                std::minus<Real>());

            std::transform( // par_unseq,
                clnmDerivativeCentralVasp[index_atom].begin() + typeOffset3x + 2 * basisSetSize,
                clnmDerivativeCentralVasp[index_atom].begin() + typeOffset3x + 3 * basisSetSize,
                clnm_derivative_pair_begin + 2 * basisSetSize,
                clnmDerivativeCentralVasp[index_atom].begin() + typeOffset3x + 2 * basisSetSize,
                std::minus<Real>());
        }
    }
    vaspFormatComputed = true;
}

void DescriptorSHS2::computeVaspCoefficientsFromPairCoefficientsGPU(void)
{

    std::for_each(VASPML_PAR_UNSEQ
                  centralAtomIndex.cbegin(),
                  centralAtomIndex.cbegin() + nAtoms,
                  [&](const std::size_t& index_atom) mutable
                  {
                      auto& clnmPair = (*this->clnmPair)[index_atom];
                      auto& clnmPairDerivative = (*this->clnmPairDerivative)[index_atom];
                      auto& clnmVasp = (*this->clnmVasp)[index_atom];
                      auto& clnmDerivativeCentralVasp =
                          (*this->clnmDerivativeCentralVasp)[index_atom];
                      const Vec1Int& neighborType = neighborList->get_typeIndex(index_atom);
                      for (std::size_t neighbor = 0; neighbor < neighborList->get_size(index_atom);
                           neighbor++)
                      {
                          const auto& clnmPair_begin = clnmPair.begin() + neighbor * basisSetSize;
                          const auto& clnm_derivative_pair_begin =
                              clnmPairDerivative.begin() + 3 * neighbor * basisSetSize;
                          std::size_t typeOffset = neighborType[neighbor] * basisSetSize;
                          std::size_t typeOffset3x = 3 * neighborType[neighbor] * basisSetSize;

                          std::transform( // par_unseq,
                              clnmPair_begin,
                              clnmPair_begin + basisSetSize,
                              clnmVasp.begin() + typeOffset,
                              clnmVasp.begin() + typeOffset,
                              std::plus<Real>());
                          std::transform( // par_unseq,
                              clnmDerivativeCentralVasp.begin() + typeOffset3x,
                              clnmDerivativeCentralVasp.begin() + typeOffset3x + basisSetSize,
                              clnm_derivative_pair_begin,
                              clnmDerivativeCentralVasp.begin() + typeOffset3x,
                              std::minus<Real>());

                          std::transform( // par_unseq,
                              clnmDerivativeCentralVasp.begin() + typeOffset3x + basisSetSize,
                              clnmDerivativeCentralVasp.begin() + typeOffset3x + 2 * basisSetSize,
                              clnm_derivative_pair_begin + basisSetSize,
                              clnmDerivativeCentralVasp.begin() + typeOffset3x + basisSetSize,
                              std::minus<Real>());

                          std::transform( // par_unseq,
                              clnmDerivativeCentralVasp.begin() + typeOffset3x + 2 * basisSetSize,
                              clnmDerivativeCentralVasp.begin() + typeOffset3x + 3 * basisSetSize,
                              clnm_derivative_pair_begin + 2 * basisSetSize,
                              clnmDerivativeCentralVasp.begin() + typeOffset3x + 2 * basisSetSize,
                              std::minus<Real>());
                      }
                  });

    vaspFormatComputed = true;
}

void DescriptorSHS2::computeThermodynamicIntegration(ThermodynIntegration& thermodynVars)
{

    const auto& atomList = thermodynVars.get_atomList();
    const auto& couplingConstants = thermodynVars.get_couplingConstants();

    for (std::size_t i = 0; i < atomList.size(); i++)
    {
        std::size_t centralAtom = atomList[i];
        std::fill((*clnmVasp)[centralAtom].begin(), (*clnmVasp)[centralAtom].end(), (Real)0);
        for (std::size_t neighbor = 0; neighbor < neighborList->get_size(centralAtom); neighbor++)
        {
            // compute basis set offset
            std::size_t pair_offset = neighbor * basisSetSize;
            // get neighbor index and type
            const Int&         neighborIndex = neighborList->get_globalIndex(centralAtom, neighbor);
            const std::size_t& neighborType = neighborList->get_typeIndex(centralAtom, neighbor);
            // type offset in clnmVasp
            std::size_t typeOffset = neighborType * basisSetSize;
            for (std::size_t nbasis = 0; nbasis < basisSetSize; nbasis++)
            {
                thermodynVars.set_clnmCoupling(centralAtom,
                                               typeOffset + nbasis,
                                               (*clnmVasp)[centralAtom][typeOffset + nbasis]);
                (*clnmVasp)[centralAtom][typeOffset + nbasis] +=
                    couplingConstants[neighborIndex] * (*clnmPair)[i][pair_offset + nbasis];
            }
        }
        // rescale parts which only depend on central atom coupling constant
        std::transform((*clnmPairDerivative)[centralAtom].begin(),
                       (*clnmPairDerivative)[centralAtom].end(),
                       (*clnmPairDerivative)[centralAtom].begin(),
                       std::bind(std::multiplies<Real>(),
                                 std::placeholders::_1,
                                 couplingConstants[centralAtom]));

        std::transform((*clnmDerivativeCentralVasp)[centralAtom].begin(),
                       (*clnmDerivativeCentralVasp)[centralAtom].end(),
                       (*clnmDerivativeCentralVasp)[centralAtom].begin(),
                       std::bind(std::multiplies<Real>(),
                                 std::placeholders::_1,
                                 couplingConstants[centralAtom]));
    }
}

void DescriptorSHS2::resize_centralAtomIndex(void)
{
    const std::size_t& nAtoms = neighborList->get_numberAtoms();
    bool               resize = centralAtomIndexSize.checkResize(nAtoms);
    if (resize)
    {
        centralAtomIndex.resize(nAtoms);
        std::iota(centralAtomIndex.begin(), centralAtomIndex.end(), 0);
    }
}

void DescriptorSHS2::resize_radialBasisFunction(void)
{
    // check if number of radial basis functions changed
    bool resize = radialBasisFunctionSize.checkResize(nRadialBasis);
    if (resize)
    {
        radialBasisFunction.resize(nRadialBasis);
        radialBasisFunctionDerivative.resize(nRadialBasis);
    }
}

void DescriptorSHS2::resize_ylm(const Size_t& ldim)
{
    const std::size_t& nAtoms = neighborList->get_numberAtoms();
    const Vec1Int&     numberNeighbors = neighborList->get_size();
    bool               resize = ylmSize.checkResize1Dim(nAtoms);
    if (resize)
    {
        ylmSize.act1Dim = nAtoms;
        ylmSize.max1Dim = nAtoms;
        ylmSize.resizeArray1Dim(ylm, numberNeighbors, ldim);
        ylmDerivativeSize.act1Dim = nAtoms;
        ylmDerivativeSize.max1Dim = nAtoms;
        ylmDerivativeSize.resizeArray1Dim(ylmDerivative, numberNeighbors, 3 * ldim);
    }
    else
    {
        resize = ylmSize.checkResize2Dim(numberNeighbors, ldim);
        if (resize)
        {
            ylmSize.resizeArray2Dim(ylm, numberNeighbors, ldim);
            ylmDerivativeSize.resizeArray2Dim(ylmDerivative, numberNeighbors, 3 * ldim);
        }
    }
}

void DescriptorSHS2::resizePairArraysGPU(const Size_t& ldim)
{

    const std::size_t& nAtoms = neighborList->get_numberAtoms();
    const Vec1Int&     numberNeighbors = neighborList->get_size();
    bool               resize = clnmPairSize.checkResize1Dim(nAtoms);
    // number of atoms increased. Reset all arrays set new max Size
    if (resize)
    {
        clnmPairSize.resizeArray1Dim(*clnmPair, numberNeighbors, basisSetSize);
        clnmPairDerivativeSize.act1Dim = nAtoms;
        clnmPairDerivativeSize.max1Dim = nAtoms;
        clnmPairDerivativeSize.resizeArray1Dim(*clnmPairDerivative,
                                               numberNeighbors,
                                               3 * basisSetSize);
        ylmSize.act1Dim = nAtoms;
        ylmSize.max1Dim = nAtoms;
        ylmSize.resizeArray1Dim(ylm, numberNeighbors, ldim);
        ylmDerivativeSize.act1Dim = nAtoms;
        ylmDerivativeSize.max1Dim = nAtoms;
        ylmDerivativeSize.resizeArray1Dim(ylmDerivative, numberNeighbors, 3 * ldim);
        centralAtomIndexSize.checkResize(nAtoms);
        centralAtomIndex.resize(nAtoms);
        std::iota(centralAtomIndex.begin(), centralAtomIndex.end(), 0);
    }
    else
    {
        // set new actual sizes of first dimension.
        clnmPairDerivativeSize.act1Dim = nAtoms;
        ylmSize.act1Dim = nAtoms;
        ylmDerivativeSize.act1Dim = nAtoms;
        // check second dimension which means number of neighbors changed or
        // basis set size; Check already valid for all arrays
        resize = clnmPairSize.checkResize2Dim(numberNeighbors, basisSetSize);
        if (resize)
        {
            clnmPairSize.resizeArray2Dim(*clnmPair, numberNeighbors, basisSetSize);
            clnmPairDerivativeSize.resizeArray2Dim(*clnmPairDerivative,
                                                   numberNeighbors,
                                                   3 * basisSetSize);
            ylmSize.resizeArray2Dim(ylm, numberNeighbors, ldim);
            ylmDerivativeSize.resizeArray2Dim(ylmDerivative, numberNeighbors, 3 * ldim);
        }
    }
}

void DescriptorSHS2::allocatePairCoefficientArrays(const Size_t& ldim)
{
    if (algoExecution == ExecutionPolicy::cpuSingleCore)
    {
        const std::size_t& nAtoms = neighborList->get_numberAtoms();
        auto&              clnmPair = *(this->clnmPair);
        auto&              clnmPairDerivative = *(this->clnmPairDerivative);
        clnmPair.resize(nAtoms);
        clnmPairDerivative.resize(nAtoms);
        ylm.resize(nAtoms);
        ylmDerivative.resize(nAtoms);
        radialBasisFunction.resize(nRadialBasis);
        radialBasisFunctionDerivative.resize(nRadialBasis);
        for (std::size_t atom = 0; atom < nAtoms; atom++)
        {
            const Size_t& numberNeighbors = neighborList->get_size(atom);
            clnmPair[atom].resize(basisSetSize * numberNeighbors);
            clnmPairDerivative[atom].resize(3 * numberNeighbors * basisSetSize);
            ylm[atom].resize(numberNeighbors * ldim);
            ylmDerivative[atom].resize(3 * numberNeighbors * ldim);
        }
    }
    else if (algoExecution == ExecutionPolicy::gpuStdLib) resizePairArraysGPU(ldim);
}

void DescriptorSHS2::resizeArraysVaspFormatGPU(void)
{
    const Size_t  nAtoms = neighborList->get_nAtoms();
    const Size_t& numberTypes = neighborList->get_numberTypes();
    bool          resize = clnmVaspSize.checkResize1Dim(nAtoms);
    // number of atoms increased. Reset all arrays set new max Size
    if (resize)
    {
        clnmVaspSize.resizeArray1Dim(*clnmVasp, basisSetSize * numberTypes);
        clnmDerivativeCentralVaspSize.act1Dim = nAtoms;
        clnmDerivativeCentralVaspSize.max1Dim = nAtoms;
        clnmDerivativeCentralVaspSize.resizeArray1Dim(*clnmDerivativeCentralVasp,
                                                      3 * basisSetSize * numberTypes);
    }
    else
    {
        clnmDerivativeCentralVaspSize.act1Dim = nAtoms;
        resize = clnmVaspSize.checkResize2Dim(basisSetSize * numberTypes);
        if (resize)
        {

            clnmVaspSize.resizeArray2Dim(*clnmVasp, basisSetSize * numberTypes);
            clnmDerivativeCentralVaspSize.resizeArray1Dim(*clnmDerivativeCentralVasp,
                                                          3 * basisSetSize * numberTypes);
        }
    }
    // send data back to GPU if it was resized. And rezero
    std::for_each(VASPML_PAR_UNSEQ
                  clnmVasp->begin(),
                  clnmVasp->begin() + clnmVaspSize.act1Dim,
                  [&](Vec1Real& slice)
                  { std::fill(slice.begin(), slice.begin() + clnmVaspSize.actSizeRec, (Real)0); });
    std::for_each(VASPML_PAR_UNSEQ
                  clnmDerivativeCentralVasp->begin(),
                  clnmDerivativeCentralVasp->begin() + clnmDerivativeCentralVaspSize.act1Dim,
                  [&](Vec1Real& slice)
                  {
                      std::fill(slice.begin(),
                                slice.begin() + clnmDerivativeCentralVaspSize.actSizeRec,
                                (Real)0);
                  });
}

void DescriptorSHS2::allocateArraysVaspFormat(void)
{
    if (algoExecution == ExecutionPolicy::cpuSingleCore)
    {
        const Size_t& numberTypes = neighborList->get_numberTypes();
        auto&         clnmVasp = *(this->clnmVasp);
        auto&         clnmDerivativeCentralVasp = *(this->clnmDerivativeCentralVasp);
        clnmVasp.resize(neighborList->get_nAtoms());
        clnmDerivativeCentralVasp.resize(neighborList->get_nAtoms());
        for (std::size_t i = 0; i < neighborList->get_nAtoms(); i++)
        {
            clnmVasp[i].resize(numberTypes * basisSetSize);
            std::fill(clnmVasp[i].begin(), clnmVasp[i].end(), (Real)0.0);
            clnmDerivativeCentralVasp[i].resize(3 * basisSetSize * numberTypes);
            std::fill(clnmDerivativeCentralVasp[i].begin(),
                      clnmDerivativeCentralVasp[i].end(),
                      (Real)0.0);
        }
    }
    else if (algoExecution == ExecutionPolicy::gpuStdLib) resizeArraysVaspFormatGPU();
}

void DescriptorSHS2::computeOffsets(void)
{

    lOffset.resize(maxOrder);
    for (std::size_t ll = 1; ll < maxOrder; ll++)
    {
        lOffset[ll] = lOffset[ll - 1] + (2 * (ll - 1) + 1) * nRoots[ll - 1];
    }
}

std::size_t DescriptorSHS2::get_lOffset(const std::size_t l) const
{
    VASPML_DEBUG_L1(
       if (l >= maxOrder)
       {    
           throw std::runtime_error ("ERROR in const std::size_t DescriptorSHS2::get_lOffset" 
                 "(const std::size_t l) const l index out of bounds");
       }
    );
    return lOffset[l];
}

std::size_t DescriptorSHS2::get_maxOrder(void) const
{
    return maxOrder - 1;
}

const std::size_t& DescriptorSHS2::get_nRootsOrder(const std::size_t l) const
{
    return nRoots[l];
}

const std::vector<std::size_t>& DescriptorSHS2::get_nRoots(void) const
{
    return nRoots;
}

const Real& DescriptorSHS2::get_clnmVasp(const std::size_t central_atom,
                                         const std::size_t type_number,
                                         const std::size_t l,
                                         const std::size_t n,
                                         const std::size_t m) const
{

    VASPML_DEBUG_L1(
      if ( !vaspFormatComputed ){
         throw std::runtime_error( "DescriptorSHS2::get_clnmVasp vasp format coefficients not computed \n"
                                   "please call computeVaspCoefficientsFromPairCoefficients() beforehand"
         );
      }
      if ( central_atom >= nAtoms ){
         throw std::runtime_error(
               "DescriptorSHS2::get_clnmVasp index central_atom out of bounds" );
      }
      if ( l >= maxOrder ){
         throw std::runtime_error(
               "DescriptorSHS2::get_clnmVasp index l out of bounds" );
      }
      if ( n >= nRoots[ l ] ){
         throw std::runtime_error(
               "DescriptorSHS2::get_clnmVasp index n out of bounds" );
      }
      if ( m >= 2*l+1 ){
         throw std::runtime_error(
               "DescriptorSHS2::get_clnmVasp index m out of bounds" );
      }
      if ( !pairsComputed ){
         throw std::runtime_error(
               "DescriptorSHS2::get_clnmVasp vasp pairs are not computed" );
      }
   );
    return (*clnmVasp)[central_atom][get_Index(type_number, l, n, m)];
}

VASPML_NV_HOST_DEVICE
Size_t DescriptorSHS2::get_Index(const std::size_t atom,
                                 const std::size_t l,
                                 const std::size_t n,
                                 const std::size_t m) const
{
    return atom * basisSetSize + lOffset[l] + n * (2 * l + 1) + m;
}

const Vec1Real& DescriptorSHS2::get_clnmPair(const std::size_t centralAtom) const
{
    return (*clnmPair)[centralAtom];
}

const Vec1Real& DescriptorSHS2::get_clnmPairDerivative(const std::size_t centralAtom) const
{
    return (*clnmPairDerivative)[centralAtom];
}

VASPML_NV_HOST_DEVICE
const Vec1Real& DescriptorSHS2::get_clnmVasp(const std::size_t centralAtom) const
{
    return (*clnmVasp)[centralAtom];
}

VASPML_NV_HOST_DEVICE
const Vec1Real& DescriptorSHS2::get_descriptor(const std::size_t centralAtom) const
{
    return (*clnmVasp)[centralAtom];
}

VASPML_NV_HOST_DEVICE
const Vec2Real& DescriptorSHS2::get_descriptor(void) const
{
    return *clnmVasp;
}

VASPML_NV_HOST_DEVICE
Int DescriptorSHS2::get_sizeDescriptor(const std::size_t centralAtom) const
{
    return (*clnmVasp)[centralAtom].size();
}

void DescriptorSHS2::rescale_descriptor(const std::size_t centralAtom, const Real scaleFactor)
{
    std::transform( // par_unseq,
        (*clnmVasp)[centralAtom].begin(),
        (*clnmVasp)[centralAtom].end(),
        (*clnmVasp)[centralAtom].begin(),
        std::bind(std::multiplies<Real>(), std::placeholders::_1, scaleFactor));
}

const Vec1Real& DescriptorSHS2::get_clnmDerivativeCentralVasp(const std::size_t centralAtom)
{
    return (*clnmDerivativeCentralVasp)[centralAtom];
}

std::size_t DescriptorSHS2::computeIndexDerivativeMatrix(std::size_t centralAtom,
                                                         std::size_t typeIndexForceField,
                                                         std::size_t radialIndex,
                                                         std::size_t totalShift)
{

    return centralAtom * totalShift + typeIndexForceField * totalShift + radialIndex;
}

void DescriptorSHS2::compute_forcePreContract(const Vec2Real& derivativeMatrix,
                                              Vec2Real&       forcePreContract,
                                              const TypeMap&  typeMap) const
{
    VASPML_PROFILING_START("DescriptorSHS2::compute_forcePreContract");

    // switch ( algoExecution ){
    //    case ExecutionPolicy::cpuSingleCore:
    compute_forcePreContractCPU(derivativeMatrix, forcePreContract, typeMap);
    //      break;
    //   case ExecutionPolicy::gpuStdLib:
    //      VASPML_PARALLEL( compute_forcePreContractGPU( derivativeMatrix,
    //                                                    forcePreContract,
    //                                                    typeMap );
    //                     );
    //      break;
    // }

    VASPML_PROFILING_STOP("DescriptorSHS2::compute_forcePreContract");
}

void DescriptorSHS2::compute_forcePreContractCPU(const Vec2Real& derivativeMatrix,
                                                 Vec2Real&       forcePreContract,
                                                 const TypeMap&  typeMap) const
{
    const Vec1Int& typeIndexCentral = neighborList->get_typeIndexCentral();
    const Vec1Int& get_centralAtomIndexPerType = neighborList->get_centralAtomIndexPerType();

    std::vector<std::size_t> structureTypes;
    structureTypes.resize(typeMap.countStructureTypes());
    std::iota(structureTypes.begin(), structureTypes.end(), 0);
    std::size_t atomCounter = 0;
    std::size_t atomShift = basisSetSize * typeMap.countForceFieldTypes();

    std::for_each(forcePreContract.begin(),
                  forcePreContract.begin() + nAtoms,
                  [&](Vec1Real& slice)
                  {
                      std::size_t centralType = typeIndexCentral[atomCounter];
                      std::size_t atomInType = get_centralAtomIndexPerType[atomCounter];

                      std::for_each(structureTypes.cbegin(),
                                    structureTypes.cend(),
                                    [&](std::size_t typeInStructure)
                                    {
                                        Int neighborTypeFF = typeMap.toType(typeInStructure);
                                        for (std::size_t nDesc = 0; nDesc < basisSetSize; nDesc++)
                                        {
                                            std::size_t indexMatrix = atomInType * atomShift
                                                                    + neighborTypeFF * basisSetSize
                                                                    + nDesc;
                                            slice[get_Index(typeInStructure, 0, nDesc, 0)] =
                                                derivativeMatrix[centralType][indexMatrix];
                                        }
                                    });

                      atomCounter++;
                  });
}

void DescriptorSHS2::compute_forcePreContractGPU(const Vec2Real& derivativeMatrix,
                                                 Vec2Real&       forcePreContract,
                                                 const TypeMap&  typeMap) const
{

    const Vec1Int& typeIndexCentral = neighborList->get_typeIndexCentral();
    const Vec1Int& centralAtomIndexPerType = neighborList->get_centralAtomIndexPerType();

    std::vector<std::size_t> structureTypes;
    structureTypes.resize(typeMap.countStructureTypes());
    std::iota(structureTypes.begin(), structureTypes.end(), 0);
    std::size_t atomShift = basisSetSize * typeMap.countForceFieldTypes();

    std::for_each( // seq,
        centralAtomIndex.cbegin(),
        centralAtomIndex.cbegin() + nAtoms,
        [&](const std::size_t& atomCounter) mutable
        {
            std::size_t centralType = typeIndexCentral[atomCounter];
            std::size_t atomInType = centralAtomIndexPerType[atomCounter];
            std::for_each(
                structureTypes.cbegin(),
                structureTypes.cend(),
                [&](std::size_t typeInStructure)
                {
                    Int neighborTypeFF = typeMap.toType(typeInStructure);
                    for (std::size_t nDesc = 0; nDesc < basisSetSize; nDesc++)
                    {
                        std::size_t indexMatrix =
                            atomInType * atomShift + neighborTypeFF * basisSetSize + nDesc;
                        forcePreContract[atomCounter][get_Index(typeInStructure, 0, nDesc, 0)] =
                            derivativeMatrix[centralType][indexMatrix];
                    }
                });
        });
}

std::size_t DescriptorSHS2::get_forcePreContractSize(const std::size_t atomIndex) const
{
    return (*clnmVasp)[atomIndex].size();
}

const std::vector<std::size_t>& DescriptorSHS2::get_lOffset(void) const
{
    return lOffset;
}

void DescriptorSHS2::computeForceTerms(const Vec2Real& forcePreContract,
                                       Vec2Real&       pairForces,
                                       Vec1Real&       centralForces) const
{
    VASPML_PROFILING_START("DescriptorSHS2::computeForceTerms");

    switch (algoExecution)
    {
    case ExecutionPolicy::cpuSingleCore:
        computeForceTermsCPU(forcePreContract, pairForces, centralForces);
        break;
    case ExecutionPolicy::gpuStdLib:
        computeForceTermsGPU(forcePreContract, pairForces, centralForces);
        break;
    }

    VASPML_PROFILING_STOP("DescriptorSHS2::computeForceTerms");
}

void DescriptorSHS2::computeForceTermsCPU(const Vec2Real& forcePreContract,
                                          Vec2Real&       pairForces,
                                          Vec1Real&       centralForces) const
{

    std::vector<std::size_t> structureTypes(neighborList->get_numberTypes());
    std::iota(structureTypes.begin(), structureTypes.end(), 0);
    const Vec2Int& neighborTypes = neighborList->get_typeIndex();
    for (std::size_t atom = 0; atom < neighborList->get_nAtoms(); atom++)
    {
        std::for_each(
            structureTypes.begin(),
            structureTypes.end(),
            [&](const std::size_t neighborType)
            {
                // x-component
                const Real* forcePreContractPtr =
                    &forcePreContract[atom][basisSetSize * neighborType];
                const Real* shecPtr =
                    &(*clnmDerivativeCentralVasp)[atom][basisSetSize * neighborType * 3];
                centralForces[3 * atom] -=
                    linalg::dotProduct(shecPtr, forcePreContractPtr, basisSetSize, linalgContext);
                // y-component
                forcePreContractPtr = &forcePreContract[atom][basisSetSize * neighborType];
                shecPtr = &(*clnmDerivativeCentralVasp)[atom][basisSetSize * neighborType * 3
                                                              + basisSetSize];
                centralForces[3 * atom + 1] -=
                    linalg::dotProduct(shecPtr, forcePreContractPtr, basisSetSize, linalgContext);
                // z-component
                forcePreContractPtr = &forcePreContract[atom][basisSetSize * neighborType];
                shecPtr = &(*clnmDerivativeCentralVasp)[atom][basisSetSize * neighborType * 3
                                                              + 2 * basisSetSize];
                centralForces[3 * atom + 2] -=
                    linalg::dotProduct(shecPtr, forcePreContractPtr, basisSetSize, linalgContext);
            });
        std::size_t counter = 0;
        std::size_t neighborAtom = 0;
        std::for_each(
            neighborTypes[atom].cbegin(),
            neighborTypes[atom].cend(),
            [&](const Int& neighborType)
            {
                const Real* forcePreContractPtr =
                    &forcePreContract[atom][basisSetSize * neighborType];
                // x component
                const Real* shecPtr = &(*clnmPairDerivative)[atom][basisSetSize * neighborAtom * 3];
                Real        dotProduct =
                    linalg::dotProduct(shecPtr, forcePreContractPtr, basisSetSize, linalgContext);
                pairForces[atom][counter] = -dotProduct;
                counter++;
                // y component
                shecPtr =
                    &(*clnmPairDerivative)[atom][basisSetSize * neighborAtom * 3 + basisSetSize];
                dotProduct =
                    linalg::dotProduct(shecPtr, forcePreContractPtr, basisSetSize, linalgContext);
                pairForces[atom][counter] = -dotProduct;
                counter++;
                // z component
                shecPtr = &(
                    *clnmPairDerivative)[atom][basisSetSize * neighborAtom * 3 + 2 * basisSetSize];
                dotProduct =
                    linalg::dotProduct(shecPtr, forcePreContractPtr, basisSetSize, linalgContext);
                pairForces[atom][counter] = -dotProduct;
                counter++;
                neighborAtom++;
            });
    }
}

void DescriptorSHS2::computeForceTermsGPU(const Vec2Real& forcePreContract,
                                          Vec2Real&       pairForces,
                                          Vec1Real&       centralForces) const
{

    std::vector<std::size_t> structureTypes(neighborList->get_numberTypes());
    std::iota(structureTypes.begin(), structureTypes.end(), 0);

    const Vec2Int& neighborTypes = neighborList->get_typeIndex();

    std::for_each(
        VASPML_SEQ
        centralAtomIndex.cbegin(),
        centralAtomIndex.cbegin() + centralAtomIndexSize.actDim,
        [&](const std::size_t atom)
        {
            std::for_each(
                structureTypes.begin(),
                structureTypes.end(),
                [&](const std::size_t neighborType)
                {
                    // x-component
                    const Real* forcePreContractPtr =
                        &forcePreContract[atom][basisSetSize * neighborType];
                    const Real* shecPtr =
                        &(*clnmDerivativeCentralVasp)[atom][basisSetSize * neighborType * 3];
                    centralForces[3 * atom] -= linalg::dotProduct(shecPtr,
                                                                  forcePreContractPtr,
                                                                  basisSetSize,
                                                                  linalgContext);
                    // y-component
                    forcePreContractPtr = &forcePreContract[atom][basisSetSize * neighborType];
                    shecPtr = &(*clnmDerivativeCentralVasp)[atom][basisSetSize * neighborType * 3
                                                                  + basisSetSize];
                    centralForces[3 * atom + 1] -= linalg::dotProduct(shecPtr,
                                                                      forcePreContractPtr,
                                                                      basisSetSize,
                                                                      linalgContext);
                    // z-component
                    forcePreContractPtr = &forcePreContract[atom][basisSetSize * neighborType];
                    shecPtr = &(*clnmDerivativeCentralVasp)[atom][basisSetSize * neighborType * 3
                                                                  + 2 * basisSetSize];
                    centralForces[3 * atom + 2] -= linalg::dotProduct(shecPtr,
                                                                      forcePreContractPtr,
                                                                      basisSetSize,
                                                                      linalgContext);
                });
            std::size_t counter = 0;
            std::size_t neighborAtom = 0;
            std::for_each(
                neighborTypes[atom].cbegin(),
                neighborTypes[atom].cend(),
                [&](const Int& neighborType)
                {
                    const Real* forcePreContractPtr =
                        &forcePreContract[atom][basisSetSize * neighborType];
                    // x component
                    const Real* shecPtr =
                        &(*clnmPairDerivative)[atom][basisSetSize * neighborAtom * 3];
                    Real dotProduct = linalg::dotProduct(shecPtr,
                                                         forcePreContractPtr,
                                                         basisSetSize,
                                                         linalgContext);
                    pairForces[atom][counter] = -dotProduct;
                    counter++;
                    // y component
                    shecPtr = &(
                        *clnmPairDerivative)[atom][basisSetSize * neighborAtom * 3 + basisSetSize];
                    dotProduct = linalg::dotProduct(shecPtr,
                                                    forcePreContractPtr,
                                                    basisSetSize,
                                                    linalgContext);
                    pairForces[atom][counter] = -dotProduct;
                    counter++;
                    // z component
                    shecPtr = &(*clnmPairDerivative)[atom][basisSetSize * neighborAtom * 3
                                                           + 2 * basisSetSize];
                    dotProduct = linalg::dotProduct(shecPtr,
                                                    forcePreContractPtr,
                                                    basisSetSize,
                                                    linalgContext);
                    pairForces[atom][counter] = -dotProduct;
                    counter++;
                    neighborAtom++;
                });
        });
}

void DescriptorSHS2::write_clnmVasp(const String& fname) const
{
    auto file = file_io::openFileO(fname);
    for (const auto& x : *clnmVasp)
        for (const auto& y : x) file << str("%24.16E ", y) << std::endl;
    file.close();
}
