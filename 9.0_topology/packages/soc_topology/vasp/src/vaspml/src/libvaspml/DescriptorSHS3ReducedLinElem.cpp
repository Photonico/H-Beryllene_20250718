#include "DescriptorSHS3ReducedLinElem.hpp"
#include "constants.hpp"
#include "nearest_neighbor.hpp"
#include "utils.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>

// #include <iostream>
// #include <fstream>

using namespace vaspml;

DescriptorSHS3ReducedLinElem::DescriptorSHS3ReducedLinElem(
    const Vec2Int&        descriptorList,
    const bool            angularFilterOn,
    const Int             angularFilterType,
    const Real            angularFilterScaling,
    const Vec1Size_t&     nRoots,
    const Int             maxOrder,
    const Real            weight,
    const bool            isNormalized,
    const ShVec2Real&     clnmVaspNoElement,
    const ShVec2Real&     spectrumVasp,
    const DescriptorType  descriptorType,
    const ExecutionPolicy algoExecution) :
    Descriptor(weight, isNormalized, descriptorType, algoExecution)
{

    if (clnmVaspNoElement == nullptr) { this->clnmVaspNoElement = std::make_shared<Vec2Real>(); }
    else { this->clnmVaspNoElement = clnmVaspNoElement; }
    if (spectrumVasp == nullptr) { this->spectrumVasp = std::make_shared<Vec2Real>(); }
    else { this->spectrumVasp = spectrumVasp; }

    this->angularFilterOn = angularFilterOn;
    this->angularFilterType = angularFilterType;
    this->angularFilterScaling = angularFilterScaling;
    descriptorSHS2 = nullptr;
    isSparseforcePreContractReady = false;
    numberElementsStructure = 0;
    numberElementsMLFF = descriptorList.size();
    make_angularFilter(maxOrder);
    make_sparseList(descriptorList, nRoots, maxOrder);
}

void DescriptorSHS3ReducedLinElem::computeSHS3(const std::shared_ptr<DescriptorSHS2>& shec,
                                               const TypeMap&                         typeMap)
{

    descriptorSHS2 = shec;
    numberAtoms = descriptorSHS2->get_nAtoms();
    set_neighborList(shec->get_neighborList_ptr());
    numberElementsStructure = neighborList->get_numberTypes();
    allocateArrays(typeMap);
    computeSHS3clnmNoElement();
    // setup sparse maps for derivatives
    if (!isSparseforcePreContractReady) prepareSparseforcePreContract(typeMap);
    // actual calculations
    switch (algoExecution)
    {
    case ExecutionPolicy::cpuSingleCore:
        computeSHS3CPU(typeMap);
        break;
    case ExecutionPolicy::gpuStdLib:
        VASPML_PARALLEL(
            computeSHS3GPU(typeMap);
        );
        break;
    }
}

void DescriptorSHS3ReducedLinElem::computeSHS3CPU(const TypeMap& typeMap)
{
    for (std::size_t atom = 0; atom < (std::size_t)numberAtoms; atom++)
    {
        computeSHS3SingleAtom(typeMap, atom);
    }
}

void DescriptorSHS3ReducedLinElem::computeSHS3GPU(const TypeMap& typeMap)
{
    std::vector<std::size_t> centralAtomIndex(numberAtoms);
    std::iota(centralAtomIndex.begin(), centralAtomIndex.end(), 0);

    std::for_each( // par_unseq,
        centralAtomIndex.cbegin(),
        centralAtomIndex.cbegin() + numberAtoms,
        [&](const std::size_t& atom) { computeSHS3SingleAtom(typeMap, atom); });
}

// Give back descriptor vector for a selected atom
const Vec1Real& DescriptorSHS3ReducedLinElem::get_SHS3Atom_vasp(const Int centralAtom) const
{
    return (*spectrumVasp)[centralAtom];
}

// compute element agnostic descriptors
void DescriptorSHS3ReducedLinElem::computeSHS3clnmNoElement(void)
{

    const Int basisSetSize = (Size_t)descriptorSHS2->get_basisSetSize();

    Vec2Real& clnmVaspNoElement = (*this->clnmVaspNoElement);

    // One dimension is number of atoms,
    // other is descriptor without the types
    for (std::size_t atom = 0; atom < (Size_t)numberAtoms; atom++)
    {
        // Get whole descriptor (l, n, m) for each atom
        const Vec1Real& cnlmVasp = descriptorSHS2->get_descriptor(atom);
        for (std::size_t type = 0; type < numberElementsStructure; type++)
        {
            linalg::scaleVectorPlusVector((Real)1.0,
                                          &cnlmVasp[type * basisSetSize],
                                          clnmVaspNoElement[atom],
                                          basisSetSize,
                                          linalgContext);
        }
    }
}

// get index of a given element agnostic Cnlm
std::size_t DescriptorSHS3ReducedLinElem::getIndexNoElement(const std::size_t l,
                                                            const std::size_t n,
                                                            const std::size_t m) const
{
    return descriptorSHS2->get_lOffset(l) + n * (2 * l + 1) + m;
}

// get array of summed clnm for a given central atom
const Vec1Real& DescriptorSHS3ReducedLinElem::get_clnmVaspNoElement(
    const std::size_t central_atom) const
{
    return (*clnmVaspNoElement)[central_atom];
}

// calculate actual reduced descriptor with respect to element type
void DescriptorSHS3ReducedLinElem::computeSHS3SingleAtom(const TypeMap&    typeMap,
                                                         const std::size_t atom)
{

    const Int centralType = neighborList->get_typeIndexCentral(atom);

    const Int centralType_ff = typeMap.toType(centralType);

    const Vec1Real& clnmVasp = descriptorSHS2->get_clnmVasp(atom);

    const Vec1Real& clnmVasp_noelement = get_clnmVaspNoElement(atom);

    std::size_t nDesc = 0;
    std::for_each(type1List[centralType_ff].cbegin(),
                  type1List[centralType_ff].cend(),
                  [&](const Int& type1_ff) mutable
                  {
                      const Int type1 = typeMap.toSubType(type1_ff);
                      if (type1 >= 0)
                      {
                          const Int   order = angularList[centralType_ff][nDesc];
                          const Int   n0 = n0List[centralType_ff][nDesc];
                          const Int   n1 = n1List[centralType_ff][nDesc];
                          Real        descTemp = 0;
                          std::size_t indx0 =
                              descriptorSHS2->get_Index(type1, order, n0, (std::size_t)0);
                          std::size_t indx1 = getIndexNoElement(order, n1, (Size_t)0);
                          for (std::size_t m = 0; m < (std::size_t)2 * order + 1; m++)
                          {
                              descTemp += clnmVasp[indx0] * clnmVasp_noelement[indx1];
                              indx0++;
                              indx1++;
                          }
                          (*spectrumVasp)[atom][nDesc] =
                              weightFactor[centralType_ff][nDesc] * descTemp;
                      }
                      nDesc++;
                  });
}

void DescriptorSHS3ReducedLinElem::make_sparseList(const Vec2Int&                  descriptorList,
                                                   const std::vector<std::size_t>& nRoots,
                                                   const Int                       maxOrder)
{

    type0List.resize(descriptorList.size());
    type1List.resize(descriptorList.size());
    angularList.resize(descriptorList.size());
    n0List.resize(descriptorList.size());
    n1List.resize(descriptorList.size());
    weightFactor.resize(descriptorList.size());

    for (std::size_t type0 = 0; type0 < descriptorList.size(); type0++)
    {
        Int descCounter = 0;
        for (std::size_t type1 = 0; type1 < descriptorList.size(); type1++)
        {
            // for ( std::size_t type2 = 0; type2 < descriptorList.size(); type2++ ){
            for (std::size_t orderL = 0; orderL < (std::size_t)maxOrder + 1; orderL++)
            {
                for (std::size_t order0 = 0; order0 < nRoots[orderL]; order0++)
                {
                    for (std::size_t order1 = order0; order1 < nRoots[orderL]; order1++)
                    {
                        descCounter++;
                        // check if descriptor list is sparsified
                        for (std::size_t nDesc = 0; nDesc < descriptorList[type0].size(); nDesc++)
                        {
                            if (descCounter == descriptorList[type0][nDesc])
                            {
                                // Frage: Ist type0List wirklich benötigt?
                                type0List[type0].push_back(type0);
                                type1List[type0].push_back(type1);
                                // type2List[ type0 ].push_back( type2 );
                                angularList[type0].push_back(orderL);
                                n0List[type0].push_back(order0);
                                n1List[type0].push_back(order1);
                                if (order0 == order1)
                                {
                                    weightFactor[type0].push_back((Real)1.0
                                                                  * angularFilter[orderL]);
                                }
                                else
                                {
                                    weightFactor[type0].push_back(constants::SQRT2
                                                                  * angularFilter[orderL]);
                                }
                                break;
                            }
                        }
                    }
                }
            }
            // }
        }
    }
}

void DescriptorSHS3ReducedLinElem::make_angularFilter(const Int maxOrder)
{
    angularFilter.resize(maxOrder + 1);
    Real factor = constants::PI * constants::PI * (Real)8.0;

    for (std::size_t order = 0; order < (std::size_t)maxOrder + 1; order++)
    {
        if (angularFilterOn)
        {
            if (angularFilterType == 1)
            {
                angularFilter[order] =
                    std::sqrt(factor / (Real)(2 * order + 1)) / std::sqrt((Real)(2 * order + 1));
            }
            else if (angularFilterType == 2)
            {
                Real factor2 =
                    angularFilterScaling * (Real)((order * (order + 1)) * (order * (order + 1)));
                factor2 = ((Real)1 + factor2);
                factor2 *= factor2;
                angularFilter[order] = std::sqrt(factor / ((Real)(2 * order + 1))) / factor2;
            }
        }
        else { angularFilter[order] = std::sqrt(factor / (Real)(2 * order + 1)); }
    }
}

void DescriptorSHS3ReducedLinElem::allocateArrays(const TypeMap& typeMap)
{

    auto& nList = *neighborList;
    spectrumVasp->resize(numberAtoms);
    if (algoExecution == ExecutionPolicy::cpuSingleCore)
    {
        for (std::size_t atom = 0; atom < (Size_t)numberAtoms; atom++)
        {
            const Int& type = nList.get_typeIndexCentral(atom);
            const Int  ff_type = typeMap.toType(type);
            (*spectrumVasp)[atom].resize(type0List[ff_type].size());
        }
        // allocate and fill with zeros
        clnmVaspNoElement->resize(numberAtoms);
        std::for_each(clnmVaspNoElement->begin(),
                      clnmVaspNoElement->end(),
                      [&](Vec1Real& slice)
                      {
                          slice.resize(descriptorSHS2->get_basisSetSize());
                          std::fill(slice.begin(), slice.end(), (Real)0.0);
                      });
    }
    else if (algoExecution == ExecutionPolicy::gpuStdLib) { resizeArraysGPU(typeMap); }
}

void DescriptorSHS3ReducedLinElem::resizeArraysGPU(const TypeMap& typeMap)
{
    bool resize = centralAtomIndexSize.checkResize(numberAtoms);
    if (resize)
    {
        auto nList = descriptorSHS2->get_neighborList();
        centralAtomIndex.resize(numberAtoms);
        std::iota(centralAtomIndex.begin(), centralAtomIndex.end(), 0);
        spectrumVaspSize.act1Dim = numberAtoms;
        spectrumVaspSize.max1Dim = numberAtoms;
        Vec1Int sizeArray = nList.get_typeIndexCentral();
        std::for_each(sizeArray.begin(),
                      sizeArray.end(),
                      [&](Int& size) { size = type0List[typeMap.toType(size)].size(); });
        spectrumVaspSize.resizeArray1Dim(*spectrumVasp, sizeArray);
        typeMapLoc = typeMap;
    }
    else
    {
        auto    nList = descriptorSHS2->get_neighborList();
        Vec1Int sizeArray = nList.get_typeIndexCentral();
        std::for_each(sizeArray.begin(),
                      sizeArray.end(),
                      [&](Int& size) { size = descriptorSize[typeMap.toType(size)]; });
        resize = spectrumVaspSize.checkResize2Dim(sizeArray);
        if (resize) { spectrumVaspSize.resizeArray2Dim(*spectrumVasp, sizeArray); }
    }
}

VASPML_NV_HOST_DEVICE
const Vec1Real& DescriptorSHS3ReducedLinElem::get_descriptor(const std::size_t centralAtom) const
{
    return (*spectrumVasp)[centralAtom];
}

VASPML_NV_HOST_DEVICE
Int DescriptorSHS3ReducedLinElem::get_sizeDescriptor(const std::size_t centralAtom) const
{
    return (*spectrumVasp)[centralAtom].size();
}

void DescriptorSHS3ReducedLinElem::rescale_descriptor(const std::size_t centralAtom,
                                                      const Real        scaleFactor)
{
    std::transform( // par_unseq,
        (*spectrumVasp)[centralAtom].begin(),
        (*spectrumVasp)[centralAtom].end(),
        (*spectrumVasp)[centralAtom].begin(),
        std::bind(std::multiplies<Real>(), std::placeholders::_1, scaleFactor));
}

void DescriptorSHS3ReducedLinElem::compute_forcePreContract(const Vec2Real& derivativeMatrix,
                                                            Vec2Real&       forcePreContract,
                                                            const TypeMap&  typeMap) const
{
    VASPML_DEBUG_L2(
        if (!isSparseforcePreContractReady)
        {
            throw std::runtime_error("ERROR: You are trying to use compute_forcePreContract(const "
                                     "ShVec2Real & derivativeMatrix, ShVec2Real& forcePreContract, "
                                     "const TypeMap& typeMap) const.\n"
                                     "Sparse index arrays were not computed yet. Make sure that "
                                     "prepareSparseforcePreContract( const TypeMap& typeMap ) is "
                                     "called before");
        }
    );
    VASPML_PROFILING_START("DescriptorSHS3ReducedLinElem::compute_forcePreContract");
    switch (algoExecution)
    {
    case ExecutionPolicy::cpuSingleCore:
        compute_forcePreContractSparseCPU(derivativeMatrix, forcePreContract, typeMap);
        break;
    case ExecutionPolicy::gpuStdLib:
        VASPML_PARALLEL(
            compute_forcePreContractSparseGPU(derivativeMatrix, forcePreContract, typeMap);
        );
        break;
    }
    VASPML_PROFILING_STOP("DescriptorSHS3ReducedLinElem::compute_forcePreContract");
}

void DescriptorSHS3ReducedLinElem::compute_forcePreContractSparseCPU(
    const Vec2Real& derivativeMatrix,
    Vec2Real&       forcePreContract,
    const TypeMap&  typeMap) const
{

    const Vec1Int& typeIndexCentral = neighborList->get_typeIndexCentral();
    const Vec1Int& centralAtomIndexPerType = neighborList->get_centralAtomIndexPerType();
    std::size_t    derivativeAtom = 0;
    std::for_each(
        forcePreContract.begin(),
        forcePreContract.begin() + numberAtoms,
        [&](Vec1Real& slice)
        {
            std::size_t     typeStruc = typeIndexCentral[derivativeAtom];
            std::size_t     typeStrucFF = typeMap.toType(typeStruc);
            std::size_t     centralAtomShift = n0List[typeStrucFF].size();
            const Vec1Real& descriptor2Body = descriptorSHS2->get_descriptor(derivativeAtom);
            std::size_t     centralAtomPerType = centralAtomIndexPerType[derivativeAtom];
            std::size_t     centralProductTmp = centralAtomShift * centralAtomPerType;

            // first contraction
            std::size_t nCombi = 0;
            std::for_each(
                preFactorPreContractB[typeStrucFF].cbegin(),
                preFactorPreContractB[typeStrucFF].cend(),
                [&](const Real& preFac)
                {
                    // calculate sum over all element types for Cnlms
                    Real presumTemp = 0;
                    for (std::size_t neighborType2 = 0; neighborType2 < numberElementsStructure;
                         neighborType2++)
                    {
                        presumTemp +=
                            descriptor2Body[sparseMap_SHS2_B[typeStrucFF][neighborType2][nCombi]];
                    }
                    std::size_t index_derivativeMatrix =
                        centralProductTmp + sparseMap_derivativeMatrix_B[typeStrucFF][nCombi];
                    std::size_t centralIndex = sparseMap_central_B[typeStrucFF][nCombi];
                    slice[centralIndex] +=
                        preFac * derivativeMatrix[typeStruc][index_derivativeMatrix] * presumTemp;
                    nCombi++;
                });

            // second contraction
            nCombi = 0;
            std::for_each(preFactorPreContract[typeStrucFF].cbegin(),
                          preFactorPreContract[typeStrucFF].cend(),
                          [&](const Real& preFac)
                          {
                              std::size_t index_derivativeMatrix =
                                  centralProductTmp
                                  + sparseMap_derivativeMatrix[typeStrucFF][nCombi];
                              std::size_t centralIndex = sparseMap_central[typeStrucFF][nCombi];
                              std::size_t shs2Index = sparseMap_SHS2[typeStrucFF][nCombi];
                              slice[centralIndex] +=
                                  preFac * derivativeMatrix[typeStruc][index_derivativeMatrix]
                                  * descriptor2Body[shs2Index];
                              nCombi++;
                          });

            derivativeAtom++;
        });
}

void DescriptorSHS3ReducedLinElem::compute_forcePreContractSparseGPU(
    const Vec2Real& derivativeMatrix,
    Vec2Real&       forcePreContract,
    const TypeMap&  typeMap) const
{

    const Vec1Int& typeIndexCentral = neighborList->get_typeIndexCentral();
    const Vec1Int& centralAtomIndexPerType = neighborList->get_centralAtomIndexPerType();

    std::vector<std::size_t> centralAtomIndex(numberAtoms);
    std::iota(centralAtomIndex.begin(), centralAtomIndex.end(), 0);

    std::for_each( // seq,
        centralAtomIndex.cbegin(),
        centralAtomIndex.cbegin() + numberAtoms,
        [this,
         derivativeMatrix,
         forcePreContract,
         typeMap,
         typeIndexCentral,
         centralAtomIndexPerType](const std::size_t& derivativeAtom) mutable
        {
            std::size_t     typeStruc = typeIndexCentral[derivativeAtom];
            std::size_t     typeStrucFF = typeMap.toType(typeStruc);
            std::size_t     centralAtomShift = n0List[typeStrucFF].size();
            const Vec1Real& descriptor2Body = descriptorSHS2->get_descriptor(derivativeAtom);
            std::size_t     centralAtomPerType = centralAtomIndexPerType[derivativeAtom];
            std::size_t     centralProductTmp = centralAtomShift * centralAtomPerType;
            // first contraction
            std::size_t nCombi = 0;
            std::for_each(
                preFactorPreContractB[typeStrucFF].cbegin(),
                preFactorPreContractB[typeStrucFF].cend(),
                [&](const Real& preFac)
                {
                    // calculate sum over all element types for Cnlms
                    Real presumTemp = 0;
                    for (std::size_t neighborType2 = 0; neighborType2 < numberElementsStructure;
                         neighborType2++)
                    {
                        presumTemp +=
                            descriptor2Body[sparseMap_SHS2_B[typeStrucFF][neighborType2][nCombi]];
                    }
                    std::size_t index_derivativeMatrix =
                        centralProductTmp + sparseMap_derivativeMatrix_B[typeStrucFF][nCombi];
                    std::size_t centralIndex = sparseMap_central_B[typeStrucFF][nCombi];
                    forcePreContract[derivativeAtom][centralIndex] +=
                        preFac * derivativeMatrix[typeStruc][index_derivativeMatrix] * presumTemp;
                    nCombi++;
                });
            // second contraction
            nCombi = 0;
            std::for_each(preFactorPreContract[typeStrucFF].cbegin(),
                          preFactorPreContract[typeStrucFF].cend(),
                          [&](const Real& preFac)
                          {
                              std::size_t index_derivativeMatrix =
                                  centralProductTmp
                                  + sparseMap_derivativeMatrix[typeStrucFF][nCombi];
                              std::size_t centralIndex = sparseMap_central[typeStrucFF][nCombi];
                              std::size_t shs2Index = sparseMap_SHS2[typeStrucFF][nCombi];
                              forcePreContract[derivativeAtom][centralIndex] +=
                                  preFac * derivativeMatrix[typeStruc][index_derivativeMatrix]
                                  * descriptor2Body[shs2Index];
                              nCombi++;
                          });
        });
}

void DescriptorSHS3ReducedLinElem::prepareSparseforcePreContract(const TypeMap& typeMap)
{
    sparseMap_derivativeMatrix.resize(n0List.size());
    preFactorPreContract.resize(n0List.size());
    sparseMap_SHS2.resize(n0List.size());
    sparseMap_central.resize(n0List.size());

    sparseMap_derivativeMatrix_B.resize(n0List.size());
    preFactorPreContractB.resize(n0List.size());
    vector_tools::allocate_vector(sparseMap_SHS2_B, n0List.size(), numberElementsStructure);
    sparseMap_central_B.resize(n0List.size());

    // loop over central atom type
    for (std::size_t centralType = 0; centralType < n0List.size(); centralType++)
    {
        // loop over central atom neighbor types for the derivative Cnlm
        for (std::size_t centralNeighborTypeFF = 0; centralNeighborTypeFF < n0List.size();
             centralNeighborTypeFF++)
        {
            Int centralNeighborType = typeMap.toSubType(centralNeighborTypeFF);
            if (centralNeighborType < 0) continue;
            // loop over the l index of derivative Cnlm.
            for (std::size_t l0 = 0; l0 < descriptorSHS2->get_maxOrder() + 1; l0++)
            {
                // loop over the the radial index of the derivative Cnlm
                for (std::size_t nRadial0 = 0; nRadial0 < descriptorSHS2->get_nRoots()[l0];
                     nRadial0++)
                {
                    // loop over the m index of the clnm of the central atom.
                    // One has to filter prefactors and index maps with respect to all possible
                    // combinations of the central atom clnm
                    for (std::size_t m0 = 0; m0 < 2 * l0 + 1; m0++)
                    {
                        std::size_t nDesc0 =
                            descriptorSHS2->get_Index((std::size_t)centralNeighborType,
                                                      l0,
                                                      nRadial0,
                                                      m0);
                        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        // loop over all possible descriptor combinations p^{iJ_{1}}_{ln_{1}n_{2}}
                        for (std::size_t nDesc1 = 0; nDesc1 < n0List[centralType].size(); nDesc1++)
                        {
                            if (l0 != angularList[centralType][nDesc1]) continue;
                            std::size_t neighborType1FF = type1List[centralType][nDesc1];
                            Int         neighborType1 = typeMap.toSubType(neighborType1FF);
                            if (neighborType1 < 0) continue;
                            std::size_t nRadial1 = n0List[centralType][nDesc1];
                            std::size_t nRadial2 = n1List[centralType][nDesc1];
                            std::size_t nDesc2 =
                                descriptorSHS2->get_Index((std::size_t)neighborType1,
                                                          l0,
                                                          nRadial1,
                                                          m0);
                            if (nDesc0 == nDesc2)
                            {
                                sparseMap_central_B[centralType].push_back(nDesc0);
                                sparseMap_derivativeMatrix_B[centralType].push_back(nDesc1);
                                if (nRadial1 == nRadial2)
                                {
                                    preFactorPreContractB[centralType].push_back(angularFilter[l0]);
                                }
                                else
                                {
                                    preFactorPreContractB[centralType].push_back(
                                        std::sqrt((Real)2.0) * angularFilter[l0]);
                                }
                                for (std::size_t neighborType2 = 0;
                                     neighborType2 < numberElementsStructure;
                                     neighborType2++)
                                {
                                    std::size_t nDesc3 =
                                        descriptorSHS2->get_Index((std::size_t)neighborType2,
                                                                  l0,
                                                                  nRadial2,
                                                                  m0);
                                    sparseMap_SHS2_B[centralType][neighborType2].push_back(nDesc3);
                                }
                            }
                            for (std::size_t neighborType2 = 0;
                                 neighborType2 < numberElementsStructure;
                                 neighborType2++)
                            {
                                std::size_t nDesc3 =
                                    descriptorSHS2->get_Index((std::size_t)neighborType2,
                                                              l0,
                                                              nRadial2,
                                                              m0);
                                if (nDesc0 == nDesc3)
                                {
                                    sparseMap_central[centralType].push_back(nDesc0);
                                    sparseMap_derivativeMatrix[centralType].push_back(nDesc1);
                                    if (nRadial1 == nRadial2)
                                    {
                                        preFactorPreContract[centralType].push_back(
                                            angularFilter[l0]);
                                    }
                                    else
                                    {
                                        preFactorPreContract[centralType].push_back(
                                            std::sqrt((Real)2.0) * angularFilter[l0]);
                                    }
                                    sparseMap_SHS2[centralType].push_back(nDesc2);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    isSparseforcePreContractReady = true;
}

std::size_t DescriptorSHS3ReducedLinElem::get_forcePreContractSize(
    const std::size_t atomIndex) const
{
    return descriptorSHS2->get_forcePreContractSize(atomIndex);
}

void DescriptorSHS3ReducedLinElem::computeForceTerms(const Vec2Real& forcePreContract,
                                                     Vec2Real&       pairForces,
                                                     Vec1Real&       centralForces) const
{
    descriptorSHS2->computeForceTerms(forcePreContract, pairForces, centralForces);
}
