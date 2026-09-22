#include "DescriptorSHS3.hpp"
#include "Timer.hpp"
#include "constants.hpp"
#include "debug.hpp"
#include "math.hpp"
#include "nearest_neighbor.hpp"
#include "utils.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>

using namespace vaspml;

DescriptorSHS3::DescriptorSHS3(const Vec2Int&                  descriptorList,
                               const bool                      angularFilterOn,
                               const Int                       angularFilterType,
                               const Real                      angularFilterScaling,
                               const std::vector<std::size_t>& nRoots,
                               const Int                       maxOrder,
                               const Real                      weight,
                               const bool                      isNormalized,
                               const ShVec2Real&               spectrumVasp,
                               const ShVec2Real&               spectrum_pair,
                               const DescriptorType            descriptorType,
                               const ExecutionPolicy           algoExecution) :
    Descriptor(weight, isNormalized, descriptorType, algoExecution)
{
    if (spectrumVasp == nullptr) { this->spectrumVasp = std::make_shared<Vec2Real>(); }
    else { this->spectrumVasp = spectrumVasp; }

    if (spectrum_pair == nullptr) { this->spectrumPair = std::make_shared<Vec2Real>(); }
    else { this->spectrumPair = spectrum_pair; }

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

void DescriptorSHS3::computeSHS3(const std::shared_ptr<DescriptorSHS2>& descriptorSHS2,
                                 const TypeMap&                         typeMap)
{

    this->descriptorSHS2 = descriptorSHS2;
    set_neighborList(this->descriptorSHS2->get_neighborList_ptr());
    numberAtoms = descriptorSHS2->get_nAtoms();
    numberElementsStructure = typeMap.countStructureTypes();
    allocateArrays(typeMap);

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

void DescriptorSHS3::computeSHS3CPU(const TypeMap& /* typeMap */)
{
    for (std::size_t atom = 0; atom < (std::size_t)numberAtoms; atom++)
    {
        computeSHS3SingleAtom(*descriptorSHS2, atom);
    }
}

void DescriptorSHS3::computeSHS3GPU(const TypeMap& /* typeMap */)
{
    std::for_each(VASPML_PAR_UNSEQ
                  centralAtomIndex.cbegin(),
                  centralAtomIndex.cbegin() + centralAtomIndexSize.actDim,
                  [&](const std::size_t& atom) { computeSHS3SingleAtom(*descriptorSHS2, atom); });
}

void DescriptorSHS3::computeSHS3SingleAtom(const DescriptorSHS2& descriptorSHS2,
                                           const std::size_t     atom)
{

    const Int centralType = neighborList->get_typeIndexCentral(atom);

    const Int       centralType_ff = typeMapLoc.toType(centralType);
    const Vec1Real& clnm_vasp = descriptorSHS2.get_clnmVasp(atom);
    Vec1Real&       spectrumVasp = (*this->spectrumVasp)[atom];

    std::size_t nDesc = 0;
    std::for_each(type1List[centralType_ff].cbegin(),
                  type1List[centralType_ff].cend(),
                  [&](const Int& type1_ff) mutable
                  {
                      const Int type2_ff = type2List[centralType_ff][nDesc];
                      const Int type1 = typeMapLoc.toSubType(type1_ff);
                      const Int type2 = typeMapLoc.toSubType(type2_ff);
                      if (type1 >= 0 && type2 >= 0)
                      {
                          const Int order = angularList[centralType_ff][nDesc];
                          const Int n0 = n0List[centralType_ff][nDesc];
                          const Int n1 = n1List[centralType_ff][nDesc];
                          Real      descTemp = 0;
                          Size_t    indx0 = descriptorSHS2.get_Index(type1, order, n0, (Size_t)0);
                          Size_t    indx1 = descriptorSHS2.get_Index(type2, order, n1, (Size_t)0);
                          for (std::size_t m = 0; m < (Size_t)2 * order + 1; m++)
                          {
                              descTemp += clnm_vasp[indx0] * clnm_vasp[indx1];
                              indx0++;
                              indx1++;
                          }
                          spectrumVasp[nDesc] = weightFactor[centralType_ff][nDesc] * descTemp;
                      }
                      nDesc++;
                  });
}

void DescriptorSHS3::make_sparseList(const Vec2Int&                  descriptorList,
                                     const std::vector<std::size_t>& nRoots,
                                     const Int                       maxOrder)
{

    type0List.resize(descriptorList.size());
    type1List.resize(descriptorList.size());
    type2List.resize(descriptorList.size());
    angularList.resize(descriptorList.size());
    n0List.resize(descriptorList.size());
    n1List.resize(descriptorList.size());
    weightFactor.resize(descriptorList.size());

    for (std::size_t type0 = 0; type0 < descriptorList.size(); type0++)
    {
        Int descCounter = 0;
        for (std::size_t type1 = 0; type1 < descriptorList.size(); type1++)
        {
            for (std::size_t type2 = 0; type2 < descriptorList.size(); type2++)
            {
                for (std::size_t orderL = 0; orderL < (std::size_t)maxOrder + 1; orderL++)
                {
                    for (std::size_t order0 = 0; order0 < nRoots[orderL]; order0++)
                    {
                        for (std::size_t order1 = order0; order1 < nRoots[orderL]; order1++)
                        {
                            descCounter++;
                            // check if descriptor list is sparsified
                            for (std::size_t nDesc = 0; nDesc < descriptorList[type0].size();
                                 nDesc++)
                            {
                                if (descCounter == descriptorList[type0][nDesc])
                                {
                                    type0List[type0].push_back(type0);
                                    type1List[type0].push_back(type1);
                                    type2List[type0].push_back(type2);
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
            }
        }
        n0ListSize.push_back(n0List[type0].size());
    }
}

void DescriptorSHS3::make_angularFilter(const Int maxOrder)
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

void DescriptorSHS3::allocateArrays(const TypeMap& typeMap)
{

    if (algoExecution == ExecutionPolicy::cpuSingleCore)
    {
        typeMapLoc = typeMap;
        auto nList = descriptorSHS2->get_neighborList();
        spectrumVasp->resize(numberAtoms);

        for (std::size_t atom = 0; atom < (std::size_t)numberAtoms; atom++)
        {
            const Int& type = nList.get_typeIndexCentral(atom);
            const Int  ff_type = typeMap.toType(type);
            (*spectrumVasp)[atom].resize(type0List[ff_type].size());
        }
    }
    else if (algoExecution == ExecutionPolicy::gpuStdLib) { resizeArraysGPU(typeMap); }
}

void DescriptorSHS3::resizeArraysGPU(const TypeMap& typeMap)
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
        spectrumVaspSize.act1Dim = numberAtoms;
        auto    nList = descriptorSHS2->get_neighborList();
        Vec1Int sizeArray = nList.get_typeIndexCentral();
        std::for_each(sizeArray.begin(),
                      sizeArray.end(),
                      [&](Int& size) { size = descriptorSize[typeMap.toType(size)]; });
        resize = spectrumVaspSize.checkResize2Dim(sizeArray);
        if (resize) { spectrumVaspSize.resizeArray2Dim(*spectrumVasp, sizeArray); }
    }
}

const Vec1Real& DescriptorSHS3::get_SHS3Atom_vasp(const Int centralAtom) const
{
    return (*spectrumVasp)[centralAtom];
}

VASPML_NV_HOST_DEVICE
const Vec1Real& DescriptorSHS3::get_descriptor(const std::size_t centralAtom) const
{
    return (*spectrumVasp)[centralAtom];
}

VASPML_NV_HOST_DEVICE
Int DescriptorSHS3::get_sizeDescriptor(const std::size_t centralAtom) const
{
    return (*spectrumVasp)[centralAtom].size();
}

void DescriptorSHS3::rescale_descriptor(const std::size_t centralAtom, const Real scaleFactor)
{
    std::transform( // par_unseq,
        (*spectrumVasp)[centralAtom].begin(),
        (*spectrumVasp)[centralAtom].end(),
        (*spectrumVasp)[centralAtom].begin(),
        std::bind(std::multiplies<Real>(), std::placeholders::_1, scaleFactor));
}

void DescriptorSHS3::compute_forcePreContract(const Vec2Real& derivativeMatrix,
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
    // #pragma acc data copy(this,derivativeMatrix,forcePreContract)
    //  {
    VASPML_PROFILING_START("DescriptorSHS3::compute_forcePreContract");
    //   switch ( algoExecution ){
    //      case ExecutionPolicy::cpuSingleCore:
    compute_forcePreContractSparseSingleCPU(derivativeMatrix, forcePreContract, typeMap);
    //         break;
    //      case ExecutionPolicy::gpuStdLib:
    //         VASPML_PARALLEL( compute_forcePreContractSparseGPU( derivativeMatrix,
    //                                                             forcePreContract,
    //                                                             typeMap );
    //         );
    //         break;
    //   }
    VASPML_PROFILING_STOP("DescriptorSHS3::compute_forcePreContract");
    // }
}

void DescriptorSHS3::compute_forcePreContractSparseSingleCPU(const Vec2Real& derivativeMatrix,
                                                             Vec2Real&       forcePreContract,
                                                             const TypeMap&  typeMap) const
{

    const Vec1Int& typeIndexCentral = neighborList->get_typeIndexCentral();
    const Vec1Int& centralAtomIndexPerType = neighborList->get_centralAtomIndexPerType();
    std::size_t    centralAtom = 0;
    std::for_each(
        forcePreContract.begin(),
        forcePreContract.begin() + numberAtoms,
        [&](Vec1Real& slice)
        {
            std::size_t     typeStruc = typeIndexCentral[centralAtom];
            std::size_t     typeStrucFF = typeMap.toType(typeStruc);
            std::size_t     centralAtomShift = n0List[typeStrucFF].size();
            const Vec1Real& descriptor2Body = descriptorSHS2->get_descriptor(centralAtom);
            std::size_t     centralAtomPerType = centralAtomIndexPerType[centralAtom];
            std::size_t     nCombi = 0;
            std::for_each(
                preFactorPreContract[typeStrucFF].cbegin(),
                preFactorPreContract[typeStrucFF].cend(),
                [&](const Real& preFac)
                {
                    std::size_t index_derivativeMatrix =
                        sparseMap_derivativeMatrix[typeStrucFF][nCombi];
                    std::size_t shs2Index = sparseMap_SHS2[typeStrucFF][nCombi];
                    std::size_t lnmCombiIndex = sparseMap_central[typeStrucFF][nCombi];
                    slice[lnmCombiIndex] +=
                        preFac
                        * derivativeMatrix[typeStruc][centralAtomShift * centralAtomPerType
                                                      + index_derivativeMatrix]
                        * descriptor2Body[shs2Index];
                    nCombi++;
                });

            centralAtom++;
        });
}

void DescriptorSHS3::compute_forcePreContractSparseGPU(const Vec2Real& derivativeMatrix,
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
         typeIndexCentral,
         centralAtomIndexPerType,
         derivativeMatrix,
         forcePreContract,
         typeMap](const std::size_t& centralAtom) mutable
        {
            std::size_t     typeStruc = typeIndexCentral[centralAtom];
            std::size_t     typeStrucFF = typeMap.toType(typeStruc);
            std::size_t     centralAtomShift = n0List[typeStrucFF].size();
            const Vec1Real& descriptor2Body = descriptorSHS2->get_descriptor(centralAtom);
            std::size_t     centralAtomPerType = centralAtomIndexPerType[centralAtom];
            std::size_t     nCombi = 0;
            std::for_each(
                preFactorPreContract[typeStrucFF].cbegin(),
                preFactorPreContract[typeStrucFF].cend(),
                [&](const Real& preFac)
                {
                    std::size_t index_derivativeMatrix =
                        sparseMap_derivativeMatrix[typeStrucFF][nCombi];
                    std::size_t shs2Index = sparseMap_SHS2[typeStrucFF][nCombi];
                    std::size_t lnmCombiIndex = sparseMap_central[typeStrucFF][nCombi];
                    forcePreContract[centralAtom][lnmCombiIndex] +=
                        preFac
                        * derivativeMatrix[typeStruc][centralAtomShift * centralAtomPerType
                                                      + index_derivativeMatrix]
                        * descriptor2Body[shs2Index];
                    nCombi++;
                });
        });
}

std::size_t DescriptorSHS3::get_forcePreContractSize(std::size_t atomIndex) const
{
    return descriptorSHS2->get_forcePreContractSize(atomIndex);
}

void DescriptorSHS3::computeForceTerms(const Vec2Real& forcePreContract,
                                       Vec2Real&       pairForces,
                                       Vec1Real&       centralForces) const
{
    descriptorSHS2->computeForceTerms(forcePreContract, pairForces, centralForces);
}

void DescriptorSHS3::prepareSparseforcePreContract(const TypeMap& typeMap)
{

    sparseMap_derivativeMatrix.resize(n0List.size());
    preFactorPreContract.resize(n0List.size());
    sparseMap_SHS2.resize(n0List.size());
    sparseMap_central.resize(n0List.size());
    // loop over central atom type
    for (std::size_t centralTypeFF = 0; centralTypeFF < n0List.size(); centralTypeFF++)
    {
        // loop over central atom neighbor types. These are neighbor types as in the clnm
        for (std::size_t centralNeighborTypeFF = 0; centralNeighborTypeFF < n0List.size();
             centralNeighborTypeFF++)
        {
            Int centralNeighborType = typeMap.toSubType(centralNeighborTypeFF);
            if (centralNeighborType < 0) continue;
            // loop over the l index of the central clnm.
            for (std::size_t l0 = 0; l0 < descriptorSHS2->get_maxOrder() + 1; l0++)
            {
                // loop over the the radial index of the clnm of the central atom
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
                        // loop over all possible descriptor combinations
                        // p^{iJ_{1}J_{2}}_{ln_{1}n_{2}}
                        for (std::size_t nDesc1 = 0; nDesc1 < n0List[centralTypeFF].size();
                             nDesc1++)
                        {
                            std::size_t neighborType1FF = type1List[centralTypeFF][nDesc1];
                            std::size_t neighborType2FF = type2List[centralTypeFF][nDesc1];
                            if (l0 != angularList[centralTypeFF][nDesc1]) continue;
                            Int neighborType1 = typeMap.toSubType(neighborType1FF);
                            Int neighborType2 = typeMap.toSubType(neighborType2FF);
                            if (neighborType1 < 0 || neighborType2 < 0) continue;
                            std::size_t nRadial1 = n0List[centralTypeFF][nDesc1];
                            std::size_t nRadial2 = n1List[centralTypeFF][nDesc1];
                            std::size_t nDesc2 =
                                descriptorSHS2->get_Index((std::size_t)neighborType1,
                                                          l0,
                                                          nRadial1,
                                                          m0);
                            std::size_t nDesc3 =
                                descriptorSHS2->get_Index((std::size_t)neighborType2,
                                                          l0,
                                                          nRadial2,
                                                          m0);
                            // check if the derivative
                            // \frac{p^{iJ_{1}J_{2}}_{ln_{1}n_{2}}}{c^{iJ_{3}}_{ln_{3}m}} is non
                            // zero.
                            if (nDesc0 == nDesc2 || nDesc0 == nDesc3)
                            {
                                sparseMap_central[centralTypeFF].push_back(nDesc0);
                                sparseMap_derivativeMatrix[centralTypeFF].push_back(nDesc1);
                                if (nRadial1 == nRadial2 && neighborType1 == neighborType2)
                                {
                                    preFactorPreContract[centralTypeFF].push_back(
                                        (Real)2.0 * angularFilter[l0]);
                                }
                                else
                                {
                                    if (nRadial1 == nRadial2 && neighborType1 != neighborType2)
                                    {
                                        preFactorPreContract[centralTypeFF].push_back(
                                            angularFilter[l0]);
                                    }
                                    else
                                    {
                                        preFactorPreContract[centralTypeFF].push_back(
                                            std::sqrt((Real)2.0) * angularFilter[l0]);
                                    }
                                }
                                // store the entries of the descriptorSHS2 which have to be
                                // multiplied by the derivativeMatrix and summed over. note always
                                // the entry which is not the one in target forcePreContract is
                                // stored.
                                if (nDesc0 == nDesc2)
                                {
                                    sparseMap_SHS2[centralTypeFF].push_back(nDesc3);
                                }
                                else { sparseMap_SHS2[centralTypeFF].push_back(nDesc2); }
                            }
                        }
                        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    }
                }
            }
        }
        descriptorSize.push_back(n0List[centralTypeFF].size());
    }
    isSparseforcePreContractReady = true;
}
