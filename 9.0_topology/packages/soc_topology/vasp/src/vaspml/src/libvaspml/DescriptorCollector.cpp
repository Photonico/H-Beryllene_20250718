#include "DescriptorCollector.hpp"
#include "DescriptorSHS2.hpp"
#include "ParallelEnvironemt.hpp"
#include "constants.hpp"
#include "debug.hpp"
#include "utils.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

using namespace vaspml;

template<>
SmartEnum<DescriptorStorage>::EnumMap SmartEnum<DescriptorStorage>::names = {
    {DescriptorStorage::Type,        "TypeOrder"       },
    {DescriptorStorage::CentralAtom, "CentralAtomOrder"},
};

DescriptorCollector::DescriptorCollector(Int                                 storageOrder,
                                         const std::map<String, ShVec2Real>& descriptorsNormalized,
                                         const ShVec1Real&                   normsAtom,
                                         const ExecutionPolicy               algoExecution) :
    algoExecution(algoExecution),
    descriptorKeyList(constants::descriptorKeyList)
{
    switch (storageOrder)
    {
    case 0:
        this->storageOrder = DescriptorStorage::Type;
        break;
    case 1:
        this->storageOrder = DescriptorStorage::CentralAtom;
        break;
    }
    if (descriptorsNormalized.empty())
    {
        this->descriptorsNormalized["SHS2-2-body"] = std::make_shared<Vec2Real>();
        this->descriptorsNormalized["SHS3-3-body"] = std::make_shared<Vec2Real>();
    }
    else { this->descriptorsNormalized = descriptorsNormalized; }
    if (normsAtom == nullptr) { this->normsAtom = std::make_shared<Vec1Real>(); }
    else { this->normsAtom = normsAtom; }
    for (const String& key : constants::descriptorKeyList)
    {
        length_descriptorsNormalized[key] = std::make_shared<Vec1Int>();
        descriptorsNormalizedSize[key] = ArrayResizing2D();
        length_descriptorsNormalizedSize[key] = ArrayResizing1D();
    }
}

void DescriptorCollector::setDescriptorMap(
    const std::map<String, std::shared_ptr<Descriptor>>& descriptors)
{
    this->descriptors = descriptors;
}

void DescriptorCollector::updateCollector(void)
{

    VASPML_DEBUG_L1(
        for (std::size_t i = 0; i < constants::descriptorKeyList.size()-1; i++)
        {
           for (std::size_t j = i+1; j < constants::descriptorKeyList.size(); j++)
           {
              if (descriptors.at(constants::descriptorKeyList[i])->get_nAtoms() != 
                  descriptors.at(constants::descriptorKeyList[j])->get_nAtoms())
              {
                global_scope::tutor.bug("ERROR: DescriptorCollector::updateCollector(): "
                                        "the number of atoms does not agree in your descriptors");
              }
           }
        }
    );
    numberAtoms = descriptors.at("SHS2-2-body")->get_nAtoms();
    allocateArrays();
    normalizeDescriptors();
}

void DescriptorCollector::normalizeDescriptors(void)
{

    switch (algoExecution)
    {
    case ExecutionPolicy::cpuSingleCore:
        normalizeDescriptorsCPU();
        break;
    case ExecutionPolicy::gpuStdLib:
        normalizeDescriptorsGPU();
        break;
    }
}

void DescriptorCollector::normalizeDescriptorsCPU(void)
{
    Vec1Real& normsAtom = (*this->normsAtom);
    for (std::size_t i = 0; i < (std::size_t)numberAtoms; i++)
    {
        Real norm = computeNormSingleAtom(i);
        if (storageOrder == DescriptorStorage::Type)
        {
            normalizeDescriptorsSingleAtomTypeOrder(i, norm);
        }
        else if (storageOrder == DescriptorStorage::CentralAtom)
        {
            normalizeDescriptorsSingleAtomCentralOrder(i, norm);
        }
        normsAtom[i] = norm;
    }
}

void DescriptorCollector::normalizeDescriptorsGPU(void)
{

    Vec1Real& normsAtom = (*this->normsAtom);
    std::for_each(VASPML_SEQ // TODO THIS SHOULD BE A VASPML_PAR_UNSEQ. Static linking for cubblas
                             // libraries needed
                  centralAtomIndex.cbegin(),
                  centralAtomIndex.cbegin() + centralAtomIndexSize.actDim,
                  [&](const std::size_t& atom) mutable
                  {
                      Real norm = computeNormSingleAtom(atom);
                      if (storageOrder == DescriptorStorage::Type)
                      {
                          normalizeDescriptorsSingleAtomTypeOrder(atom, norm);
                      }
                      else if (storageOrder == DescriptorStorage::CentralAtom)
                      {
                          normalizeDescriptorsSingleAtomCentralOrder(atom, norm);
                      }
                      normsAtom[atom] = norm;
                  });
}

Real DescriptorCollector::computeNormSingleAtom(std::size_t atomIndex)
{
    Real normTotal = 0;
    for (const String& key : descriptorKeyList)
    {
        Real weight = descriptors.at(key)->get_weight();
        if (weight > 0)
        {
            // const Vec1Real& data = descriptors.at( key ) -> get_descriptor( atomIndex );
            Real normDescriptor = linalg::l2Norm(descriptors.at(key)->get_descriptor(atomIndex),
                                                 descriptors.at(key)->get_sizeDescriptor(atomIndex),
                                                 linalgContext);
            normTotal += weight * normDescriptor * normDescriptor;
        }
    }
    return sqrt(normTotal);
}

void DescriptorCollector::normalizeDescriptorsSingleAtomTypeOrder(std::size_t atomIndex,
                                                                  const Real  norm)
{
    for (const String& key : descriptorKeyList)
    {
        if (descriptors.at(key)->get_isNormalized())
        {
            Real weight = descriptors.at(key)->get_weight();
            if (weight > (Real)0)
            {
                if (norm > constants::EPS_TOL)
                {
                    Real            factor = sqrt(weight) / norm;
                    const Vec1Real& desc = descriptors.at(key)->get_descriptor(atomIndex);
                    const Int&      type = descriptors.at(key)->get_typeIndexCentral(atomIndex);
                    Vec1Real& descriptorsNormalized = (*this->descriptorsNormalized[key])[type];
                    Int&      length_descriptorsNormalized =
                        (*this->length_descriptorsNormalized[key])[type];

                    std::transform( // par_unseq,
                        desc.cbegin(),
                        desc.cend(),
                        descriptorsNormalized.begin() + length_descriptorsNormalized,
                        [&factor](const Real& x) { return x * factor; });
                    length_descriptorsNormalized += desc.size();
                }
            }
        }
    }
}

void DescriptorCollector::normalizeDescriptorsSingleAtomCentralOrder(std::size_t atomIndex,
                                                                     const Real  norm)
{
    for (const String& key : descriptorKeyList)
    {
        if (descriptors.at(key)->get_isNormalized())
        {
            Real weight = descriptors.at(key)->get_weight();
            if (weight > (Real)0)
            {
                if (norm > constants::EPS_TOL)
                {
                    Real            factor = sqrt(weight) / norm;
                    const Vec1Real& desc = descriptors.at(key)->get_descriptor(atomIndex);
                    // rescale descriptors
                    Vec1Real& descriptorsNormalized =
                        (*this->descriptorsNormalized[key])[atomIndex];
                    Int& length_descriptorsNormalized =
                        (*this->length_descriptorsNormalized[key])[atomIndex];
                    std::for_each( // par_unseq,
                        desc.cbegin(),
                        desc.cend(),
                        [&](const Real& desc)
                        {
                            descriptorsNormalized[length_descriptorsNormalized] = desc * factor;
                            length_descriptorsNormalized++;
                        });
                }
            }
        }
    }
}

void DescriptorCollector::allocateArrays(void)
{
    if (algoExecution == ExecutionPolicy::cpuSingleCore) { allocateArraysCPU(); }
    else if (algoExecution == ExecutionPolicy::gpuStdLib) { allocateArraysGPU(); }
}

void DescriptorCollector::allocateArraysCPU(void)
{
    switch (storageOrder)
    {
    case DescriptorStorage::Type:
        for (const String& key : constants::descriptorKeyList)
        {
            const Vec1Int& nAtomsPerType = descriptors.at(key)->get_nAtomsType();
            length_descriptorsNormalized[key]->resize(nAtomsPerType.size());
            std::fill(length_descriptorsNormalized[key]->begin(),
                      length_descriptorsNormalized[key]->end(),
                      (Int)0.0);
            descriptorsNormalized[key]->resize(nAtomsPerType.size());
            // number of descriptors is same for all atoms for certain descriptor
            const Int numberDescriptors = descriptors.at(key)->get_sizeDescriptor(0);
            for (std::size_t type = 0; type < nAtomsPerType.size(); type++)
            {
                (*descriptorsNormalized[key])[type].resize(nAtomsPerType[type] * numberDescriptors);
                std::fill((*descriptorsNormalized[key])[type].begin(),
                          (*descriptorsNormalized[key])[type].end(),
                          (Real)0.0);
            }
        }
        break;
    case DescriptorStorage::CentralAtom:
        for (const String& key : constants::descriptorKeyList)
        {
            const std::size_t& nAtoms = descriptors.at(key)->get_nAtoms();
            descriptorsNormalized[key]->resize(nAtoms);
            length_descriptorsNormalized[key]->resize(nAtoms);
            std::fill(length_descriptorsNormalized[key]->begin(),
                      length_descriptorsNormalized[key]->end(),
                      (Int)0.0);
            const Int numberDescriptors = descriptors.at(key)->get_sizeDescriptor(0);
            for (std::size_t atom = 0; atom < nAtoms; atom++)
            {
                // number of descriptors is same for all atoms for certain descriptor
                (*descriptorsNormalized[key])[atom].resize(numberDescriptors);
            }
        }
        break;
    }
    Vec1Real& normsAtom = (*this->normsAtom);
    normsAtom.resize(numberAtoms);
}

void DescriptorCollector::allocateArraysGPU(void)
{

    bool resize = centralAtomIndexSize.checkResize(numberAtoms);
    if (resize)
    {
        normsAtom->resize(numberAtoms);
        centralAtomIndex.resize(numberAtoms);
        std::iota(centralAtomIndex.begin(),
                  centralAtomIndex.begin() + centralAtomIndexSize.actDim,
                  0);
    }

    switch (storageOrder)
    {
    case DescriptorStorage::Type:
        for (const String& key : constants::descriptorKeyList)
        {
            const Vec1Int& nAtomsPerType = descriptors.at(key)->get_nAtomsType();
            length_descriptorsNormalized[key]->resize(nAtomsPerType.size());
            // zero filling has to be shifted to gpu
            // number of descriptors is same for all atoms for certain descriptor
            const Int numberDescriptors = descriptors.at(key)->get_sizeDescriptor(0);
            resize = descriptorsNormalizedSize[key].checkResize1Dim(nAtomsPerType.size());
            if (resize)
            {
                descriptorsNormalizedSize[key].resizeArray1Dim(*(descriptorsNormalized[key]),
                                                               nAtomsPerType,
                                                               numberDescriptors);
                length_descriptorsNormalizedSize[key].maxDim = nAtomsPerType.size();
                length_descriptorsNormalizedSize[key].actDim = nAtomsPerType.size();
                length_descriptorsNormalized[key]->resize(nAtomsPerType.size());
            }
            else
            {
                length_descriptorsNormalizedSize[key].actDim = nAtomsPerType.size();
                resize = descriptorsNormalizedSize[key].checkResize2Dim(nAtomsPerType,
                                                                        numberDescriptors);
                if (resize)
                {
                    descriptorsNormalizedSize[key].resizeArray2Dim(*(descriptorsNormalized[key]),
                                                                   nAtomsPerType,
                                                                   numberDescriptors);
                }
            }
            std::fill(VASPML_PAR_UNSEQ
                      length_descriptorsNormalized[key]->begin(),
                      length_descriptorsNormalized[key]->begin()
                          + length_descriptorsNormalizedSize[key].actDim,
                      (Int)0.0);
        }
        break;
    case DescriptorStorage::CentralAtom:
        for (const String& key : constants::descriptorKeyList)
        {
            const std::size_t& nAtoms = descriptors.at(key)->get_nAtoms();
            resize = descriptorsNormalizedSize[key].checkResize1Dim(nAtoms);
            const Int numberDescriptors = descriptors.at(key)->get_sizeDescriptor(0);
            if (resize)
            {
                descriptorsNormalizedSize[key].resizeArray1Dim(*(descriptorsNormalized[key]),
                                                               numberDescriptors);
                length_descriptorsNormalizedSize[key].maxDim = nAtoms;
                length_descriptorsNormalizedSize[key].actDim = nAtoms;
            }
            else
            {
                resize = descriptorsNormalizedSize[key].checkResize2Dim(numberDescriptors);
                length_descriptorsNormalizedSize[key].actDim = nAtoms;
                if (resize)
                {
                    descriptorsNormalizedSize[key].resizeArray2Dim((*descriptorsNormalized[key]),
                                                                   numberDescriptors);
                }
            }
            length_descriptorsNormalized[key]->resize(nAtoms);
            std::fill(VASPML_PAR_UNSEQ
                      length_descriptorsNormalized[key]->begin(),
                      length_descriptorsNormalized[key]->begin()
                          + length_descriptorsNormalizedSize[key].actDim,
                      (Int)0.0);
        }
        break;
    }
}

const std::map<String, ShVec2Real>& DescriptorCollector::get_descriptorsNormalized(void) const
{
    return descriptorsNormalized;
}

Real DescriptorCollector::get_normAtom(std::size_t atomIndx) const
{
    return (*normsAtom)[atomIndx];
}

const Vec1Real& DescriptorCollector::get_normAtom(void) const
{
    return *normsAtom;
}

void DescriptorCollector::rearrangeSHS2Body(const TypeMap& typeMap)
{
    switch (algoExecution)
    {
    case ExecutionPolicy::cpuSingleCore:
        rearrangeSHS2BodyCPU(typeMap);
        break;
    case ExecutionPolicy::gpuStdLib:
        VASPML_PARALLEL( 
               rearrangeSHS2BodyGPU( typeMap );
        );
        break;
    }
}

void DescriptorCollector::rearrangeSHS2BodyCPU(const TypeMap& typeMap)
{

    if (typeMap.countStructureTypes() == typeMap.countForceFieldTypes()) { return; }
    // need to make copy here otherwise data will be overwritten before use
    const Vec2Real descriptorSHS2body = *descriptorsNormalized["SHS2-2-body"];
    const Int&     numberDescriptors =
        static_cast<DescriptorSHS2*>(descriptors["SHS2-2-body"].get())->get_nRootsOrder(0);
    const Vec1Int& nAtomsPerType = descriptors.at("SHS2-2-body")->get_nAtomsType();
    Vec2Real&      descriptorsNormalized = (*this->descriptorsNormalized["SHS2-2-body"]);

    std::size_t centralType = 0;
    std::for_each(
        descriptorsNormalized.begin(),
        descriptorsNormalized.begin() + typeMap.countStructureTypes(),
        [&](Vec1Real& slice)
        {
            slice.resize(nAtomsPerType[centralType] * typeMap.countForceFieldTypes()
                         * numberDescriptors);
            std::fill(slice.begin(), slice.end(), (Real)0);
            for (std::size_t centralAtom = 0; centralAtom < (std::size_t)nAtomsPerType[centralType];
                 centralAtom++)
            {
                for (std::size_t neighborTypeStruc = 0;
                     neighborTypeStruc < typeMap.countStructureTypes();
                     neighborTypeStruc++)
                {
                    std::size_t neighborTypeFF = typeMap.toType(neighborTypeStruc);
                    std::size_t neighborOffsetFF =
                        neighborTypeFF * numberDescriptors
                        + centralAtom * numberDescriptors * typeMap.countForceFieldTypes();
                    std::size_t neighborOffsetStruc =
                        neighborTypeStruc * numberDescriptors
                        + centralAtom * numberDescriptors * typeMap.countStructureTypes();
                    for (std::size_t desc = 0; desc < (std::size_t)numberDescriptors; desc++)
                    {
                        slice[neighborOffsetFF + desc] =
                            descriptorSHS2body[centralType][neighborOffsetStruc + desc];
                    }
                }
            }
            centralType++;
        });
}

void DescriptorCollector::rearrangeSHS2BodyGPU(const TypeMap& typeMap)
{

    if (typeMap.countStructureTypes() == typeMap.countForceFieldTypes()) { return; }
    // need to make copy here otherwise data will be overwritten before use
    const Vec2Real descriptorSHS2body = *descriptorsNormalized["SHS2-2-body"];
    const Int&     numberDescriptors =
        static_cast<DescriptorSHS2*>(descriptors["SHS2-2-body"].get())->get_nRootsOrder(0);
    Vec2Real&      descriptorsNormalized = (*this->descriptorsNormalized["SHS2-2-body"]);
    const Vec1Int& typeIndexCentral =
        descriptors.at("SHS2-2-body")->get_neighborList().get_typeIndexCentral();
    const Vec1Int& centralAtomIndexPerType =
        descriptors.at("SHS2-2-body")->get_neighborList().get_centralAtomIndexPerType();

    // allocate descriptorsNormalized for SHS2-body
    allocateDescriptorSHS2BodyNormalized(typeMap);

    // Now do the actual calculations over numberAtoms
    std::for_each( // seq,
        centralAtomIndex.cbegin(),
        centralAtomIndex.cbegin() + centralAtomIndexSize.actDim,
        [&](const std::size_t& centralAtom) mutable
        {
            std::size_t centralType = typeIndexCentral[centralAtom];
            std::size_t centralAtomInType = centralAtomIndexPerType[centralAtom];
            for (std::size_t neighborTypeStruc = 0;
                 neighborTypeStruc < typeMap.countStructureTypes();
                 neighborTypeStruc++)
            {
                std::size_t neighborTypeFF = typeMap.toType(neighborTypeStruc);
                std::size_t neighborOffsetFF =
                    neighborTypeFF * numberDescriptors
                    + centralAtomInType * numberDescriptors * typeMap.countForceFieldTypes();
                std::size_t neighborOffsetStruc =
                    neighborTypeStruc * numberDescriptors
                    + centralAtomInType * numberDescriptors * typeMap.countStructureTypes();
                for (std::size_t desc = 0; desc < (std::size_t)numberDescriptors; desc++)
                {
                    descriptorsNormalized[centralType][neighborOffsetFF + desc] =
                        descriptorSHS2body[centralType][neighborOffsetStruc + desc];
                }
            }
        });
}

void DescriptorCollector::allocateDescriptorSHS2BodyNormalized(const TypeMap& typeMap)
{

    const Vec1Int& nAtomsPerType = descriptors.at("SHS2-2-body")->get_nAtomsType();
    Vec2Real&      descriptorsNormalized = (*this->descriptorsNormalized["SHS2-2-body"]);
    const Int&     numberDescriptors =
        static_cast<DescriptorSHS2*>(descriptors["SHS2-2-body"].get())->get_nRootsOrder(0);
    std::size_t centralType = 0;
    std::for_each(descriptorsNormalized.begin(),
                  descriptorsNormalized.begin() + typeMap.countStructureTypes(),
                  [&](Vec1Real& slice)
                  {
                      slice.resize(nAtomsPerType[centralType] * typeMap.countForceFieldTypes()
                                   * numberDescriptors);
                      std::fill(slice.begin(), slice.end(), (Real)0);
                      centralType++;
                  });
}

const Descriptor& DescriptorCollector::getDescriptor(const std::string& key) const
{
    return (*descriptors.at(key));
}

void DescriptorCollector::writeDescriptorCollector(void) const
{

    for (const auto& [key, data] : descriptorsNormalized)
    {
        auto file = file_io::openFileO("DescriptorCollector_" + key + ".dat");
        for (const auto& x : *data)
        {
            for (const auto& y : x) file << str("%24.16E ", y) << std::endl;
        }
        file.close();
    }
}
