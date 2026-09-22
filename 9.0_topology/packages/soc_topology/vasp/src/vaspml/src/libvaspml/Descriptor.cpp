#include "Descriptor.hpp"
#include "DescriptorSHS2.hpp"
#include "DescriptorSHS3.hpp"
#include "DescriptorSHS3ReducedLinElem.hpp"

#include <limits>
#include <stdexcept>

using namespace vaspml;

Descriptor::Descriptor(void) :
    descriptorType(DescriptorType::none),
    algoExecution(ExecutionPolicy::cpuSingleCore)
{
    this->weight = (Real)0;
    this->isNormalized = false;
    this->derivativeMapComputed = false;
    this->neighborList = nullptr;
}

Descriptor::Descriptor(const Real            weight,
                       const bool            isNormalized,
                       const DescriptorType  descriptorType,
                       const ExecutionPolicy algoExecution) :
    descriptorType(descriptorType),
    algoExecution(algoExecution)
{
    this->weight = weight;
    this->isNormalized = isNormalized;
    this->derivativeMapComputed = false;
    this->neighborList = nullptr;
}

VASPML_NV_HOST_DEVICE
Real Descriptor::get_weight(void) const
{
    return weight;
}

VASPML_NV_HOST_DEVICE
bool Descriptor::get_isNormalized(void) const
{
    return isNormalized;
}

void Descriptor::set_neighborList(const std::shared_ptr<const NearestNeighborNSquare>& neighborList)
{
    this->neighborList = neighborList;
}

VASPML_NV_HOST_DEVICE
const Int& Descriptor::get_typeIndexCentral(const std::size_t atomIndex) const
{
    return neighborList->get_typeIndexCentral(atomIndex);
}

VASPML_NV_HOST_DEVICE
const Vec1Int& Descriptor::get_typeIndexCentral(void) const
{
    return neighborList->get_typeIndexCentral();
}

const Int& Descriptor::get_typeIndex(const std::size_t atomIndex, const std::size_t indxNeigh) const
{
    return neighborList->get_typeIndex(atomIndex, indxNeigh);
}

Int Descriptor::get_nAtoms(void) const
{
    return neighborList->get_nAtoms();
}

const Vec1Int& Descriptor::get_nAtomsType(void) const
{
    return neighborList->get_nAtomsType();
}

const NearestNeighborNSquare& Descriptor::get_neighborList(void) const
{
    return *neighborList;
}

std::shared_ptr<const NearestNeighborNSquare> Descriptor::get_neighborList_ptr(void) const
{
    return neighborList;
}

//
// implementation of member functions which were virtual
//
VASPML_NV_HOST_DEVICE
const Vec1Real& Descriptor::get_descriptor(const Size_t centralAtom) const
{
    if (descriptorType == DescriptorType::bodyOrder2)
    {
        return static_cast<const DescriptorSHS2*>(this)->get_descriptor(centralAtom);
    }
    else if (descriptorType == DescriptorType::bodyOrder3)
    {
        return static_cast<const DescriptorSHS3*>(this)->get_descriptor(centralAtom);
    }
    else if (descriptorType == DescriptorType::bodyOrder3LinearElement)
    {
        return static_cast<const DescriptorSHS3ReducedLinElem*>(this)->get_descriptor(centralAtom);
    }
#ifndef VASPML_USE_PSTL
    else throw std::runtime_error("ERROR: Invalid descriptor type.");
#endif
}

VASPML_NV_HOST_DEVICE
Int Descriptor::get_sizeDescriptor(const std::size_t centralAtom) const
{
    if (descriptorType == DescriptorType::bodyOrder2)
    {
        return static_cast<const DescriptorSHS2*>(this)->get_sizeDescriptor(centralAtom);
    }
    else if (descriptorType == DescriptorType::bodyOrder3)
    {
        return static_cast<const DescriptorSHS3*>(this)->get_sizeDescriptor(centralAtom);
    }
    else if (descriptorType == DescriptorType::bodyOrder3LinearElement)
    {
        return static_cast<const DescriptorSHS3ReducedLinElem*>(this)->get_sizeDescriptor(
            centralAtom);
    }
    else return -1;
}

void Descriptor::rescale_descriptor(const std::size_t centralAtom, const Real scaleFactor)
{
    if (descriptorType == DescriptorType::bodyOrder2)
    {
        return static_cast<DescriptorSHS2*>(this)->rescale_descriptor(centralAtom, scaleFactor);
    }
    else if (descriptorType == DescriptorType::bodyOrder3)
    {
        return static_cast<DescriptorSHS3*>(this)->rescale_descriptor(centralAtom, scaleFactor);
    }
    else if (descriptorType == DescriptorType::bodyOrder3LinearElement)
    {
        return static_cast<DescriptorSHS3ReducedLinElem*>(this)->rescale_descriptor(centralAtom,
                                                                                    scaleFactor);
    }
    else return;
}

void Descriptor::compute_forcePreContract(const Vec2Real& derivativeMatrix,
                                          Vec2Real&       forcePreContract,
                                          const TypeMap&  typeMap) const
{
    if (descriptorType == DescriptorType::bodyOrder2)
    {
        return static_cast<const DescriptorSHS2*>(this)->compute_forcePreContract(derivativeMatrix,
                                                                                  forcePreContract,
                                                                                  typeMap);
    }
    else if (descriptorType == DescriptorType::bodyOrder3)
    {
        return static_cast<const DescriptorSHS3*>(this)->compute_forcePreContract(derivativeMatrix,
                                                                                  forcePreContract,
                                                                                  typeMap);
    }
    else if (descriptorType == DescriptorType::bodyOrder3LinearElement)
    {
        return static_cast<const DescriptorSHS3ReducedLinElem*>(this)->compute_forcePreContract(
            derivativeMatrix,
            forcePreContract,
            typeMap);
    }
    else return;
}

Size_t Descriptor::get_forcePreContractSize(const std::size_t atomIndex) const
{
    if (descriptorType == DescriptorType::bodyOrder2)
    {
        return static_cast<const DescriptorSHS2*>(this)->get_forcePreContractSize(atomIndex);
    }
    else if (descriptorType == DescriptorType::bodyOrder3)
    {
        return static_cast<const DescriptorSHS3*>(this)->get_forcePreContractSize(atomIndex);
    }
    else if (descriptorType == DescriptorType::bodyOrder3LinearElement)
    {
        return static_cast<const DescriptorSHS3ReducedLinElem*>(this)->get_forcePreContractSize(
            atomIndex);
    }
    else return std::numeric_limits<Size_t>::max();
}

void Descriptor::computeForceTerms(const Vec2Real& forcePreContract,
                                   Vec2Real&       pairForces,
                                   Vec1Real&       centralForces) const
{
    if (descriptorType == DescriptorType::bodyOrder2)
    {
        return static_cast<const DescriptorSHS2*>(this)->computeForceTerms(forcePreContract,
                                                                           pairForces,
                                                                           centralForces);
    }
    else if (descriptorType == DescriptorType::bodyOrder3)
    {
        return static_cast<const DescriptorSHS3*>(this)->computeForceTerms(forcePreContract,
                                                                           pairForces,
                                                                           centralForces);
    }
    else if (descriptorType == DescriptorType::bodyOrder3LinearElement)
    {
        return static_cast<const DescriptorSHS3ReducedLinElem*>(this)->computeForceTerms(
            forcePreContract,
            pairForces,
            centralForces);
    }
    else return;
}
