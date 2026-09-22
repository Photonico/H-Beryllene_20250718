#include "nearest_neighbor.hpp"
#include "ParallelEnvironemt.hpp"
#include "utils.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>

using namespace vaspml;

NearestNeighborNSquare::NearestNeighborNSquare(void)
{
    //   globalIndex = nullptr;
}

NearestNeighborNSquare::NearestNeighborNSquare(Real              Input_CutOff,
                                               bool              isTypeSorted_in,
                                               bool              isDistSorted_in,
                                               const ShVec2Int&  globalIndex,
                                               const ShVec2Int&  typeIndex,
                                               const ShVec1Int&  typeIndexCentral,
                                               const ShVec2Real& distances,
                                               const ShVec2Real& connectionVector,
                                               const ShVec2Real& connectionVectorNormalized,
                                               const ShVec1Int&  numberNeighbors,
                                               const ShVec2Int&  numberNeighborsType,
                                               const ShVec1Int&  nAtomsType,
                                               const ShVec1Int&  centralAtomIndexPerType)
{
    if (globalIndex == nullptr) { this->globalIndex = std::make_shared<Vec2Int>(); }
    else { this->globalIndex = globalIndex; }
    if (typeIndex == nullptr) { this->typeIndex = std::make_shared<Vec2Int>(); }
    else { this->typeIndex = typeIndex; }
    if (typeIndexCentral == nullptr) { this->typeIndexCentral = std::make_shared<Vec1Int>(); }
    else { this->typeIndexCentral = typeIndexCentral; }
    if (distances == nullptr) { this->distances = std::make_shared<Vec2Real>(); }
    else { this->distances = distances; }
    if (connectionVector == nullptr) { this->connectionVector = std::make_shared<Vec2Real>(); }
    else { this->connectionVector = connectionVector; }
    if (connectionVectorNormalized == nullptr)
    {
        this->connectionVectorNormalized = std::make_shared<Vec2Real>();
    }
    else { this->connectionVectorNormalized = connectionVectorNormalized; }
    if (numberNeighbors == nullptr) { this->numberNeighbors = std::make_shared<Vec1Int>(); }
    else { this->numberNeighbors = numberNeighbors; }
    if (numberNeighborsType == nullptr) { this->numberNeighborsType = std::make_shared<Vec2Int>(); }
    else { this->numberNeighborsType = numberNeighborsType; }
    if (nAtomsType == nullptr) { this->nAtomsType = std::make_shared<Vec1Int>(); }
    else { this->nAtomsType = nAtomsType; }
    if (centralAtomIndexPerType == nullptr)
    {
        this->centralAtomIndexPerType = std::make_shared<Vec1Int>();
    }
    else { this->centralAtomIndexPerType = centralAtomIndexPerType; }

    cutOff = Input_CutOff;
    cutOffSquared = Input_CutOff * Input_CutOff;
    isTypeSorted = isTypeSorted_in;
    isDistSorted = isDistSorted_in;
    eps = (Real)1e-6;
    if (isTypeSorted && isDistSorted)
    {
        //        global_scope::tutor.bug( "ERROR: NearestNeighborNSquare::NearestNeighborNSquare
        //        \n"
        //                                 "I STILL HAVE TO IMPLEMENT THIS " );
    }
}

void NearestNeighborNSquare::computeNearestNeighborsDirectCoordinates(const Structure& positions)
{
    nAtoms = positions.get_numAtoms();
    allocateStorageArrays(nAtoms);
    *typeIndexCentral = *positions.get_types();
    compute_unique_types(*typeIndexCentral);
    for (std::size_t i = 0; i < nAtoms; i++)
    {
        computeNearestNeighborsSingleAtomDirect(i, positions);
    }
    compute_centralAtomIndexPerType();
}

void NearestNeighborNSquare::computeNearestNeighborsCartesianCoordinates(const Structure& positions)
{
    nAtoms = positions.get_numAtoms();
    allocateStorageArrays(nAtoms);
    typeIndexCentral = positions.get_types();
    compute_unique_types(*typeIndexCentral);
    for (std::size_t i = 0; i < nAtoms; i++)
    {
        computeNearestNeighborsSingleAtomCartesian(i, positions);
    }
    compute_centralAtomIndexPerType();
}

void NearestNeighborNSquare::compute_centralAtomIndexPerType(void)
{
    Vec1Int&    nAtomsType = *(this->nAtomsType);
    Vec1Int&    centralAtomIndexPerType = *(this->centralAtomIndexPerType);
    std::size_t counter = 0;
    for (std::size_t type = 0; type < (std::size_t)nTypes; type++)
    {
        for (std::size_t atomPerType = 0; atomPerType < (std::size_t)nAtomsType[type];
             atomPerType++)
        {
            centralAtomIndexPerType[counter] = atomPerType;
            counter++;
        }
    }
}

void NearestNeighborNSquare::allocateStorageArrays(std::size_t size)
{
    globalIndex->resize(size);
    typeIndex->resize(size);
    distances->resize(size);
    connectionVector->resize(size);
    connectionVectorNormalized->resize(size);
    numberNeighbors->resize(size);

    typeIndexCentral->resize(size);
    typeEnd = std::make_shared<Vec2Int>();
    typeEnd->resize(size);
    typeStart = std::make_shared<Vec2Int>();
    typeStart->resize(size);
    centralAtomIndexPerType->resize(size);
}

void NearestNeighborNSquare::allocateWorkArrays(Vec1Int&    nn_list,
                                                Vec1Int&    types_index,
                                                Vec1Real&   dist,
                                                Vec1Real&   dist_vecs,
                                                Vec1Real&   normed_dist,
                                                std::size_t N)
{
    nn_list.resize(N);
    types_index.resize(N);
    dist.resize(N);
    dist_vecs.resize(3 * N);
    normed_dist.resize(3 * N);
}

void NearestNeighborNSquare::computeNearestNeighborsSingleAtomDirect(const std::size_t indx,
                                                                     const Structure&  positions)
{
    Vec1Int xyz_min(3);
    Vec1Int xyz_max(3);

    auto& [pos_x, pos_y, pos_z] = positions.getAtom(indx);
    Vec1Real pos_central(3);
    pos_central[0] = pos_x;
    pos_central[1] = pos_y;
    pos_central[2] = pos_z;

    const Lattice& lattice = positions.get_lattice();
    const Lattice& inverse_lattice = positions.get_inverseLattice();
    latticeVolume = lattice.get_volume();

    lattice.timesVectorInPlace(pos_central);

    set_periodic_size_images(pos_central, inverse_lattice, xyz_min, xyz_max);

    // position central atom
    // positions neighbor atom
    Real neighbor_x;
    Real neighbor_y;
    Real neighbor_z;
    // distance vector between central atom and neighbor
    // nearest image convention
    Real delta_x;
    Real delta_y;
    Real delta_z;
    // compute distance vector in cartesian coordintes
    Real distance;

    // define temporrary arrays to store nearest neighbor data
    // for certain atom
    Vec1Int  nn_list;
    Vec1Int  types_index;
    Vec1Real dist;
    Vec1Real dist_vecs;
    Vec1Real normed_dist;

    std::size_t multi = (xyz_max[0] - xyz_min[0] + 1) * (xyz_max[1] - xyz_min[1] + 1)
                      * (xyz_max[2] - xyz_min[2] + 1);
    // allocate work arrays
    allocateWorkArrays(nn_list, types_index, dist, dist_vecs, normed_dist, multi * nAtoms);

    std::size_t counter = 0;
    for (int ix = xyz_min[0]; ix <= xyz_max[0]; ix++)
    {
        for (int iy = xyz_min[1]; iy <= xyz_max[1]; iy++)
        {
            for (int iz = xyz_min[2]; iz <= xyz_max[2]; iz++)
            {
                for (std::size_t n_indx = 0; n_indx < nAtoms; n_indx++)
                {
                    //NOTE following workaround is needed to get
                    //compiled code with nc++ (NCC) 5.0.1 (Build 16:03:10 May  8 2023)
                    //auto [ x, y, z ] = positions.get_Atom( n_indx );
                    std::tuple<const Real&, const Real&, const Real&> tmp =
                        positions.getAtom(n_indx);
                    auto [x, y, z] = tmp;

                    neighbor_x = x + (Real)ix;
                    neighbor_y = y + (Real)iy;
                    neighbor_z = z + (Real)iz;

                    delta_x = neighbor_x - pos_x;
                    delta_y = neighbor_y - pos_y;
                    delta_z = neighbor_z - pos_z;
                    lattice.timesVectorInPlace(delta_x, delta_y, delta_z);

                    distance = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;

                    // collecting data
                    if (distance <= cutOffSquared && distance > eps)
                    {
                        nn_list[counter] = n_indx;
                        types_index[counter] = (*typeIndexCentral)[n_indx];
                        dist[counter] = std::sqrt(distance);
                        dist_vecs[3 * counter] = delta_x;
                        dist_vecs[3 * counter + 1] = delta_y;
                        dist_vecs[3 * counter + 2] = delta_z;
                        normed_dist[3 * counter] = delta_x / dist[counter];
                        normed_dist[3 * counter + 1] = delta_y / dist[counter];
                        normed_dist[3 * counter + 2] = delta_z / dist[counter];
                        counter++;
                    }
                }
            }
        }
    }
    fillNeighborArrays(indx, counter, nn_list, types_index, dist, dist_vecs, normed_dist);
}

void NearestNeighborNSquare::computeNearestNeighborsSingleAtomCartesian(const std::size_t indx,
                                                                        const Structure&  positions)
{
    Vec1Int xyz_min(3);
    Vec1Int xyz_max(3);
    auto& [pos_x, pos_y, pos_z] = positions.getAtom(indx);

    Vec1Real pos_central(3);
    pos_central[0] = pos_x;
    pos_central[1] = pos_y;
    pos_central[2] = pos_z;
    const Lattice& lattice = positions.get_lattice();
    const Lattice& inverse_lattice = positions.get_inverseLattice();

    set_periodic_size_images(pos_central, inverse_lattice, xyz_min, xyz_max);

    Real neighbor_x;
    Real neighbor_y;
    Real neighbor_z;
    // distance vector between central atom and neighbor
    // nearest image convention
    Real delta_x;
    Real delta_y;
    Real delta_z;
    // compute distance vector in cartesian coordintes
    Real distance;

    // define temporrary arrays to store nearest neighbor data
    // for certain atom
    Vec1Int  nn_list;
    Vec1Int  types_index;
    Vec1Real dist;
    Vec1Real dist_vecs;
    Vec1Real normed_dist;

    Real shift_x;
    Real shift_y;
    Real shift_z;

    std::size_t multi = (xyz_max[0] - xyz_min[0] + 1) * (xyz_max[1] - xyz_min[1] + 1)
                      * (xyz_max[2] - xyz_min[2] + 1);

    // allocate work arrays
    allocateWorkArrays(nn_list, types_index, dist, dist_vecs, normed_dist, multi * nAtoms);
    std::size_t counter = 0;
    for (int ix = xyz_min[0]; ix <= xyz_max[0]; ix++)
    {
        for (int iy = xyz_min[1]; iy <= xyz_max[1]; iy++)
        {
            for (int iz = xyz_min[2]; iz <= xyz_max[2]; iz++)
            {
                for (std::size_t n_indx = 0; n_indx < nAtoms; n_indx++)
                {

                    //NOTE following workaround is needed to get
                    //compiled code with nc++ (NCC) 5.0.1 (Build 16:03:10 May  8 2023)
                    //auto [ x, y, z ] = positions.get_Atom( n_indx );
                    std::tuple<const Real&, const Real&, const Real&> tmp =
                        positions.getAtom(n_indx);
                    auto [x, y, z] = tmp;

                    shift_x = (Real)ix;
                    shift_y = (Real)iy;
                    shift_z = (Real)iz;

                    lattice.timesVectorInPlace(shift_x, shift_y, shift_z);

                    neighbor_x = x + shift_x;
                    neighbor_y = y + shift_y;
                    neighbor_z = z + shift_z;

                    delta_x = neighbor_x - pos_x;
                    delta_y = neighbor_y - pos_y;
                    delta_z = neighbor_z - pos_z;
                    distance = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;
                    if (distance <= cutOffSquared && distance > eps)
                    {
                        nn_list[counter] = n_indx;
                        types_index[counter] = (*typeIndexCentral)[n_indx];
                        dist[counter] = std::sqrt(distance);
                        dist_vecs[3 * counter] = delta_x;
                        dist_vecs[3 * counter + 1] = delta_y;
                        dist_vecs[3 * counter + 2] = delta_z;
                        normed_dist[3 * counter] = delta_x / dist[counter];
                        normed_dist[3 * counter + 1] = delta_y / dist[counter];
                        normed_dist[3 * counter + 2] = delta_z / dist[counter];
                        counter++;
                    }
                }
            }
        }
    }

    fillNeighborArrays(indx, counter, nn_list, types_index, dist, dist_vecs, normed_dist);
}

void NearestNeighborNSquare::set_periodic_size_images(const Vec1Real& pos,
                                                      const Lattice&  inverse_lattice,
                                                      Vec1Int&        imin,
                                                      Vec1Int&        imax)
{
    Vec1Real rr(3);
    Vec1Real xyz_max(9);
    Vec1Real xyz_min(9);

    for (std::size_t ixyz = 0; ixyz < 3; ixyz++)
    {
        rr[ixyz] =
            std::sqrt(inverse_lattice.component(0, ixyz) * inverse_lattice.component(0, ixyz)
                      + inverse_lattice.component(1, ixyz) * inverse_lattice.component(1, ixyz)
                      + inverse_lattice.component(2, ixyz) * inverse_lattice.component(2, ixyz));
    }
    std::size_t counter = 0;
    for (std::size_t ixyz = 0; ixyz < 3; ixyz++)
    {
        for (std::size_t jxyz = 0; jxyz < 3; jxyz++)
        {
            xyz_max[counter] =
                pos[jxyz] + cutOff * inverse_lattice.component(jxyz, ixyz) / rr[jxyz];
            xyz_min[counter] =
                pos[jxyz] - cutOff * inverse_lattice.component(jxyz, ixyz) / rr[jxyz];
            counter++;
        }
    }

    for (std::size_t ixyz = 0; ixyz < 3; ixyz++)
    {
        imax[ixyz] = ceil(inverse_lattice.component(0, ixyz) * xyz_max[3 * ixyz]
                          + inverse_lattice.component(1, ixyz) * xyz_max[3 * ixyz + 1]
                          + inverse_lattice.component(2, ixyz) * xyz_max[3 * ixyz + 2])
                   - 1;
        imin[ixyz] = ceil(inverse_lattice.component(0, ixyz) * xyz_min[3 * ixyz]
                          + inverse_lattice.component(1, ixyz) * xyz_min[3 * ixyz + 1]
                          + inverse_lattice.component(2, ixyz) * xyz_min[3 * ixyz + 2])
                   - 1;
    }

    return;
}

void NearestNeighborNSquare::fillNeighborArrays(const std::size_t index,
                                                const std::size_t nelements,
                                                const Vec1Int&    nn_list,
                                                const Vec1Int&    types_index,
                                                const Vec1Real&   dist,
                                                const Vec1Real&   dist_vecs,
                                                const Vec1Real&   normed_dist)
{
    std::vector<std::size_t> indx_array;
    if (isTypeSorted) { indx_array = vector_tools::argSortSlice(types_index, 0, nelements); }
    else if (isDistSorted) { indx_array = vector_tools::argSortSlice(dist, 0, nelements); }
    else { indx_array = vector_tools::generateIntSequence((std::size_t)0, (std::size_t)nelements); }

    (*globalIndex)[index].resize(nelements);
    (*typeIndex)[index].resize(nelements);
    (*distances)[index].resize(nelements);
    (*connectionVector)[index].resize(3 * nelements);
    (*connectionVectorNormalized)[index].resize(3 * nelements);

    int type = types_index[indx_array[0]];
    for (std::size_t i = 0; i < (std::size_t)indx_array.size(); i++)
    {
        std::size_t myElement = indx_array[i];
        (*globalIndex)[index][i] = nn_list[myElement];
        (*typeIndex)[index][i] = types_index[myElement];
        (*distances)[index][i] = dist[myElement];
        (*connectionVector)[index][i * 3] = dist_vecs[myElement * 3];
        (*connectionVector)[index][i * 3 + 1] = dist_vecs[myElement * 3 + 1];
        (*connectionVector)[index][i * 3 + 2] = dist_vecs[myElement * 3 + 2];
        (*connectionVectorNormalized)[index][i * 3] = normed_dist[myElement * 3];
        (*connectionVectorNormalized)[index][i * 3 + 1] = normed_dist[myElement * 3 + 1];
        (*connectionVectorNormalized)[index][i * 3 + 2] = normed_dist[myElement * 3 + 2];
        (*numberNeighborsType)[index][unique_types[(*typeIndex)[index][i]]]++;
        if (isTypeSorted)
        {
            if ((*typeIndex)[index][i] != type)
            {
                (*typeEnd)[index][unique_types[type]] = i;
                type = (*typeIndex)[index][i];
            }
        }
    }

    (*numberNeighbors)[index] = nelements;
    for (std::size_t i = 1; i < (*typeStart)[index].size(); i++)
    {
        (*typeStart)[index][i] = (*typeEnd)[index][i - 1];
    }
}

void NearestNeighborNSquare::compute_unique_types(const Vec1Int& types)
{
    Vec1Int copy_types(types.size());
    for (std::size_t i = 0; i < types.size(); i++) { copy_types[i] = types[i]; }
    vector_tools::get_unique(copy_types);
    numberNeighborsType->resize(types.size());
    unique_types.clear();
    types_unique.clear();

    for (std::size_t i = 0; i < copy_types.size(); i++)
    {
        unique_types.insert({copy_types[i], (Int)i});
        types_unique.insert({(Int)i, copy_types[i]});
    }

    nTypes = copy_types.size();
    for (std::size_t i = 0; i < numberNeighborsType->size(); i++)
    {
        (*numberNeighborsType)[i].resize(nTypes);
        (*typeEnd)[i].resize(nTypes);
        (*typeStart)[i].resize(nTypes);
    }

    nAtomsType->resize(nTypes);
    // count number of atoms per type
    copy_types.resize(types.size());
    for (std::size_t i = 0; i < (std::size_t)types.size(); i++) { copy_types[i] = types[i]; }

    for (std::size_t i = 0; i < (std::size_t)nTypes; i++)
    {
        (*nAtomsType)[i] = vector_tools::count_element(copy_types, unique_types[i]);
    }
}

void NearestNeighborNSquare::writeListToScreen(void) const
{
    for (std::size_t i = 0; i < globalIndex->size(); i++)
    {
        //Vec1Real d = (*distances)[i];
        //Vec1Int d = (*typeIndex)[i];
        //std::sort(d.begin(), d.end());
        std::cout << "Central atom " << std::setw(7) << i << " |   ";
        for (std::size_t j = 0; j < (*globalIndex)[i].size(); j++)
        {
            std::cout << std::setw(7) << (*globalIndex)[i][j];
            //std::cout << str(" %5.3f ", (*distances)[i][j]);
            //std::cout << str(" %24.16E ", d[j]);
            //std::cout << str(" %5.3f ", d[j]);
            //std::cout << str(" %3d ", d[j]);
        }
        std::cout << std::endl;
    }
}

const Int& NearestNeighborNSquare::get_globalIndex(const std::size_t atomIndx,
                                                   const std::size_t indxNeigh) const
{
    return (*globalIndex)[atomIndx][indxNeigh];
}

const Vec2Int& NearestNeighborNSquare::get_globalIndex(void) const
{
    return *globalIndex;
}

const Real& NearestNeighborNSquare::get_distances(const std::size_t atomIndx,
                                                  const std::size_t indxNeigh) const
{
    return (*distances)[atomIndx][indxNeigh];
}

const Vec2Real& NearestNeighborNSquare::get_distances(void) const
{
    return *distances;
}

VASPML_NV_HOST_DEVICE
const Vec1Real& NearestNeighborNSquare::get_distances(std::size_t atomIndx) const
{
    return (*distances)[atomIndx];
}

const Real& NearestNeighborNSquare::get_connectionVector_x(const std::size_t atomIndx,
                                                           const std::size_t indxNeigh) const
{
    return (*connectionVector)[atomIndx][3 * indxNeigh];
}

const Real& NearestNeighborNSquare::get_connectionVector_y(const std::size_t atomIndx,
                                                           const std::size_t indxNeigh) const
{
    return (*connectionVector)[atomIndx][3 * indxNeigh + 1];
}

const Real& NearestNeighborNSquare::get_connectionVector_z(const std::size_t atomIndx,
                                                           const std::size_t indxNeigh) const
{
    return (*connectionVector)[atomIndx][3 * indxNeigh + 2];
}

const std::tuple<const Real&, const Real&, const Real&>
NearestNeighborNSquare::get_connectionVector(const std::size_t atomIndx,
                                             const std::size_t indxNeigh) const
{
    return std::tie(get_connectionVector_x(atomIndx, indxNeigh),
                    get_connectionVector_y(atomIndx, indxNeigh),
                    get_connectionVector_z(atomIndx, indxNeigh));
}

const Vec2Real& NearestNeighborNSquare::get_connectionVector(void) const
{
    return *connectionVector;
}

const Vec1Real& NearestNeighborNSquare::get_connectionVector(std::size_t atomIndx) const
{
    return (*connectionVector)[atomIndx];
}

const Real& NearestNeighborNSquare::get_connectionVectorNormalized_x(const std::size_t atomIndx,
                                                                     std::size_t indxNeigh) const
{
    return (*connectionVectorNormalized)[atomIndx][3 * indxNeigh];
}

const Real& NearestNeighborNSquare::get_connectionVectorNormalized_y(const std::size_t atomIndx,
                                                                     std::size_t indxNeigh) const
{
    return (*connectionVectorNormalized)[atomIndx][3 * indxNeigh + 1];
}

const Real& NearestNeighborNSquare::get_connectionVectorNormalized_z(const std::size_t atomIndx,
                                                                     std::size_t indxNeigh) const
{
    return (*connectionVectorNormalized)[atomIndx][3 * indxNeigh + 2];
}

const std::tuple<const Real&, const Real&, const Real&>
NearestNeighborNSquare::get_connectionVectorNormalized(const std::size_t atomIndx,
                                                       const std::size_t indxNeigh) const
{
    return std::tie(get_connectionVectorNormalized_x(atomIndx, indxNeigh),
                    get_connectionVectorNormalized_y(atomIndx, indxNeigh),
                    get_connectionVectorNormalized_z(atomIndx, indxNeigh));
}

const Vec2Real& NearestNeighborNSquare::get_connectionVectorNormalized(void) const
{
    return *connectionVectorNormalized;
}

VASPML_NV_HOST_DEVICE
const Vec1Real& NearestNeighborNSquare::get_connectionVectorNormalized(std::size_t atomIndx) const
{
    return (*connectionVectorNormalized)[atomIndx];
}

VASPML_NV_HOST_DEVICE
std::size_t NearestNeighborNSquare::get_size(const std::size_t atomIndx) const
{
    return (*numberNeighbors)[atomIndx];
}

const Vec1Int& NearestNeighborNSquare::get_size(void) const
{
    return *numberNeighbors;
}

VASPML_NV_HOST_DEVICE
const Int& NearestNeighborNSquare::get_typeIndexCentral(const std::size_t atomIndex) const
{
    return (*typeIndexCentral)[atomIndex];
}

const Vec1Int& NearestNeighborNSquare::get_typeIndexCentral(void) const
{
    return *typeIndexCentral;
}

VASPML_NV_HOST_DEVICE
const Int& NearestNeighborNSquare::get_typeIndex(const std::size_t atomIndex,
                                                 const std::size_t indxNeigh) const
{
    return (*typeIndex)[atomIndex][indxNeigh];
}

VASPML_NV_HOST_DEVICE
const Vec1Int& NearestNeighborNSquare::get_typeIndex(const std::size_t atomIndex) const
{
    return (*typeIndex)[atomIndex];
}

const Vec2Int& NearestNeighborNSquare::get_typeIndex(void) const
{
    return *typeIndex;
}

const std::tuple<const Int&,
                 const Real&,
                 const Real&,
                 const Real&,
                 const Real&,
                 const Real&,
                 const Real&,
                 const Real&>
NearestNeighborNSquare::get_neighborData(const std::size_t atomIndex,
                                         const std::size_t indxNeigh) const
{

    return std::tie((*typeIndex)[atomIndex][indxNeigh],
                    (*distances)[atomIndex][indxNeigh],
                    get_connectionVector_x(atomIndex, indxNeigh),
                    get_connectionVector_y(atomIndex, indxNeigh),
                    get_connectionVector_z(atomIndex, indxNeigh),
                    get_connectionVectorNormalized_x(atomIndex, indxNeigh),
                    get_connectionVectorNormalized_y(atomIndex, indxNeigh),
                    get_connectionVectorNormalized_z(atomIndex, indxNeigh));
}

std::size_t NearestNeighborNSquare::get_numberAtoms(void) const
{
    return typeIndex->size();
}

bool NearestNeighborNSquare::is_typeSorted(void) const
{
    return isTypeSorted;
}

std::size_t NearestNeighborNSquare::get_numberTypes(void) const
{
    return nTypes;
}

const Vec1Int& NearestNeighborNSquare::get_nAtomsType(void) const
{
    return *nAtomsType;
}

const Int& NearestNeighborNSquare::get_nAtomsType(const std::size_t type) const
{
    return (*nAtomsType)[type];
}

const std::size_t& NearestNeighborNSquare::get_nAtoms(void) const
{
    return nAtoms;
}

const Real& NearestNeighborNSquare::get_latticeVolume(void) const
{
    return latticeVolume;
}

const Vec1Int& NearestNeighborNSquare::get_centralAtomIndexPerType(void) const
{
    return *centralAtomIndexPerType;
}

Real NearestNeighborNSquare::get_cutOff() const
{
    return cutOff;
}

void NearestNeighborNSquare::set_nAtoms(const Int nAtoms)
{
    this->nAtoms = nAtoms;
}

void NearestNeighborNSquare::set_nTypes(const Int nTypes)
{
    this->nTypes = nTypes;
}
