#ifndef NEAREST_NEIGHBOR_HPP
#define NEAREST_NEIGHBOR_HPP

#include "Lattice.hpp"
#include "Structure.hpp"
#include "types.hpp"

#include <cstddef>
#include <map>
#include <tuple>

namespace vaspml
{

/// Nsquare nearest neighbor algorithm
class NearestNeighborNSquare
{

  public:
    NearestNeighborNSquare(void);
    NearestNeighborNSquare(Real              Cutoff,
                           bool              isTypeSorted_in = false,
                           bool              isDistSorted_in = false,
                           const ShVec2Int&  globalIndex = nullptr,
                           const ShVec2Int&  typeIndex = nullptr,
                           const ShVec1Int&  typeIndexCentral = nullptr,
                           const ShVec2Real& distances = nullptr,
                           const ShVec2Real& connectionVector = nullptr,
                           const ShVec2Real& connectionVectorNormalized = nullptr,
                           const ShVec1Int&  numberNeighbors = nullptr,
                           const ShVec2Int&  numberNeighborsType = nullptr,
                           const ShVec1Int&  nAtomsType = nullptr,
                           const ShVec1Int&  centralAtomIndexPerType = nullptr);

    /** compute nearest neighbor arrays in terms of direct coordinates
     *  the results, distance vectors, distance, normalized distance
     *  will still be in cartesian coordinates
     *
     * @param positions strucure containing information as positions, lattice, atom types
     */
    void computeNearestNeighborsDirectCoordinates(const Structure& positions);

    /** compute nearest neighbor arrays in terms of Cartesian coordinates
     *  the results, distance vectors, distance, normalized distance
     *  will be in cartesian coordinates
     *
     * @param positions strucure containing information as positions, lattice, atom types
     */
    void computeNearestNeighborsCartesianCoordinates(const Structure& positions);

    /// write nearest neighbor list to screen
    void writeListToScreen(void) const;

    /*******************************************************************************************
     * get global index of neighbor in position array
     * @param atomIndx    central atom of which a neigbor is requested
     * @param indxNeigh   number of neighor in nearest neighbor array to return
     *******************************************************************************************/
    const Int& get_globalIndex(const std::size_t atomIndx, const std::size_t indxNeigh) const;

    /**  get whole global index array of whole nearest neighbor list
     */
    const Vec2Int& get_globalIndex(void) const;

    /** get distance between atoms atomIndx and atom get_index_list( atomIndx, inxdNeigh )
     *  distances are ordered in the same way as index_list
     *@param atomIndx   -> central atom from which a distance to
     *@param indxNeigh  -> points to the list entry of the neighbor
     */
    const Real& get_distances(const std::size_t atomIndx, const std::size_t indxNeigh) const;
    /**  get all distances between atoms in neighbor list
     *
     * @note format of return array is [central_atom][ vector of neighbor distances ]
     */
    const Vec2Real& get_distances(void) const;

    /**
     * get distances of nearest neighbors for given central atom
     *
     * @param atomIndx
     */
    const Vec1Real& get_distances(std::size_t atomIndx) const;

    /**
     * get x component of distance vector between atom pair
     *
     * distances are ordered in the same way as index_list
     *
     *@param atomIndx  central atom from which a distance to
     *@param indxNeigh points to the list entry of the neighbor
     */
    const Real& get_connectionVector_x(const std::size_t atomIndx, std::size_t indxNeigh) const;

    /**  get y component of distance vector between atoms atomIndx and atom get_index_list(
     * atomIndx, inxdNeigh ) distances are ordered in the same way as index_list
     * @param atomIndx  central atom from which a distance to
     * @param indxNeigh points to the list entry of the neighbor
     */
    const Real& get_connectionVector_y(const std::size_t atomIndx, std::size_t indxNeigh) const;

    /**  get z componeent of distance vector between atoms atomIndx and atom get_index_list(
     *atomIndx, inxdNeigh ) distances are ordered in the same way as index_list
     *@param atomIndx   -> central atom from which a distance to
     *@param indxNeigh  -> points to the list entry of the neighbor
     */
    const Real& get_connectionVector_z(const std::size_t atomIndx, std::size_t indxNeigh) const;

    /**
     * get x,y,z component of distance vector between atom pair
     *
     * distances are ordered in the same way as index_list for optimal
     * performance call as auto & [ x, y, z ] = get_connectionVector( atomIndx, indxNeigh )
     *
     *@param atomIndx  central atom from which a distance to
     *@param indxNeigh points to the list entry of the neighbor
     */
    const std::tuple<const Real&, const Real&, const Real&> get_connectionVector(
        std::size_t atomIndx,
        std::size_t indxNeigh) const;

    /**  get x,y,z component of distance vector between all atoms
     *
     * first index is the central atom
     * format is xyz atom1 neighbor 1, xyz atom1 neighbor 2,...., xyz atom2, neighbor1...
     */
    const Vec2Real& get_connectionVector(void) const;

    /**
     * get connection vectors for cetrain central atom defined by atomIndx
     *
     * data is of format xyz neighbo1, xyz neighbor2, xyz neighbor3,...
     * @param atomIndx list index of central atom
     */
    const Vec1Real& get_connectionVector(std::size_t atomIndx) const;

    /** get x component of normalized distance vector between atom pair
     *
     *  distances are ordered in the same way as index_list
     *
     *@param atomIndx  central atom from which a distance to
     *@param indxNeigh points to the list entry of the neighbor
     */
    const Real& get_connectionVectorNormalized_x(const std::size_t atomIndx,
                                                 std::size_t       indxNeigh) const;

    /**  get y component of normalized distance vector between atom pairs
     *
     * distances are ordered in the same way as index_list
     * @param atomIndx  central atom from which a distance to
     * @param indxNeigh points to the list entry of the neighbor
     */
    const Real& get_connectionVectorNormalized_y(const std::size_t atomIndx,
                                                 std::size_t       indxNeigh) const;

    /**  get z component of normalized distance vector between atom pairs
     *
     *   distances are ordered in the same way as index_list
     *@param atomIndx  central atom from which a distance to
     *@param indxNeigh points to the list entry of the neighbor
     */
    const Real& get_connectionVectorNormalized_z(const std::size_t atomIndx,
                                                 std::size_t       indxNeigh) const;

    /**  get x,y,z component of distance vector between atom pairs
     *
     * distances are ordered in the same way as index_list for optimal
     * performance call as auto & [ x, y, z ] = get_connectionVector( atomIndx, indxNeigh )
     *
     * @param atomIndx  central atom from which a distance to
     * @param indxNeigh points to the list entry of the neighbor
     */
    const std::tuple<const Real&, const Real&, const Real&> get_connectionVectorNormalized(
        std::size_t atomIndx,
        std::size_t indxNeigh) const;

    /**  get x,y,z component of distance vector between all atoms
     *  format is xyz atom1 neighbor 1, xyz atom1 neighbor 2,...., xyz atom2, neighbor1...
     */
    const Vec2Real& get_connectionVectorNormalized(void) const;

    /**
     * get normalized connection vectors for certain central atom
     *
     * @param atomIndx atom index of desired central atom
     */
    const Vec1Real& get_connectionVectorNormalized(std::size_t atomIndx) const;

    /**
     * get number of nearest neighbors for atom atomIndx
     *@param atomIndx central atom index
     */
    std::size_t get_size(const std::size_t atomIndx) const;

    /**
     * get number of nearest neighbors whole array
     */
    const Vec1Int& get_size(void) const;

    /*******************************************************************************************
     * get type index of central atom in nearest neighbor list
     *
     * @param atomIndex of central atom from which type has to be retrieved
     *******************************************************************************************/
    const Int& get_typeIndexCentral(const std::size_t atomIndex) const;
    /*******************************************************************************************
     * get types of all central atoms
     *
     * for performance call with const std::vector<Int>& typeIndexCentral = get_typeIndexCentral();
     *******************************************************************************************/
    const Vec1Int& get_typeIndexCentral(void) const;

    /**
     * get type index of neighbor atom in nearest neighbor list
     *@param atomIndex central atom index from which neighbor the types is retrieved
     *@param indxNeigh number of neighbor atom from which type is wanted
     */
    const Int& get_typeIndex(const std::size_t atomIndex, const std::size_t indxNeigh) const;
    /*******************************************************************************************
     * get type index of neighbors for given central atom
     * @param atomIndex central atom index from which neighbor the types is retrieved
     *******************************************************************************************/
    const Vec1Int& get_typeIndex(const std::size_t atomIndex) const;
    /**
     * get types of all neighbor atoms
     *
     * for performance call with const std::vector<std::vector<Int>>& typeIndex = get_typeIndex();
     */
    const Vec2Int& get_typeIndex(void) const;

    /**
     * get total number of atoms for types
     *
     */
    const Vec1Int& get_nAtomsType(void) const;

    /**
     * get total number of atoms for certain type
     *@param type index of type for which atom number should be obtained
     */
    const Int& get_nAtomsType(const std::size_t type) const;

    /// get total number of atoms
    const std::size_t& get_nAtoms(void) const;

    /**
     * get all nearest neighbor data from one atom pair
     *
     * call function as auto &[ type, r, x, y, z, x_norm, y_norm, z_norm ] = get_neighborData(
     * atomIndx, neighIndx )
     * @param atomIndex integer defining the central atom
     * @param indxNeigh integer defining entry of neighbor atom
     */
    const std::tuple<const Int&,
                     const Real&,
                     const Real&,
                     const Real&,
                     const Real&,
                     const Real&,
                     const Real&,
                     const Real&>
    get_neighborData(const std::size_t atomIndex, const std::size_t indxNeigh) const;

    /**
     * return the number of atoms in current neighbor list
     */
    std::size_t get_numberAtoms(void) const;

    /**
     *
     * get the number of different types in current structure
     */
    std::size_t get_numberTypes(void) const;

    /**
     * getter of control variable if neighbor list is type sorted
     */
    bool is_typeSorted(void) const;

    /**
     * get number of atoms per type
     */
    std::vector<std::size_t> get_numberAtomsPerType(void) const;

    /*******************************************************************************************
     * get array which stores the index of a central atom per type
     *******************************************************************************************/
    const Vec1Int& get_centralAtomIndexPerType(void) const;

    /**
     * get unique atom types
     *
     *
     */
    const Vec1Int& get_uniqueTypes(void) const;

    /**
     *returns the latticeVolume
     */
    const Real& get_latticeVolume(void) const;
    /// Get cutoff radius used for this neighbor list.
    Real get_cutOff() const;
    /*******************************************************************************************
     * setting the total number of atoms central atoms in the neighbor.
     *
     * @param nAtoms number of central atoms in current neighbor list
     *
     * @note this variable has to be set when the neighbor list is filled via an interface
     * for example to VASP. If the neighbor list is computed with it's own routines
     * the variable will be automatically set.
     *******************************************************************************************/
    void set_nAtoms(const Int nAtoms);
    /*******************************************************************************************
     * setting the total number of atom types in neighbor list
     *
     * @param nTypes number of unique atom types in neighbor list
     *
     * @note this variable has to be set when the neighbor list is filled via an interface
     * for example to VASP. If the neighbor list is computed with it's own routines
     * the variable will be automatically set.
     *******************************************************************************************/
    void set_nTypes(const Int nTypes);
    /*******************************************************************************************
     * setting up the centralAtomIndexPerType array
     *
     * Array stores for given global atom index the number of the atom within it's type
     *******************************************************************************************/
    void compute_centralAtomIndexPerType(void);

  private:
    /// check for unique atom types
    void compute_unique_types(const Vec1Int& types);

    /** compute nearest neighbor data for a single atom starting from
     *  direct coordinates, results will be given in Cartesian coordinates
     *
     * @param indx  unsigned integer denoting the actual central atom
     * @param positions structure containing containing information about atomic structure
     */
    void computeNearestNeighborsSingleAtomDirect(const std::size_t indx,
                                                 const Structure&  positions);
    /** compute nearest neighbor data for a single atom starting from
     *  Cartesian coordinates, results will be given in Cartesian coordinates
     *
     * @param indx  unsigned integer denoting the actual central atom
     * @param positions structure containing containing information about atomic structure
     */
    void computeNearestNeighborsSingleAtomCartesian(const std::size_t indx,
                                                    const Structure&  positions);
    /** allocate work arrays for computation of nearest neighbors single atom
     *
     * @param nn_list      temporary nearest (N) neighbor list for single atom
     * @param types_index  temporary integer (N) array containing atom types associated to nn_list
     * @param dist         temporary real (N) array containing distances in same order as nn_list
     * @param dist_vecs    temporary real (3*N) array containint distance vectors in same order as
     * nn_list first dimension stores x,y,z
     * @param normed_dist  temporary real (3*N) array containint nromalized distance vectors
     *                     in same order as nn_list first dimension stores x,y,z
     */
    void allocateWorkArrays(Vec1Int&    nn_list,
                            Vec1Int&    types_index,
                            Vec1Real&   dist,
                            Vec1Real&   dist_vecs,
                            Vec1Real&   normed_dist,
                            std::size_t N);

    /** fill the nearest neighbor arrays as defined at VARIABLES of this object
     *  for every atom; index defines the central atom
     *
     * @param index        size_t storing the index of the actual central atom
     * @param nelements    size_t denoting how many nearest neighbors have been found
     * @param nn_list      temporary nearest (N) neighbor list for single atom
     * @param types_index  temporary integer (N) array containing atom types associated to nn_list
     * @param dist         temporary real (N) array containing distances in same order as nn_list
     * @param dist_vecs    temporary real (3*N) array containint distance vectors in same order as
     * nn_list first dimension stores x,y,z
     * @param normed_dist  temporary real (3*N) array containint nromalized distance vectors
     *                     in same order as nn_list first dimension stores x,y,z
     */
    void fillNeighborArrays(const std::size_t index,
                            const std::size_t nelements,
                            const Vec1Int&    nn_list,
                            const Vec1Int&    types_index,
                            const Vec1Real&   dist,
                            const Vec1Real&   dist_vecs,
                            const Vec1Real&   normed_dist);
    /**
     * compute how many periodic boxes have to be searched for
     * nearest neighbors of the given central atom
     *
     * @param pos 3d-vector containing Cartesian coordinates of actual central atom
     * @param inverse_lattice inverse of the supplied bravais lattice
     * @param xmin boxes in negative x,y,z direction that have to be chcked for nearest neighbors
     * @param xmax boxes in positive x,y,z direction that have to be chcked for nearest neighbors
     */
    void set_periodic_size_images(const Vec1Real& pos_central,
                                  const Lattice&  inverse_lattice,
                                  Vec1Int&        xmin,
                                  Vec1Int&        xmax);

    /**
     * allocate storage arrays
     *@param size primary size of nearest neighbor array
     * number of central atoms in nearest neighbor list
     */
    void allocateStorageArrays(std::size_t size);

    // VARIABLES
    /// nearest neighbor array containing index to position array
    ShVec2Int globalIndex;
    /// Type index of neighbor atom
    ShVec2Int typeIndex;
    /// Type index of central atom
    ShVec1Int typeIndexCentral;
    /// nearest neighbor distances
    ShVec2Real distances;
    /// nearest neighbor connection vectors
    ShVec2Real connectionVector;
    /// nearest neighbor normalized connection vectors
    ShVec2Real connectionVectorNormalized;
    /// number of nearest neighbors per central atom
    ShVec1Int numberNeighbors;
    /// number of nearest neighbors per type for central atom
    ShVec2Int numberNeighborsType;
    /// total number of atoms per type
    ShVec1Int nAtomsType;
    /*******************************************************************************************
     * gives for for global central atom index the atom number within a type
     *******************************************************************************************/
    ShVec1Int centralAtomIndexPerType;

    // parameters needed for computation
    /// cutoff radius in units of input coordinates
    Real cutOff;
    /// squared cutoff radius
    Real cutOffSquared;
    /// minimum distance atoms have to be apart
    /// to avoid self neighbors
    Real eps;

    /// point to array start indices for every atom type / both start and end point to zero
    /// if no atoms of type are present
    ShVec2Int typeStart;
    /// point to array end indices for every atom type
    ShVec2Int typeEnd;
    /// point to array start indices for every atom type
    ShVec2Int neighborType;

    /// unique types in supplied structure
    std::map<Int, Int> unique_types;
    /// unique types in supplied structure
    std::map<Int, Int> types_unique;
    /// number of different types in supplied atoms structure
    Int nTypes;

    /// total number of atoms
    std::size_t nAtoms;

    /// decides the storage order of the nearest neighbor array
    /// sorts nearest neighbors according to types
    bool isTypeSorted;
    /// decides the storage order of the nearest neighbor array
    /// sorts nearest neighbors according to distances
    bool isDistSorted;
    /**
     * volume of box of atoms from which neighbor list is computed
     *
     *the volume is needed for the stress tensor compuation in
     *DescriptorSHS2
     */
    Real latticeVolume;
};

} //namespace vaspml

#endif
