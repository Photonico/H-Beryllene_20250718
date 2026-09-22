#ifndef STRUCTURE_HPP
#define STRUCTURE_HPP

#include "Lattice.hpp"
#include "types.hpp"

#include <cstddef>
#include <fstream>
#include <string>
#include <tuple>

namespace vaspml
{

class Structure
{
  public:
    /***********************************************************************************************
     * Constructor, optionally pass smart pointers to existing memory.
     **********************************************************************************************/
    Structure(const ShVec1Real&   lattice = nullptr,
              const ShVec1Real&   positions = nullptr,
              const ShVec1Int&    types = nullptr,
              const ShVec1String& typeNames = nullptr,
              const ShVec1Int&    numAtomsPerType = nullptr);
    /***********************************************************************************************
     * Read in structure data from POSCAR file.
     *
     * @param fname Name of input file in POSCAR format.
     **********************************************************************************************/
    void readPoscar(const String& fname);
    /***********************************************************************************************
     * Set lattice vectors and compute inverse lattice matrix.
     *
     * @param lattice_in   3x3 matrix in linearized form as xx,xy,xz,yx,yy,yz,zx,zy,zz
     **********************************************************************************************/
    void set_lattice(const ShVec1Real& lattice_in);
    /* set positions and shift them to superell
     * @param positions    -> 3*N array containing the atomic positions in either direct or
     * cartesian coordinates
     * @param isDirect_in  -> specify in which type the coordinates are supplied; default is direct
     * coordinates
     *
     * @note positions is supplied as 3*N array in format x1,y1,z1,x2,y2,z2 where the number
     * refers to the atom index
     */
    void set_positions(const ShVec1Real& positions, const bool isDirect_in = true);
    /* set atom types and atom numbers; set atom_type_index array
     *
     * @param atom_types_in      Vector containing string determining the atom types
     * @param number_atoms_in    integer array of same length as atom types containing number of
     * atoms per type
     */
    void set_types(const ShVec1String& atom_types_in, const ShVec1Int& number_atoms_in);

    /* set the whole poscar structure at once
     * @param lattice_in   Lattice to which lattice has to be set, format xx,xy,xz,yx,yy,yz,zx,zy,zz
     * @param positions    -> 3*N array containing the atomic positions in either direct or
     * cartesian coordinates
     * @param isDirect_in  -> specify in which type the coordinates are supplied; default is direct
     * coordinates
     * @param atom_types_in    ->  Vector containing string determining the atom types
     * @param number_atoms_in  ->  integer array of same length as atom types containing number of
     * atoms per type
     *
     * @note positions is supplied as 3*N array in format x1,y1,z1,x2,y2,z2 where the number
     * refers to the atom index
     */
    void importPoscarData(const ShVec1Real&   lattice_in,
                          const ShVec1String& atom_types_in,
                          const ShVec1Int&    number_atoms_in,
                          const ShVec1Real&   positions,
                          const bool          isDirect_in = true);

    /// make cartesian positions from direct
    void directToCartesian(void);
    /// make direct positions from cartesian
    void cartesianToDirect(void);
    /// check if coordinates are direct or cartesian
    bool isDirect(void) const;

    /**
     * get xyz coordinates of atom
     *
     * @param atomIndx index of atom in positions array
     *
     * @return xyz as const tuple
     */
    const std::tuple<const Real&, const Real&, const Real&> getAtom(
        const std::size_t atomIndx) const;

    /**
     * get xyz coordinates of atom
     *
     * @param atomIndx index of atom in positions array
     *
     * @return xyz as rvalue refernce tuple
     */
    std::tuple<Real&&, Real&&, Real&&> getAtom(const std::size_t atomIndx);

    // getters
    /* get positions as a const reference
     * for max performance call with const Matrix<Real>& positions = get_positions()
     * to not call copy constructor
     *
     * @return positions 3*N positions array
     *
     * @note format is in x1,y1,z1,x2,y2,z2...
     */
    const ShVec1Real& get_positions(void) const;
    /* get lattice as a const reference
     * for max performance call with const Lattice& lattice = get_lattice()
     * to not call copy constructor
     */
    const Lattice& get_lattice(void) const;
    /* get inverse lattice as a const reference
     */
    const Lattice& get_inverseLattice(void) const;
    /* get atom_types as Integer indices
     */
    const ShVec1Int& get_types(void) const;
    /* get atom_types as strings
     */
    const ShVec1String& get_typeNames(void) const;
    /* get total number of atom in current structure
     */
    const ShVec1Int& get_numAtomsPerType(void) const;
    /// get the total number of atoms
    std::size_t get_numAtoms(void) const;
    /// Write structure information to screen.
    void writeToScreen() const;

  private:
    /// read lattice from poscar
    void read_lattice(std::ifstream& infile_stream);
    /// read positions and types from poscar file
    void read_positions(std::ifstream& infile_stream);
    /// read atom types and numbers
    void read_types(std::ifstream& infile_stream);
    /// shift atoms to primitive cell
    void shift_atoms_primitive_cell(void);
    /// direct to cartesian for supplied array
    void rescale_coordinates(const Lattice& lattice);
    /**
     * set up an integer array for atom types
     * types can be used to point to
     */
    void build_types(void);

    /// check if positions are in cartesian or direct coordinates
    bool direct;
    /// Lattice vectors of given structure
    Lattice lattice;
    /// inverse lattice vectors of given structure
    Lattice inverseLattice;
    /// number of atoms
    std::size_t numAtoms;
    /**
     * positions stored in format x1,y1,z1,x2,y2,z2,x3,y3,z3...
     */
    ShVec1Real positions;
    /// stores atom types
    ShVec1String typeNames;
    /// type index starting at 0 and going up to N number atom types
    ShVec1Int types;
    /// number of atoms per type
    ShVec1Int numAtomsPerType;
    /// Scaling factor.
    Real scale;
};

} // namespace vaspml

#endif
