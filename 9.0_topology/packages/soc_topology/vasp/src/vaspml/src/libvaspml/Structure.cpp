#include "Structure.hpp"

#include "constants.hpp"
#include "math.hpp"
#include "utils.hpp"

#include <cmath>
#include <iostream>

using namespace vaspml;

Structure::Structure(const ShVec1Real&   lattice,
                     const ShVec1Real&   positions,
                     const ShVec1Int&    types,
                     const ShVec1String& typeNames,
                     const ShVec1Int&    numAtomsPerType) :
    lattice(lattice),
    positions(assignOrMakeShared(positions)),
    typeNames(assignOrMakeShared(typeNames)),
    types(assignOrMakeShared(types)),
    numAtomsPerType(assignOrMakeShared(numAtomsPerType))
{}

void Structure::readPoscar(const std::string& fname)
{

    std::ifstream infile_stream = file_io::openFileI(fname);

    std::string dummy;
    // read header from POSCAR file
    getline(infile_stream, dummy);
    read_lattice(infile_stream);
    read_types(infile_stream);
    getline(infile_stream, dummy);

    direct = string_tools::checkForSubstring(string_tools::makeLowerCase(dummy), "d");
    read_positions(infile_stream);
    if (direct) { shift_atoms_primitive_cell(); }
    else
    {
        cartesianToDirect();
        shift_atoms_primitive_cell();
        directToCartesian();
    }

    build_types();
    infile_stream.close();
}

void Structure::set_lattice(const ShVec1Real& lattice_in)
{
    lattice.set_components((*lattice_in)[0],
                           (*lattice_in)[1],
                           (*lattice_in)[2],
                           (*lattice_in)[3],
                           (*lattice_in)[4],
                           (*lattice_in)[5],
                           (*lattice_in)[6],
                           (*lattice_in)[7],
                           (*lattice_in)[8]);
    inverseLattice = lattice.computeInverse();
}

void Structure::set_positions(const ShVec1Real& positions, const bool isDirect_in)
{

    direct = isDirect_in;

    this->positions = positions;
    if (direct) { shift_atoms_primitive_cell(); }
    else
    {
        cartesianToDirect();
        shift_atoms_primitive_cell();
        directToCartesian();
    }
}

void Structure::set_types(const ShVec1String& atom_types_in, const ShVec1Int& number_atoms_in)
{
    typeNames = atom_types_in;
    numAtomsPerType = number_atoms_in;
    numAtoms = math::sumVector(*number_atoms_in);
    build_types();
}

void Structure::importPoscarData(const ShVec1Real&   lattice_in,
                                 const ShVec1String& atom_types_in,
                                 const ShVec1Int&    number_atoms_in,
                                 const ShVec1Real&   positions,
                                 const bool          isDirect_in)
{

    set_lattice(lattice_in);
    set_types(atom_types_in, number_atoms_in);
    set_positions(positions, isDirect_in);
}

void Structure::directToCartesian(void)
{

    if (direct) { rescale_coordinates(lattice); }
    else
    {
        std::cout << "Warning in Structure::position_direct_to_cartesian" << std::endl;
        std::cout << "You already have direct coordinates. Doing nothing" << std::endl;
    }
    direct = false;
}

void Structure::cartesianToDirect(void)
{

    if (!direct) { rescale_coordinates(inverseLattice); }
    else
    {
        std::cout << "Warning in Structure::position_cartesian_to_direct" << std::endl;
        std::cout << "You already have direct coordinates. Doing nothing" << std::endl;
    }
    direct = true;
}

const std::tuple<const Real&, const Real&, const Real&> Structure::getAtom(
    const size_t atomIndx) const
{
    return std::tie((*positions)[3 * atomIndx],
                    (*positions)[3 * atomIndx + 1],
                    (*positions)[3 * atomIndx + 2]);
}

bool Structure::isDirect(void) const
{
    return direct;
}

const ShVec1Real& Structure::get_positions(void) const
{
    return positions;
}
const Lattice& Structure::get_lattice(void) const
{
    return lattice;
}
const Lattice& Structure::get_inverseLattice(void) const
{
    return inverseLattice;
}
const ShVec1Int& Structure::get_types(void) const
{
    return types;
}
const ShVec1String& Structure::get_typeNames(void) const
{
    return typeNames;
}
const ShVec1Int& Structure::get_numAtomsPerType(void) const
{
    return numAtomsPerType;
}

size_t Structure::get_numAtoms(void) const
{
    return numAtoms;
}

void Structure::writeToScreen() const
{
    std::cout << "-------------------------------------------\n";
    std::cout << "WARNING: All output uses Bohr length units.\n";
    std::cout << "-------------------------------------------\n";
    std::cout << str("Scaling factor: %24.16E\n", scale);
    std::cout << "Lattice:\n";
    lattice.writeToScreen();
    std::cout << "Inverse lattice:\n";
    inverseLattice.writeToScreen();
    // std::cout << "Inverse inverse lattice:\n";
    // Lattice invinv_lattice = inverse_lattice.computeInverse();
    // invinv_lattice.writeToScreen();
    std::cout << "(lattice) * (inverse lattice):\n";
    for (std::size_t i = 0; i < 3; ++i)
    {
        for (std::size_t j = 0; j < 3; ++j)
        {
            Real mij = 0.0;
            for (std::size_t k = 0; k < 3; ++k)
            {
                mij += lattice.component(i, k) * inverseLattice.component(k, j);
            }
            std::cout << str(" %10.2E", mij);
        }
        std::cout << "\n";
    }
    std::cout << str("Number of atoms: %zu\n", numAtoms);
    std::cout << "Atom types:";
    for (auto t : *typeNames) std::cout << " " + t;
    std::cout << "\n";
    std::cout << "Number of atoms per type:";
    for (auto n : *numAtomsPerType) std::cout << str(" %d", n);
    std::cout << "\n";
    if (direct) std::cout << "Direct\n";
    else std::cout << "Cartesian\n";
    for (std::size_t i = 0; i < numAtoms; ++i)
    {
        std::cout << str("%24.16E %24.16E %24.16E %d %s\n",
                         (*positions)[i],
                         (*positions)[i + 1],
                         (*positions)[i + 2],
                         (*types)[i],
                         (*typeNames)[(*types)[i]].c_str());
    }

    return;
}

//**************************************************************************************************
// Private member functions
//**************************************************************************************************

void Structure::read_lattice(std::ifstream& infile_stream)
{

    std::string line;
    getline(infile_stream, line);
    line = string_tools::trim(line);
    Vec1Real scaleVector = string_tools::extractData<Real>(line);
    if (scaleVector.size() != 1)
    {
        global_scope::tutor.error("Unsupported scaling factor format encountered while reading "
                                  "POSCAR, only a single value is supported.");
    }
    scale = scaleVector.at(0);
    Vec1Real temp_lattice;

    for (auto i = 0; i < 3; i++)
    {
        getline(infile_stream, line);
        Vec1Real tmp = string_tools::extractData<Real>(line);
        math::vectorTimesScalarNoCopy(tmp, scale / constants::AUTOA);
        temp_lattice.insert(temp_lattice.end(), tmp.begin(), tmp.end());
    }

    lattice.set_components(temp_lattice);

    inverseLattice = lattice.computeInverse();
}

void Structure::read_positions(std::ifstream& infile_stream)
{

    numAtoms = math::sumVector(*numAtomsPerType);
    positions->resize(3 * numAtoms);
    size_t counter = 0;
    for (size_t i = 0; i < numAtoms; i++)
    {
        std::string line;
        getline(infile_stream, line);
        Vec1Real atom = string_tools::extractData<Real>(line);
        if (!direct) math::vectorTimesScalarNoCopy(atom, scale / constants::AUTOA);
        (*positions)[counter] = atom[0];
        counter++;
        (*positions)[counter] = atom[1];
        counter++;
        (*positions)[counter] = atom[2];
        counter++;
    }
}

void Structure::read_types(std::ifstream& infile_stream)
{

    std::string line;
    getline(infile_stream, line);
    *typeNames = string_tools::extractData<std::string>(line);
    getline(infile_stream, line);
    *numAtomsPerType = string_tools::extractData<Int>(line);
}

void Structure::shift_atoms_primitive_cell(void)
{

    Real one = (Real)1;
    if (direct)
    {
        for (size_t i = 0; i < positions->size(); i++)
        {
            (*positions)[i] = std::modf((*positions)[i] + (Real)100, &one);
        }
    }
    else
    {
        std::cout << "Error in Structure::shift_atoms_primitive_cell " << std::endl;
        std::cout << "Function only works on direct coordinates" << std::endl;
    }
}

void Structure::rescale_coordinates(const Lattice& lattice)
{

    for (size_t i = 0; i < 3 * numAtoms; i += 3)
    {
        lattice.timesVectorInPlace((*positions)[i], (*positions)[i + 1], (*positions)[i + 2]);
    }
}

void Structure::build_types(void)
{

    types->resize(numAtoms);
    size_t indx = 0;
    for (size_t i = 0; i < numAtomsPerType->size(); i++)
    {
        for (size_t j = 0; j < (size_t)(*numAtomsPerType)[i]; j++)
        {
            (*types)[indx] = i;
            indx++;
        }
    }
}
