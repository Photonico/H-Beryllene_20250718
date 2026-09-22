#ifndef INTERFACEVASP_HPP
#define INTERFACEVASP_HPP

#include "BasisFunctions.hpp"
#include "Frame.hpp"
#include "KernelPolynomial.hpp"
#include "MlMPI.hpp"
#include "Predictor.hpp"
#include "TypeMap.hpp"
#include "types.hpp"

#include <map>
#include <memory>
#include <string>

namespace vaspml
{

/***************************************************************************************************
 * Specialized Frame class for the VASP interface.
 *
 * This struct is intended to be created from VASP via Fortran to C function calls, see
 * createFrame(), destroyFrame(), etc.
 **************************************************************************************************/
struct InterfaceVASP : public Frame
{

    /*******************************************************************************************
     * manages the input of the ML_FF file. Stores and distributes the input data
     * needed for the force field computation
     *******************************************************************************************/
    IoHandlerML_FF forceFieldData;
    /*******************************************************************************************
     * main mpi class for vaspml
     *
     * all other mpi communicators have to be created from this instance.
     * From this instance we can create shared memory mpi communicators
     *******************************************************************************************/
    std::shared_ptr<MlMPI> mpiMain = nullptr;
    /*******************************************************************************************
     * stores the basis functions as spherical harmonics and modified spherical Bessel functions.
     *
     * This data structure handles the basis set functions for DescriptorSHS2 and DescriptorSHS3
     *******************************************************************************************/
    std::map<std::string, std::shared_ptr<BasisFunctions>> basisFunctions;
    /*******************************************************************************************
     * class takes the descriptors DescriptorSHS2 and DescriptorSHS3 as input
     * and computes the kernel matrix and derivatives of the specific kernel matrix.
     *
     * The kernel matrix is computed as a similarity meassure between local reference
     * configurations and local configurations of the current structure.
     * A local reference configuration is a local representation ( within some cutoff ) of the
     * 2 and 3-body interaction in terms of the classes DescriptorSHS2 and DescriptorSHS3
     *******************************************************************************************/
    KernelPolynomial polynomialKernel;
    /*******************************************************************************************
     * class takes the polyonimalKernel as an input and computes energy, forces and stress tensor
     *
     * The class uses the polynomialKernel and the fitting weights obtained from the IoHandlerML_FF
     * to compute the total energy of the system
     * The forces are computed from the derivative of the polyonomialKernel.
     * The stress tensor is obtained from the computed pairForces with the help
     * of the virial theorem
     *******************************************************************************************/
    Predictor predictor;
    /*******************************************************************************************
     * this class manages the correspondence between atom types in the machine learning force field
     * and the structure which is currently analysed
     *
     * The types in the structure always have to be a subset of the types present in the machine
     * learning force field. Otherwise, it is not possible to make any predictions.
     * Because this would indicate there are types in the structure not considered during training
     *******************************************************************************************/
    std::shared_ptr<TypeMap> typeMapPtr;
    /*******************************************************************************************
     * allocating and resizing neighbor list arrays
     *
     * @param nions number of ions
     * @param key   descriptor key can either be 2-body or 3-body
     *******************************************************************************************/
    void resizeNeighborArrays(const Int nions, const String& key);
    /*******************************************************************************************
     * setting the number of atoms per type in the neighbor list array
     *
     * @param ntypes number of different atom types
     * @param nAtomsType number of atoms per type
     * @param keyIn the key defining to which neighbor list data is written. can be 2-body or 3-body
     *
     * will write to 2-body-n_nAtomsType and 3-body-n_nAtomsType
     *
     *******************************************************************************************/
    void set_nAtomsType(const Int ntypes, const Int* nAtomsType, const String& keyIn);
    /*******************************************************************************************
     * filling the neighbor array for one specific central atom defined by atomNumber
     *
     * @param numberNeighbors number of neighbor atoms withi the cutoff radius
     * @param atomNumber is the index of the central atom in the neighbor list arrays
     * @param centralType atom type of the central atom around which the neighbors are computed
     * @param neighborIndex index of neighbor atom in the POSION structure of vasp -1
     * @param neighborTypes atom type of neighbor atom coded as an integer
     * @param neighborDist distances between neighbor atoms and current central atom
     * @param neighborConnect normalized connection vector between central atom and neighbor atom
     * @param key key which describes to which descriptor the neighbor list belongs. Can either be
     * 2-body or 3-body
     *
     * @note the array neighborConnect is stored as neighbor 1 x neighbor 1 y, neighbor 1 z,
     * neighbor 2 x, neighbor 2 y, neighbor 3 z,...
     *******************************************************************************************/
    void fillNeighhborArrays(const Int          numberNeighbors,
                             const Int          atomNumber,
                             const Int          centralType,
                             const Int*         neighborIndex,
                             const Int*         neighborTypes,
                             const Real*        neighborDist,
                             const Real*        neighborConnect,
                             const std::string& key);
    /*******************************************************************************************
     * setting up the type map which gives the relationship between force field types and structure
     * types
     *
     * @param types string array which contains the atom types
     *******************************************************************************************/
    void set_typeMap(const Vec1String& types);
    /*******************************************************************************************
     * compute the total ionic force from the central and pair derivatives of the
     * polynomialKernel matrix computed in the predictor class
     *******************************************************************************************/
    void computeForces(void** ptrFrame);
};

} // namespace vaspml

extern "C"
{
using namespace vaspml;
/***************************************************************************************************
 * Create an instance of InterfaceVASP and return void pointer to it.
 *
 * @return Pointer to newly created InterfaceVASP instance.
 *
 * The pointer is received and stored on the Fortran side to recall the instance.
 **************************************************************************************************/
void* createFrame();
/***************************************************************************************************
 * Destroy an instance of InterfaceVASP.
 *
 * @param ptrFrame Pointer to InterfaceVASP (passed on with call by reference).
 *
 * @pre The pointer has to be created beforehand via createFrame().
 **************************************************************************************************/
void destroyFrame(void** ptrFrame);
/***************************************************************************************************
 * Set up a force field from information in ML_FF file.
 *
 * @param ptrFrame Pointer to InterfaceVASP (passed on with call by reference).
 * @param mpiComm mpi communicator which will be used as COMM_WORLD by vaspml
 *
 * @pre The pointer has to be created beforehand via createFrame().
 **************************************************************************************************/
void setupForceField(void** ptrFrame, MPI_Fint* mpiComm = nullptr);
/*******************************************************************************************
 * getting weighting factor of DescriptorSHS2
 *
 * @param ptrFrame pointer to frame class will be cast to InterfaceVASP
 *******************************************************************************************/
Real get_W1(void** ptrFrame);
/*******************************************************************************************
 * getting weighting factor of DescriptorSHS3
 *
 * @param ptrFrame pointer to frame class will be cast to InterfaceVASP
 *******************************************************************************************/
Real get_W2(void** ptrFrame);
/*******************************************************************************************
 * getting cutoff radius of DescriptorSHS2
 *
 * @param ptrFrame pointer to frame class will be cast to InterfaceVASP
 *******************************************************************************************/
Real get_RCUT1(void** ptrFrame);
/*******************************************************************************************
 * getting cutoff radius of DescriptorSHS3
 *
 * @param ptrFrame pointer to frame class will be cast to InterfaceVASP
 *******************************************************************************************/
Real get_RCUT2(void** ptrFrame);
/*******************************************************************************************
 * resize the neighbor list arrays
 *
 * @param ptrFrame pointer to frame class will be cast to InterfaceVASP
 * @param nions number of ions
 * @param key   descriptor key can either be 2-body or 3-body
 *
 * calls resizeNeighborArrays from InterfaceVASP structure
 *******************************************************************************************/
void resizeNeighborArrays(void** ptrFrame, const Int* nions, const char* keyIn);
/*******************************************************************************************
 * set the number of atoms in the neighorlist array
 *
 * @param ptrFrame pointer to frame class will be cast to InterfaceVASP
 * @param ntypes number of different atom types
 * @param nAtomsType number of atoms per type
 * @param keyIn the key defining to which neighbor list data is written. can be 2-body or 3-body
 *
 * calls the function set_nAtomsType from the InterfaceVASP structure
 *******************************************************************************************/
void set_nAtomsType(void** ptrFrame, const Int* ntypes, const Int* nAtomsType, const char* keyIn);
/*******************************************************************************************
 * filling the neighbor array for a single function
 *
 * @param ptrFrame pointer to frame class will be cast to InterfaceVASP
 * @param numberNeighbors number of neighbor atoms withi the cutoff radius
 * @param atomNumber is the index of the central atom in the neighbor list arrays
 * @param centralType atom type of the central atom around which the neighbors are computed
 * @param neighborIndex index of neighbor atom in the POSION structure of vasp -1
 * @param neighborTypes atom type of neighbor atom coded as an integer
 * @param neighborDist distances between neighbor atoms and current central atom
 * @param neighborConnect normalized connection vector between central atom and neighbor atom
 * @param key key which describes to which descriptor the neighbor list belongs. Can either be
 * 2-body or 3-body
 *
 * function will call the fillNeighhborArrays arrays of the InterfaceVASP structure
 *******************************************************************************************/
void fillNeighhborArrays(void**      ptrFrame,
                         const Int*  numberNeighbors,
                         const Int*  atomNumber,
                         const Int*  centralType,
                         const Int*  neighborIndex,
                         const Int*  neighborTypes,
                         const Real* neighborDist,
                         const Real* neighborConnect,
                         const char* keyIn);
/*******************************************************************************************
 * setting up type map to give relationship between force field and structure types
 * @param ptrFrame pointer to frame class will be cast to InterfaceVASP
 * @param types string array which contains the atom types
 *
 * calls set_typeMap from InterfaceVASP structure
 *******************************************************************************************/
void set_typeMap(void** ptrFrame, const char* types);
/*******************************************************************************************
 * compute the total ionic force given by the machine leraning model on a single central atom
 * has to be allreduced on fortran side when using several MPI ranks
 *
 * @param ptrFrame pointer to frame class will be cast to InterfaceVASP
 * @param centralVASP indicates the central atom index in the TIFOR array from vasp
 * @param element gives the entry of the atom in the force array on the c++ side of the interface
 * @param force array where the force is stored to.
 *******************************************************************************************/
void fillForceSingleAtom(void**     ptrFrame,
                         const Int* nions,
                         const Int* centralVasp,
                         const Int* element,
                         Real*      force);
/*******************************************************************************************
 * compute the stress tensor with the virial theorem
 *
 * @param ptrFrame pointer to frame class will be cast to InterfaceVASP
 * @param volume is the volume of the considered simulation box
 *******************************************************************************************/
void update(void** ptrFrame, const Real* volume);
/*******************************************************************************************
 * retrieve the total energy from the system. Note this energy has to be allreduced on fortran side
 * when using several mpi ranks
 *
 * @param ptrFrame pointer to frame class will be cast to InterfaceVASP
 * @param totalEnergy on output array stores the total energy
 *******************************************************************************************/
void getPotentialEnergy(void** ptrFrame, Real* totalEnergy);
/*******************************************************************************************
 * getter for the total stress tensor. Has to be allreduced on fortran side when
 * several MPI ranks are used
 * @param ptrFrame pointer to frame class will be cast to InterfaceVASP
 * @param stressTensor on output this will be a flattened 3x3 array storing the stress tensor
 *******************************************************************************************/
void getStressTensor(void** ptrFrame, Real* stressTensor);

} // extern "C"

#endif
