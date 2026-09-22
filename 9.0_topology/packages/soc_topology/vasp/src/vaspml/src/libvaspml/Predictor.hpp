#ifndef PREDICTOR_HPP
#define PREDICTOR_HPP

#include "ArrayResizing.hpp"
#include "DescriptorCollector.hpp"
#include "IoHandlerML_FF.hpp"
#include "Kernel.hpp"
#include "Linalg.hpp"
#include "MlMPI.hpp"
#include "ParallelEnvironemt.hpp"
#include "ShmemArray.hpp"
#include "SmartEnum.hpp"
#include "TypeMap.hpp"
#include "types.hpp"

#include <cstddef>
#include <iostream>
#include <map>
#include <memory>
#include <variant>

namespace vaspml
{

std::map<String, std::shared_ptr<ShmemArray2DVariableLen<Real>>> make_descriptorsRefConfsWeighted(
    const std::map<String, ShVec1Int>& featureSpaceSize,
    const std::map<String, Real>&      weights,
    const std::shared_ptr<MlMPI>&      mpiIn);

std::map<String, ShVec2Real> makeMaps2D(void);

std::map<String, ShVec1Real>                   makeMaps1D(void);
std::map<String, Vec1Size_t>                   makeMaps1DSize_t(void);
std::map<String, ArrayResizing1D>              makeResizeMap1D(void);
std::map<String, ArrayResizing2D>              makeResizeMap2D(void);
std::shared_ptr<ShmemArray2DVariableLen<Real>> make_regressionCoefficients(
    const IoHandlerML_FF&         inputParameters,
    const std::shared_ptr<MlMPI>& mpiIn);

enum class TotalEnergyType
{
    IsolatedAtom,
    AverageTrainEnergy,
};

/*******************************************************************************************
 * functor for computing the total atomic force acting on an ion
 *
 * Create by SumFunctorForce functor( atomicForces, pairForces[ centralAtom ] ); and use
 * with a std algorithm as
 * std::for_each( neighborIndex[ centralAtom ].cbegin(), neighborIndex[ centralAtom ].cend(),
 *functor ); This procedure has to be called from an loop over central atoms. One central atom per
 *time can be processed with this routines.
 *******************************************************************************************/
class SumFunctorForce
{
  public:
    SumFunctorForce(Vec1Real& atomicForce, const Vec1Real& pairForce) :
        atomicForce(atomicForce),
        pairForce(pairForce),
        counter(0)
    {}
    void operator()(const Int& neighborIndex)
    {
        atomicForce[3 * neighborIndex] += pairForce[3 * counter];
        atomicForce[3 * neighborIndex + 1] += pairForce[3 * counter + 1];
        atomicForce[3 * neighborIndex + 2] += pairForce[3 * counter + 2];
        counter++;
    }

  private:
    /*******************************************************************************************
     * total atomic force.
     *
     * Array of size 3*number atoms. The storage order is atom1 x atom1 y atom1 z,
     * atom2 x, atom2 y, atom3 z... atomicForce array will always point to the atomicForce
     * memory of the components of Predictor::atomicForce
     *******************************************************************************************/
    Vec1Real& atomicForce;
    /*******************************************************************************************
     * pair atomic force.
     *
     * Array of size 3*number atoms. The storage order is atom1 x atom1 y atom1 z,
     * atom2 x, atom2 y, atom3 z... pairForces array will always point to the entries of pairForces
     * memory of the components of Predictor::pairForces
     *******************************************************************************************/
    const Vec1Real& pairForce;
    /*******************************************************************************************
     * counting the number of neighbor atoms which have been processed already
     *******************************************************************************************/
    Int counter = 0;
};
/*******************************************************************************************
 * functor computing the stress tensor contribution due to a single atom
 *
 * Create by
 * SumFunctorStress functor( stressTensor, pairForces[ centralAtom ],connectionVector[ centralAtom
 *], inverseVolume ); Then use by a std algorithm as std::for_each( neighborIndex[ centralAtom
 *].cbegin(),neighborIndex[ centralAtom ].cend(), functor ); This will compute the contribution of a
 *single central atom to the stress tensor. To get total stress repeat procedure for every central
 *atom.
 *******************************************************************************************/
class SumFunctorStress
{
  public:
    SumFunctorStress(Vec1Real&       stressTensor,
                     const Vec1Real& pairForces,
                     const Vec1Real& connectionVector,
                     const Real&     inverseVolume) :
        stressTensor(stressTensor),
        pairForces(pairForces),
        connectionVector(connectionVector),
        inverseVolume(inverseVolume),
        counter(0)
    {}
    void operator()(const std::size_t)
    {
        // compute xx term
        stressTensor[0] -= inverseVolume * pairForces[3 * counter] * connectionVector[3 * counter];
        // compute xy term
        stressTensor[1] -=
            inverseVolume * pairForces[3 * counter] * connectionVector[3 * counter + 1];
        // compute xz term
        stressTensor[2] -=
            inverseVolume * pairForces[3 * counter] * connectionVector[3 * counter + 2];
        // compute yy term
        stressTensor[4] -=
            inverseVolume * pairForces[3 * counter + 1] * connectionVector[3 * counter + 1];
        // compute yz term
        stressTensor[5] -=
            inverseVolume * pairForces[3 * counter + 1] * connectionVector[3 * counter + 2];
        // compute zz term
        stressTensor[8] -=
            inverseVolume * pairForces[3 * counter + 2] * connectionVector[3 * counter + 2];
        counter++;
    }

  private:
    Vec1Real& stressTensor;
    /*******************************************************************************************
     * Pair atomic force. Stores force between pair of atoms
     *
     * Array of size 3*number atoms. The storage order is atom1 x atom1 y atom1 z,
     * atom2 x, atom2 y, atom3 z... pairForces array will always point to the entries of pairForces
     * memory of the components of Predictor::pairForces
     *******************************************************************************************/
    const Vec1Real& pairForces;
    /*******************************************************************************************
     * connection vector array between atom pairs
     *
     * Array of size 3*number atoms. The storage order is atom1 x atom1 y atom1 z,
     * atom2 x, atom2 y, atom3 z... pairForces array will always point to the entries of
     * connectionVector memory of the components of NearestNeighborNSquare::connectionVector
     *******************************************************************************************/
    const Vec1Real& connectionVector;
    /*******************************************************************************************
     * inverse volume of the current box used for the simulation
     *******************************************************************************************/
    const Real& inverseVolume;
    /*******************************************************************************************
     * counting the number of neighbor atoms which have been processed already
     *******************************************************************************************/
    Int counter = 0;
};

/*******************************************************************************************
 * @class Predictor
 * @brief Compute energy, forces and stress tensor from some supplied Kernel
 *
 * Predictor is a class which stores the fitting weights (regression coeffcients) obtained
 * during Kernel training. When a Kernel is supplied to the Predictor class the programmer
 * will obtain energy, forces and stresses for the structure of which the supplied
 * Kernel was computed. The Predictor class is agnostic to the type of the supplied Kernel.
 * The class is also agnostic to the types of the Descriptors from which the Kernel
 * classes are made. To keep the agnosticism any Descriptors have to inherit from Descriptor
 * and any used Kernel has to inherit from the virtual Kernel base class. When done
 * in this way the Predictor class does not be touched to implement new Kernel or
 * Descriptor functionality.
 *******************************************************************************************/
class Predictor
{

  public:
    /*******************************************************************************************
     * default construction of Predictor class.
     *
     * @note class is not usable in this state, but can be transformed into
     * an usable state by calling Predictor::Init( const IoHandlerML_FF& inputParameters,
     * const Kernel& kernel, const std::shared_ptr<MlMPI>& mpiIn = nullptr,
     * const ShVec2Real& energyArray = nullptr, const ShVec1Real& atomicForces = nullptr,
     * const ShVec1Real& totalStressTensor = nullptr )
     *******************************************************************************************/
    Predictor(void) = default;
    /*******************************************************************************************
     * @param inputParameters stores force field parameters
     * @param mpiIn is used to store reference structure descriptors in shared memory
     * @param energyArray used to compute energy of structure
     * @param atomicForces used to store total force acting on ion
     * @param totalStressTensor used to store stress tensor summed over descriptor contributions
     *******************************************************************************************/
    Predictor(const IoHandlerML_FF&         inputParameters,
              const std::shared_ptr<MlMPI>& mpiIn = nullptr,
              const ShVec2Real& energyArray = nullptr,
              const ShVec1Real& atomicForces = nullptr,
              const ShVec1Real& totalStressTensor = nullptr,
              const ExecutionPolicy& algoExecution = ExecutionPolicy::cpuSingleCore);
    /*******************************************************************************************
     * transfer class from default constructed state into a state which can be used to predict
     * energies and forces
     *
     * @param inputParameters stores force field parameters
     * @param mpiIn is used to store reference structure descriptors in shared memory
     * @param energyArray used to compute energy of structure
     * @param atomicForces used to store total force acting on ion
     * @param totalStressTensor used to store stress tensor summed over descriptor contributions
     *******************************************************************************************/
    void init(const IoHandlerML_FF&         inputParameters,
              const std::shared_ptr<MlMPI>& mpiIn = nullptr,
              const ShVec2Real& energyArray = nullptr,
              const ShVec1Real& atomicForces = nullptr,
              const ShVec1Real& totalStressTensor = nullptr,
              const ExecutionPolicy& algoExecution = ExecutionPolicy::cpuSingleCore);
    /*******************************************************************************************
     * Compute total energy and forces for current structure supplied by kernel.
     *
     * For the forces only the arrays centralForces and pairForces will be computed.
     * From those the total ionic force can be computed with member function
     * compute_atomicForces( const DescriptorCollector& descriptorCollection );
     *
     * @param kernel is storing the similarity measure of the current structure
     * to the reference configurations
     * kernel also has to supply the derivatives ( derivativeMatrix ) of the kernel matrix
     *******************************************************************************************/
    void update(const Kernel& kernel);
    /**********************************************************************************
     *compute total energy of system and fill energyArray
     *
     @param kernel supplies the kernel matrix form which atomic energy is computed
     **********************************************************************************/
    void updateEnergy(const Kernel& kernel);
    /*******************************************************************************************
     * compute force arrays centralForces and pairForces
     *
     * @param kernel supplies the kernel matrix and kernel derivatives needed
     * for computation.
     * @note the kernel has to be updated outside of the current class
     *
     * @f[ F^{i}_{\alpha}=-\sum_{B}w_{B}\frac{\partial
     *K(\hat{\mathbf{X}}_{i},\hat{\mathbf{X}}_{B})}{\partial r^{i}_{\alpha}} @f] centralForces
     * @f[ F^{j}_{\alpha}=\sum_{B}\sum_{i=1}^{N_{atoms}}w_{B}\frac{\partial
     *K(\hat{\mathbf{X}}_{i},\hat{\mathbf{X}}_{B})}{\partial r^{i}_{\alpha}} @f] pairForces
     *******************************************************************************************/
    void updateForces(const Kernel& kernel);
    /*******************************************************************************************
     * Compute total forces acting on ions.
     *
     * Routine sums up the centralForces and the pairForces terms over the descriptors
     * to obtain the total force acting onto an ion. Depending on the supplied execution
     * policy the function will choose between CPU and GPU execution
     *
     *
     * @param descriptorCollection storing normalized descriptors. From there neighbor lists
     * can be retrieved to compute the total force acting on ion
     *******************************************************************************************/
    void compute_atomicForces(const DescriptorCollector& descriptorCollection);
    /*******************************************************************************************
     * return x component of pairForces
     *
     * @param centralAtom index of central Atom in neighbor list ( have same format )
     * @param neighborAtom index of neighbor atom for central atom centralAtom
     *******************************************************************************************/
    const Real& get_pairForcesX(const String&     key,
                                const std::size_t centralAtom,
                                const std::size_t neighborAtom) const;
    /*******************************************************************************************
     * return y component of pairForces
     *
     * @param centralAtom index of central Atom in neighbor list ( have same format )
     * @param neighborAtom index of neighbor atom for central atom centralAtom
     *******************************************************************************************/
    const Real& get_pairForcesY(const String&     key,
                                const std::size_t centralAtom,
                                const std::size_t neighborAtom) const;
    /*******************************************************************************************
     * return z component of pairForces
     *
     * @param centralAtom index of central Atom in neighbor list ( have same format )
     * @param neighborAtom index of neighbor atom for central atom centralAtom
     *******************************************************************************************/
    const Real& get_pairForcesZ(const String&     key,
                                const std::size_t centralAtom,
                                const std::size_t neighborAtom) const;
    /*******************************************************************************************
     * return x,y,z component of pairForces for certain ion-ion pair
     *
     * @param centralAtom index of central Atom in neighbor list ( have same format )
     * @param neighborAtom index of neighbor atom for central atom centralAtom
     * @note for performance call with auto& [ x, y, z ] = get_pairForce( centralAtom, neighborAtom
     *)
     *******************************************************************************************/
    const std::tuple<const Real&, const Real&, const Real&> get_pairForces(
        const String&     key,
        const std::size_t centralAtom,
        const std::size_t neighborAtom) const;
    /*******************************************************************************************
     * return x component of centralForces
     *
     * @param centralAtom index of central Atom in neighbor list ( have same order )
     *******************************************************************************************/
    const Real& get_centralForcesX(const String& key, const std::size_t centralAtom) const;
    /*******************************************************************************************
     * return y component of centralForces
     *
     * @param centralAtom index of central Atom in neighbor list ( have same order )
     *******************************************************************************************/
    const Real& get_centralForcesY(const String& key, const std::size_t centralAtom) const;
    /*******************************************************************************************
     * return z component of centralForces
     *
     * @param centralAtom index of central Atom in neighbor list ( have same order )
     *******************************************************************************************/
    const Real& get_centralForcesZ(const String& key, const std::size_t centralAtom) const;
    /*******************************************************************************************
     * return x,y,z component of centralForces
     *
     * @param centralAtom index of central Atom in neighbor list ( have same order )
     * @note for performance call with
     * auto& [ x, y, z ] = get_centralForces( const std::size_t centralAtom, const std::size_t
     *neighborAtom )
     *******************************************************************************************/
    const std::tuple<const Real&, const Real&, const Real&> get_centralForces(
        const String&     key,
        const std::size_t centralAtom) const;
    /*******************************************************************************************
     * Return x component of atomicForces ( total force acting on ion ).
     *
     * @param centralAtom index of central Atom in neighbor list ( have same order )
     *******************************************************************************************/
    const Real& get_atomicForcesX(const std::size_t centralAtom) const;
    /*******************************************************************************************
     * Return y component of atomicForces ( total force acting on ion ).
     *
     * @param centralAtom index of central Atom in neighbor list ( have same order )
     *******************************************************************************************/
    const Real& get_atomicForcesY(const std::size_t centralAtom) const;
    /*******************************************************************************************
     * Return z component of atomicForces ( total force acting on ion ).
     *
     * @param centralAtom index of central Atom in neighbor list ( have same order )
     *******************************************************************************************/
    const Real& get_atomicForcesZ(const std::size_t centralAtom) const;
    /*******************************************************************************************
     * Return x,y,z component of atomicForces ( total force acting on ion )
     *
     * @param centralAtom index of central Atom in neighbor list ( have same order )
     * @note for performance call with
     * auto& [ x, y, z ] = get_atomicForces( const std::size_t centralAtom, const std::size_t
     *neighborAtom )
     *******************************************************************************************/
    const std::tuple<const Real&, const Real&, const Real&> get_atomicForces(
        const std::size_t centralAtom) const;
    /*******************************************************************************************
     * Computes stress tensor via the use of the Virial-theorem from the pairForces aray.
     *
     * Computes the stressTensor map which stores the stress tensor per descriptor
     * and the totalStressTensor which is the sum over the descriptors.
     *
     * computations are done either in
     * compute_stressTensorCPU( const DescriptorCollector& descriptorCollection, const Real& volume
     *); or compute_stressTensorGPU( const DescriptorCollector& descriptorCollection, const Real&
     *volume ); depending on the chosen execution policy
     *
     * @param descriptorCollection class from which the distances between the atoms can be retrieved
     * @pram volume supply the volume of the lattice in units which agree with the units of the pair
     * forces
     *
     * @f[ \tau_{\alpha,\beta}=\frac{1}{V}\sum_{i=1}^{N_{atoms}}\sum_{j\in
     *\mathcal{N}_{i}}r^{ij}_{\alpha}F^{ij}_{\beta} @f]
     *******************************************************************************************/
    void compute_stressTensor(const DescriptorCollector& descriptorCollection, const Real& volume);
    /*******************************************************************************************
     * Getter for stressTensor component per descriptor.
     *
     * @param key keyword defining for which descriptor the partial stress tensor should be obtained
     * @param indx0 first component of stress tensor for which value is retrieved
     * @param indx1 first component of stress tensor for which value is retrieved
     *
     * @note the storage order is xx, xy, xz, yx, yy, yz, zx, zy, zz
     *                            00  01  02  10  11  12  20  21  22
     *******************************************************************************************/
    const Real& get_stressTensor(const String&     key,
                                 const std::size_t indx0,
                                 const std::size_t indx1) const;
    /*******************************************************************************************
     * Getter for total stress tensor component (totalStressTensor).
     *
     * @param indx0 first component of stress tensor for which value is retrieved
     * @param indx1 first component of stress tensor for which value is retrieved
     *
     * @note the storage order is xx, xy, xz, yx, yy, yz, zx, zy, zz
     *                            00  01  02  10  11  12  20  21  22
     * the total stress tensor was summed over the descriptors
     *******************************************************************************************/
    const Real& get_totalStressTensor(const std::size_t indx0, const std::size_t indx1) const;
    /*******************************************************************************************
     * get total energy of the of considered system
     *******************************************************************************************/
    const Real& get_totalEnergy(void) const;
    /*******************************************************************************************
     * writing atomic forces to screen
     *******************************************************************************************/
    void write_atomicForceToScreen(void) const;

  private:
    /*******************************************************************************************
     * Compute total forces acting on ions.
     *
     * Routine sums up the centralForces and the pairForces terms over the descriptors
     * to obtain the total force acting onto an ion. Routine will be executed on CPU
     *
     *
     * @param descriptorCollection storing normalized descriptors. From there neighbor lists
     * can be retrieved to compute the total force acting on ion
     *******************************************************************************************/
    void compute_atomicForcesCPU(const DescriptorCollector& descriptorCollection);
    /*******************************************************************************************
     * Compute total forces acting on ions.
     *
     * Routine sums up the centralForces and the pairForces terms over the descriptors
     * to obtain the total force acting onto an ion. Routine will be executed on GPU
     *
     *
     * @param descriptorCollection storing normalized descriptors. From there neighbor lists
     * can be retrieved to compute the total force acting on ion
     *******************************************************************************************/
    void compute_atomicForcesGPU(const DescriptorCollector& descriptorCollection);
    /*******************************************************************************************
     * Computes stress tensor via the use of the Virial-theorem from the pairForces aray.
     *
     * Computes the stressTensor map which stores the stress tensor per descriptor
     * and the totalStressTensor which is the sum over the descriptors.
     *
     * This algorithm will be executed on a CPU
     *
     * @param descriptorCollection class from which the distances between the atoms can be retrieved
     * @pram volume supply the volume of the lattice in units which agree with the units of the pair
     * forces
     *
     * @f[ \tau_{\alpha,\beta}=\frac{1}{V}\sum_{i=1}^{N_{atoms}}\sum_{j\in
     *\mathcal{N}_{i}}r^{ij}_{\alpha}F^{ij}_{\beta} @f]
     *******************************************************************************************/
    void compute_stressTensorCPU(const DescriptorCollector& descriptorCollection,
                                 const Real&                volume);
    /*******************************************************************************************
     * Computes stress tensor via the use of the Virial-theorem from the pairForces aray.
     *
     * Computes the stressTensor map which stores the stress tensor per descriptor
     * and the totalStressTensor which is the sum over the descriptors.
     *
     * This algorithm will be executed on a GPU
     *
     * @param descriptorCollection class from which the distances between the atoms can be retrieved
     * @pram volume supply the volume of the lattice in units which agree with the units of the pair
     * forces
     *
     * @f[ \tau_{\alpha,\beta}=\frac{1}{V}\sum_{i=1}^{N_{atoms}}\sum_{j\in
     *\mathcal{N}_{i}}r^{ij}_{\alpha}F^{ij}_{\beta} @f]
     *******************************************************************************************/
    void compute_stressTensorGPU(const DescriptorCollector& descriptorCollection,
                                 const Real&                volume);
    /*******************************************************************************************
     * Sum over the atom in the structure per type and normalize by 1/number ions of the kernel
     matrix.
     *
     *@param kernel used to retrieve kernel matrix @f$K_{iB}@f$
     *@param typeStruc analysed atom type in structure. Atom type of atom number i
     *@param typeForceField analysed atom type in force field. Atom type of atom number i in FF
     *@param atomsPerType number of atom per type.
     *
     *
     *@f[ K( \hat{\mathbf{X}}_{i},\hat{\mathbf{X}}_{B} )=K_{iB}\in
          \mathbb{R}^{N_{atoms\ per\ type}\times N_{reference\ config}}@f]
     *@f[\mathbf{E}_{type}=\sum_{i=1}^{N_{atoms\ per\ type}}\frac{K_{iB}}{N_{atoms}}@f]
     *energyArray has storage form number of types x number reference configurations
     *@f[\mathbf{E} \in \mathbb{R}^{N_{types}\times N_{reference\ config}}@f]
     *******************************************************************************************/
    void compute_energyArray(const Kernel&     kernel,
                             const std::size_t typeStruc,
                             const std::size_t typeForceField,
                             const std::size_t atomsPerType,
                             const Real        scaleFactor);
    /********************************************************************************************
     * Compute local reference configurations times regression coefficients.
     *
     * Product is stored in array descriptorsRefConfsWeighted.
     *
     * @param inputParameters stores force field parameters and is used to retrieve
     * fitting weights @f$\mathbf{w}_{B}@f$ and reference configurations
     * @f$\hat{\mathbf{X}}_{B}@f$
     *
     * computes the Hadamard product over the index iB of
     * @f[\mathbf{X}_{B} \in \mathbb{R}^{number-types\times feature-space-size}@f]
     * @f[\mathbf{w}_{B} \in \mathbb{R}^{number-types\times local-reference-configurations}@f]
     * @f[
         \mathbf{X}_{B}^{'}=\mathbf{X}_{B} \circ \mathbf{w}_{B}
     @f]
     *******************************************************************************************/
    void compute_descriptorsRefConfsWeighted(const IoHandlerML_FF& inputParameters);
    /*******************************************************************************************
     * Compute the total reference energy of considered system.
     *
     * Result is stored to referenceEnergyTotal
     *
     * @f[E_{ref}=\sum_{n=0}^{number\ types}N_{n}*E_{n}^{ref}@f]
     *******************************************************************************************/
    void compute_referenceEnergyTotal(const ShVec1Int& atomsPerType);
    /*******************************************************************************************
     * allocating the energy array
     *
     * first dimension of energy array will be the number of types
     * the second dimension is the number of local reference configurations
     *******************************************************************************************/
    void allocate_energyArray(const std::size_t numberTypesStruc);
    void allocate_energyArrayCPU(const std::size_t numberTypesStruc);
    void allocate_energyArrayGPU(const std::size_t numberTypesStruc);
    void compute_energyArrayDim(const Size_t numberTypesStruc);
    /*******************************************************************************************
     * allocating the derivativeMatrix which stores the derivative of the specific kernel
     *
     * the dimension of the derivativeMatrix is the number of types in the first dimension
     * the second dimension is the number of atoms per type times the number of descriptors
     * the number of descriptors changes faster in the second dimension
     *******************************************************************************************/
    void allocate_derivativeMatrix(const ShVec1Int& numberAtomType);
    void allocate_derivativeMatrixCPU(const ShVec1Int& numberAtomType);
    void allocate_derivativeMatrixGPU(const ShVec1Int& numberAtomType);
    void compute_derivativeMatrixDim(const ShVec1Int& numberAtomType);
    /*******************************************************************************************
     * allocating the derivativeMatrixDescriptorDerivativeProduct
     *
     * the first dimension will be the number of atom in the supplied structure
     * the second dimension is the number of descriptors for the specific atom
     *******************************************************************************************/
    void allocate_forcePreContract(const DescriptorCollector& descriptorCollection,
                                   const std::size_t          numberAtoms);
    void allocate_forcePreContractCPU(const std::size_t numberAtoms);
    void allocate_forcePreContractGPU(const std::size_t numberAtoms);
    void compute_forcePreContractDim(const DescriptorCollector& descriptorCollection,
                                     const std::size_t          numberAtoms);
    /*******************************************************************************************
     * allocating the arrays for the pair forces
     *
     * @param descriptorCollection retrieve neighbor list to get number of neighbors
     * per descriptor and central atom
     * size is chosen in the same way as the connectionVector in neighborList
     * first index is neighbor second index is combination of neighbor xyz
     * neighbor 1 x
     * neighbor 1 y
     * neighbor 1 z
     * neighbor 2 x
     * neighbor 2 y
     * neighbor 2 z
     * .
     * .
     * .
     *******************************************************************************************/
    void allocate_pairForces(const DescriptorCollector& descriptorCollection,
                             const Size_t&              numberAtoms);
    void allocate_pairForcesCPU(const Size_t& numberAtoms);
    void allocate_pairForcesGPU(const Size_t& numberAtoms);
    void compute_pairForcesDim(const DescriptorCollector& descriptorCollection);
    /*******************************************************************************************
     * allocate central force arrays
     *
     * allocating in size number central atoms times 3 (for xyz)
     *******************************************************************************************/
    void allocate_centralForces(const Int numberAtoms);
    void allocate_centralForcesCPU(const Int numberAtoms);
    void allocate_centralForcesGPU(const Int numberAtoms);
    void allocate_centralAtomIndex(const Size_t numberAtoms);
    /*******************************************************************************************
     *allocate stressTensor and totalStressTensor
     *
     * allocation will be done as flattened 3x3 array of size 9
     *******************************************************************************************/
    void allocate_stressTensor(void);
    void allocate_stressTensorCPU(void);
    void allocate_stressTensorGPU(void);
    /*******************************************************************************************
     * allocate atomicForces
     *
     * allocation will be done as flattened Natoms x 3 array
     *******************************************************************************************/
    void allocate_atomicForces(void);
    void allocate_atomicForcesCPU(void);
    void allocate_atomicForcesGPU(void);

    void allocate_tempForceVector(const ShVec1Int& nAtomsType);
    void allocate_tempForceVectorCPU(void);
    void allocate_tempForceVectorGPU(void);
    void compute_tempForceVectorDim(const ShVec1Int& nAtomsType);
    /*******************************************************************************************
     * Computing the energy contribution to the total energy of a specific atom type.
     *
     * computes the dot product between the regression coefficients and energyArray
     @f$\mathbf{E}_{type}@f$
     * @f[\sum\limits_{B}w_{B}\left [ \sum\limits_{i} \frac{1}{N_{\mathrm{atom}}} K(\mathbf{X}_{i},
         \mathbf{X}_{B}) \right ]=\sum\limits_{B}w_{B}E_{type}^{B}@f]
     *******************************************************************************************/
    Real computeEnergyPerType(const std::size_t typeStruc, const std::size_t typeForceField);
    /*******************************************************************************************
     * Adding reference energies to the total energies.
     *
     * depending on the tag TotalEnergyType the energy will be scaled to the isolated atoms
     * or the average energy of the system
     *******************************************************************************************/
    void finalizeEnergyComputation(const ShVec1Int& atomsPerType);
    /*******************************************************************************************
     * Computing the derivative of the kernel matrix
     *
     * derivative matrix @f$\mathbf{L}^{i}=\{L_{ln_{1}n_{2}}^{iJ_{1}J_{2}}\}@f$
     * @f[
        \mathbf{L}^{i} = -\sum_{B} w_{B}
               \frac{\partial K(\hat{\mathbf{X}}_{i},\hat{\mathbf{X}}_{i})}
                    {\partial \hat{\mathbf{X}}_{i}}
       @f]

     * @note for details of the computation of the derivative matrix please check
     * the documentation of the used Kernel class, as for example KernelPolynomial.
     * Function which implements functionality will be called computeDerivativeKernel
     *******************************************************************************************/
    void compute_derivativeMatrix(const Kernel&     kernel,
                                  Vec1Real&         tempForceVector,
                                  const std::size_t typeStruc,
                                  const std::size_t typeForceField);
    /*******************************************************************************************
     * Computes the derivative matrix times the derivatives of the descriptors w.r.t SHS2 and SHS3.
     *
     * @f[\mathbf{\tilde{L}}=\sum\limits_{i}\mathbf{L}_{i} \frac{d\mathbf{X}_{i}}{dp_{nn'l}^{iJJ'}}
     * \frac{dp_{nn'l}^{iJJ'}}{dc_{n''lm}^{iJ''}}@f]
     *
     * dimensions of derivativeMatrix @f$\mathbf{\tilde{L}}\in \mathbb{R}^{N_{atoms}\times lnm }@f$
     *
     * @note implementation is descriptor specifix and the implementing code can be found
     * in the descriptor classes DescriptorSHS2 and DescriptorSHS2
     *******************************************************************************************/
    void compute_forcePreContract(const DescriptorCollector& descriptorCollection,
                                  const TypeMap&             typeMap);
    /*******************************************************************************************
     * Computing force arrays pairForces and centralForces
     *
     * centralForces
     * @f[-\frac{\partial K(\mathbf{X}_{i},\mathbf{X}_{B})}{\partial \mathbf{X}_{i}}*
     *    \frac{d\ \mathbf{X}_{i}}{d\ \mathbf{r}_{i}}@f]
     *
     * pairForces
     * @f[-\sum_{k\in \mathcal{N}_{i}}\frac{\partial K(\mathbf{X}_{i},\mathbf{X}_{B})}{\partial
     *\mathbf{X}_{i}}
     * \frac{d\ \mathbf{X}_{i}}{d\ \mathbf{r}_{k}}@f]
     *
     * @note implementation details can be found in DescriptorSHS2. There the final step is done
     * which is the contraction of forcePreContract and the pair descriptors.
     *******************************************************************************************/
    void computeForceArrays(const DescriptorCollector& descriptorCollection);
    /// number of different atom types stored in force field
    Int numberTypesForceField;
    /*******************************************************************************************
     * number of local reference configurations per force field type
     *******************************************************************************************/
    ShVec1Int numberLocalRefConfs;
    /*******************************************************************************************
     * storing the number of descriptors for 2-body, 3-body... terms
     *******************************************************************************************/
    std::map<String, ShVec1Int> numberDescriptors;
    /*******************************************************************************************
     *  feature space size for 2-body, 3-body... descriptors;
     *  feature space is number of local configurations times number of basis functions
     *******************************************************************************************/
    std::map<String, ShVec1Int> featureSpaceSize;
    /*******************************************************************************************
     * storing weights ( not regression coefficients ) for 2-body,3-body... descriptors
     *******************************************************************************************/
    std::map<String, Real> weights;
    /*******************************************************************************************
     * regression coefficients stored in format (types,local-reference-configurations)
     *******************************************************************************************/
    std::shared_ptr<ShmemArray2DVariableLen<Real>> regressionCoefficients;
    /*******************************************************************************************
     * reference energies of different atom types.
     *
     * @note this can be considered as the energies of a certain atom type in vaccuum
     *******************************************************************************************/
    ShCVec1Real referenceEnergyPerType;
    /*******************************************************************************************
     * average total energy which was ecountered during training of the force field
     *******************************************************************************************/
    Real averageTrainEnergy;
    /*******************************************************************************************
     * SmartEnum to determine if total energy is expressed as isolated atom energy or
     * with respect to average train energy
     *******************************************************************************************/
    TotalEnergyType energyType;
    /*******************************************************************************************
     * total reference energy for current system defined by the supplied Kernel class
     *******************************************************************************************/
    Real referenceEnergyTotal;
    /*******************************************************************************************
     * stores total energy of system after the update routine finished successfully
     *******************************************************************************************/
    Real totalEnergy;
    /*******************************************************************************************
     * typeMap which assigns the force field atom types to the types in the structure and vice versa
     *******************************************************************************************/
    std::shared_ptr<const TypeMap> typeMap;
    /*******************************************************************************************
     * Storing the local reference configurations expressed as descriptors times the regression
     *coefficients.
     *
     * @f[\mathbf{X}_{i_{B}} \in \mathbb{R}^{number-types\times feature-space-size}@f]
     * @f[\mathbf{w}_{i_{B}} \in \mathbb{R}^{number-types\times local-reference-configurations}@f]
     * @f[\mathbf{X}^{'}_{B}=\mathbf{X}_{B} \circ \mathbf{w}_{B}@f]
     *
     * the storage order is types, local reference configurations x numberDescriptors.
     * descriptorsRefConfsWeighted is @f$\mathbf{X}^{'}_{B}\in\mathbb{R}^{N_{types}\times
     * N_{local\ ref\ conf\times N_{descriptors}}}=\mathbb{R}^{N_{types}\times
     *N_{feature\ space}}@f$ numberDescriptors changes fastest in the squashed index
     *******************************************************************************************/
    std::map<String, std::shared_ptr<ShmemArray2DVariableLen<Real>>> descriptorsRefConfsWeighted;
    /*******************************************************************************************
     * Array to store the contraction of the kernel matrix over the atoms times 1/number-ions.
     *
     * @f[
        \sum\limits_{B}w_{B}\sum\limits_{i} \frac{1}{N_{\mathrm{atom}}}
     K(\mathbf{X}_{i},\mathbf{X}_{B})
       @f]

     * storage order energyArray: first index types second index is reference configurations
     * @f$\mathbf{E}\in \mathbb{R}^{N_{types}\times N_{local\ ref\ confs}}@f$
     *******************************************************************************************/
    ShVec2Real      energyArray;
    Vec1Size_t      energyArrayDim;
    ArrayResizing2D energyArraySize;
    /*******************************************************************************************
     * Array to store first part of the Kernel derivative w.r.t to descriptor SHS2/SHS3.
     *
     * for the polynomial Kernel this looks like
     *
     * @f[
        \sum_{B} \Lambda_{iB} \hat{\mathbf{X}}_{B}' - \sum_{B}
     w_{B}\Lambda'_{iB}\hat{\mathbf{X}}_{i}
       @f]
     * with
     * @f[\hat{\mathbf{X}}_{B}' = w_{B} \hat{\mathbf{X}}_{B}@f]
     * @f[\Lambda_{iB} = \frac{\xi}{||\mathbf{X}_{i}||} \left(\hat{\mathbf{X}}_{i} \cdot
     \hat{\mathbf{X}}_{B} \right)^{\xi-1}
       @f]
     * @f[\Lambda'_{iB}= \frac{\xi}{||\mathbf{X}_{i}||} \left(\hat{\mathbf{X}}_{i} \cdot
     \hat{\mathbf{X}}_{B} \right)^{\xi}@f]
     *
     * @note order of array is type index central atom and then squashed index of central atom,
     type1, type2, l, n0, n1
     * is computed in compute_derivativeMatrix
     * @note this matrix is used to compute the forcePreContract array
     *******************************************************************************************/
    std::map<String, ShVec2Real>      derivativeMatrix;
    std::map<String, Vec1Size_t>      derivativeMatrixDim;
    std::map<String, ArrayResizing2D> derivativeMatrixSize;
    /*******************************************************************************************
     * matrix storing the derivative of the kernel times the derivatives of the descriptors
     *
     * storage order is central atom as first index; second index neighbor type times
     * angular number times radial number x second times spherical harmonic parameter m.
     * As symbol @f$\mathbf{\tilde{L}}\in \mathbb{R}^{N_{atoms}\times lnm }@f$
     * @note the derivative of the descriptor with respect to the position is not included
     *******************************************************************************************/
    std::map<String, ShVec2Real>      forcePreContract;
    std::map<String, Vec1Size_t>      forcePreContractDim;
    std::map<String, ArrayResizing2D> forcePreContractSize;
    /*******************************************************************************************
     * stores the pairwise forces obtained by
     *
     * array is computed in routine computeForceArrays( const DescriptorCollector&
     descriptorCollection )
     *
     * @f[
        \sum_{j\in
     \mathcal{N}_{i}}\sum_{B}w_{B}\frac{K(\mathbf{X}_{i},\mathbf{K}_{B})}{d\mathbf{r}_{j}}
       @f]
     *
     * forces are stored in order first index is central atom index;
     * second index is
     * 1st neighbor x
     * 1st neighbor y
     * 1st neighbor z
     * 2nd neighbor x
     * 2nd neighbor y
     * 2nd neighbor z
     * ...
     *******************************************************************************************/
    std::map<String, ShVec2Real>      pairForces;
    std::map<String, Vec1Size_t>      pairForcesDim;
    std::map<String, ArrayResizing2D> pairForcesSize;
    /*******************************************************************************************
     * stores forces obatined by central derivative
     *
     * @f[
       -\sum_{B}w_{B}\frac{K(\mathbf{X}_{i},\mathbf{K}_{B})}{d\mathbf{r}_{i}}
       @f]
     *
     * ordering is atom1 x atom1 y atom1 z, atom2 x atom2 y atom2 z...
     *******************************************************************************************/
    std::map<String, ShVec1Real>      centralForces;
    std::map<String, ArrayResizing1D> centralForcesSize;
    /*******************************************************************************************
     * stress tensor obtained by  virial theorem
     *
     * @f[
          \tau_{\alpha,\beta}=\frac{1}{\Omega}\sum_{i=1}\sum_{j\in \mathcal{N}_{i}}\sum_{B}w_{B}
                              \frac{d\ K(\mathbf{X}_{i},\mathbf{K}_{B})}{dr_{j,\alpha}} *
     r^{ij}_{\beta} =
                              \frac{1}{\Omega}\sum_{i=1}\sum_{j\in \mathcal{N}_{i}}
                              F^{ij}_{\alpha}r^{ij}_{\beta}
       @f]
     *
     * stored in order xx, xy, xz, yx, yy, yz, zx, zy, zz
     *******************************************************************************************/
    std::map<String, ShVec1Real> stressTensor;
    /*******************************************************************************************
     * stress tensor obtained by  virial theorem and summed over various descriptors
     *
     * is sum of stressTensor variable over descriptors
     *
     * stored in order xx, xy, xz, yx, yy, yz, zx, zy, zz
     *******************************************************************************************/
    ShVec1Real      totalStressTensor;
    ArrayResizing1D totalStressTensorSize;
    /*******************************************************************************************
     * stores the total force acting on a single ion
     *
     * the order of the array is atom1 x, atom1 y, atom1 z, atom2 x, atom2 y, atom2 z,...
     *******************************************************************************************/
    ShVec1Real      atomicForces;
    ArrayResizing1D atomicForcesSize;
    /*******************************************************************************************
     * check if force arrays pairForces and centralforces were computed.
     *
     * this is a necessary condition if the atomicForces array,
     * the stressTensor or the totalStressTensor are computed
     *******************************************************************************************/
    bool                  forceArraysComputed;
    ExecutionPolicy       algoExecution;
    linalg::LinalgContext linalgContext;
    Vec1Size_t            centralAtomIndex;
    ArrayResizing1D       centralAtomIndexSize;

    Vec2Real        tempForceVector;
    ArrayResizing2D tempForceVectorSize;
    Vec1Size_t      tempForceVectorDim;
};

} //namespace vaspml

#endif
