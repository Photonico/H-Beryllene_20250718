#ifndef KERNELPOLYNOMIAL_HPP
#define KERNELPOLYNOMIAL_HPP

#include "ArrayResizing.hpp"
#include "IoHandlerML_FF.hpp"
#include "Kernel.hpp"
#include "ParallelEnvironemt.hpp"
#include "TypeMap.hpp"
#include "types.hpp"

namespace vaspml
{
/*******************************************************************************************
 * @class KernelPolynomial
 * @brief Inherits from Kernel class and implements PolynomialKernel and needed derivatives
 *
 * This class inherits from Kernel and implements the functionality for a polynomial
 * Kernel. The polynomial kernel needs to implement the derivative terms related purely to
 * the polynomial part of the kernel. These derivative terms are important in force
 * computations. The polynomial kernel will be computed from a Frame class instance
 * which supplies the needed descriptors as DescriptorSHS2, DescriptorSHS3 or DescriptorSHS3Lin.
 * These are used and combined to a DescriptorCollector class from which the kernel
 * with respect to some reference configurations supplied from the IoHandlerML_FF class.
 * The KernelPolynomial can be used in the Predictor class to predict energy, forces and
 * stress tensor
 *******************************************************************************************/
class KernelPolynomial : public Kernel
{

  public:
    /*******************************************************************************************
     * default construction of KernelPolynomial class.
     *
     * @note class is not usable in this state, but can be transformed into
     * an usable state by calling KernelPolynomial::Init(const IoHandlerML_FF& inputParameters,
     * const std::shared_ptr<MlMPI>& mpiIn = nullptr, const ShVec2Real& polynomialDerivative =
     *nullptr, const ShVec2Real& polynomialKernelNorm = nullptr )
     *******************************************************************************************/
    KernelPolynomial(void) = default;
    /*******************************************************************************************
     * @param inputParameters read in from a previously defined force field
     * @param mpiIn when supplied the reference data from the inputParameters can be stored in
     *shared memory
     * @param polynomialDerivative supply when memory for this array should be managed from the
     *outside
     * @param polynomialKernelNorm supply when memory for this array should be managed from the
     *outside
     *******************************************************************************************/
    KernelPolynomial(const IoHandlerML_FF&         inputParameters,
                     const std::shared_ptr<MlMPI>& mpiIn = nullptr,
                     const ShVec2Real&             polynomialDerivative = nullptr,
                     const ShVec2Real&             polynomialKernelNorm = nullptr,
                     const ExecutionPolicy         algoExecution = ExecutionPolicy::cpuSingleCore);
    /*******************************************************************************************
     * initializes the KernelPolynomial when it was first created with the default constructor
     *
     * @param inputParameters read in from a previously defined force field
     * @param mpiIn when supplied the reference data from the inputParameters can be stored in
     *shared memory
     * @param polynomialDerivative supply when memory for this array should be managed from the
     *outside
     * @param polynomialKernelNorm supply when memory for this array should be managed from the
     *outside
     *******************************************************************************************/
    void init(const IoHandlerML_FF&         inputParameters,
              const std::shared_ptr<MlMPI>& mpiIn = nullptr,
              const ShVec2Real&             polynomialDerivative = nullptr,
              const ShVec2Real&             polynomialKernelNorm = nullptr,
              const ExecutionPolicy         algoExecution = ExecutionPolicy::cpuSingleCore);
    /*******************************************************************************************
     * update the polynomial Kernel
     *
     * compute kernelMatrix in the Parent class compute_kernelMatrix
     * allocate or resize the needed working arrays
     * compute the polynomialKernel and terms needed for derivatives
     *
     * @param frame contains the description of the structure to compute kernel of
     *******************************************************************************************/
    void updatePolynomialKernel(const Frame& frame);
    /*******************************************************************************************
     * computing terms polynomialDerivative, polynomialKernelNorm, polynomialKernel
     *******************************************************************************************/
    void computePolynomialKernelTerms(const TypeMap& typeMap);
    /*******************************************************************************************
     * Computing the derivative of the kernel with respect to descriptor vector.
     *
     * @param derivativeMatrix stores the derivative of the kernel
     * @f$K(\hat{\mathbf{X}}_{i},\hat{\mathbf{X}}_{B})=
     * (\hat{\mathbf{X}}_{i}\hat{\mathbf{X}}_{B})^{\xi}@f$ with respect to
     * the descriptor vector @f$\hat{\mathbf{X}}_{i}@f$
     * @param forceVector is a temporary array used to compute the derivativeMatrix which
     * has to be allocated with proper size from call routine. Size is number of atoms per type
     * @f$-\sqrt{w}\Lambda'_{iB}\mathbf{w}_{B}@f$
     * @param descriptorsRefConfsWeighted reference configurations weighted with regression
     coefficients
     * @f$\hat{\mathbf{X}}_{B}^{'}=\mathbf{w}_{B}\circ \hat{\mathbf{X}}_{B}@f$
     * @param typeStruc index defining type of analysed atoms in structure
     * @param typeForceField index defining type of analysed atoms in force field
     *
     * @f[
       \mathbf{L}^{i} = \sum_{B} \Lambda_{iB} \hat{\mathbf{X}}_{B}' - \sum_{B} w_{B}
                        \Lambda'_{iB}\hat{\mathbf{X}}_{i}
       @f]
       with
       @f[\hat{\mathbf{X}}_{B}' = w_{B} \hat{\mathbf{X}}_{B}@f]
       @f[\Lambda_{iB} = \frac{\xi}{||\mathbf{X}_{i}||}
          \left(\hat{\mathbf{X}}_{i} \cdot \hat{\mathbf{X}}_{B} \right)^{\xi-1} @f]
       @f[\Lambda'_{iB}= \frac{\xi}{||\mathbf{X}_{i}||} \left(\hat{\mathbf{X}}_{i}
          \cdot \hat{\mathbf{X}}_{B} \right)^{\xi}@f]
     *******************************************************************************************/
    void computeDerivativeKernel(Vec1Real&          derivativeMatrix,
                                 Vec1Real&          forceVector,
                                 const Real*        descriptorsRefConfsWeighted,
                                 const Real*        regressionCoefficients,
                                 const std::size_t& typeStruc,
                                 const std::size_t& typeForceField,
                                 const String&      key) const;

  private:
    /*******************************************************************************************
     * Actual calculations in computing terms polynomialDerivative, polynomialKernelNorm,
     *polynomialKernel
     *
     * Single cpu version.
     *
     * @param typeMap Map from types of structure to types in force field
     * @param polynomialDerivative derivative of polynomial kernel
     *******************************************************************************************/
    void computePolynomialKernelTermsCPU(
        const TypeMap&                                            typeMap,
        std::vector<std::vector<std::tuple<Real&, Real&, Real&>>> data);

    /*******************************************************************************************
     * Actual calculations in computing terms polynomialDerivative, polynomialKernelNorm,
     *polynomialKernel
     *
     * GPU version, involves standard library algorithms.
     *
     * @param typeMap Map from types of structure to types in force field
     * @param data is a flattened version of { polynomialKernel, polynomialDerivative and
     * polynomialKernelNorm.
     *******************************************************************************************/
    void computePolynomialKernelTermsGPU(
        const TypeMap&                                            typeMap,
        std::vector<std::vector<std::tuple<Real&, Real&, Real&>>> data);

    /*******************************************************************************************
     * allocating the arrays where the terms needed for computing the kernel derivatives in
     * computeDerivativeKernel are allocated.
     *******************************************************************************************/
    void allocateDerivativeArrays(void);
    void allocateDerivativeArraysGPU(void);
    /*******************************************************************************************
     * polynomial power @f$\xi@f$ which is taken of the kernelMatrix
     * @f$K(\hat{\mathbf{X}}_{i},\hat{\mathbf{X}}_{B})=\left (
     *        \hat{\mathbf{X}}_{i}\hat{\mathbf{X}}_{B}\right )^{\xi}@f$
     *******************************************************************************************/
    Int kernelPower;
    /*******************************************************************************************
     * derivative of polynomial part of kernel function
     *
     * @f[ \frac{\xi}{\lVert X_{i}\rVert}K^{\xi-1}_{ij} @f]
     *******************************************************************************************/
    ShVec2Real      polynomialDerivative;
    ArrayResizing2D polynomialDerivativeSize;
    /*******************************************************************************************
     * polynomial kernel times norm factor
     *
     * @f[ \frac{\xi}{\lVert X_{i}\rVert}K^{\xi}_{ij} @f]
     *******************************************************************************************/
    ShVec2Real polynomialKernelNorm;
    /*******************************************************************************************
     * stores the polynomial power of the kernelMatrix
     *
     * points to the memory location of the parent class
     * where kernelMatrix is stored and overwrites it with the polynomial kernel
       @f[K(\hat{\mathbf{X}}_{i},\hat{\mathbf{X}}_{B})=
       \left (\hat{\mathbf{X}}_{i}\hat{\mathbf{X}}_{B}\right )^{\xi}@f]
     *******************************************************************************************/
    ShVec2Real      polynomialKernel;
    Vec1Size_t      centralAtomIndex;
    ArrayResizing1D centralAtomIndexSize;
};

} //namespace vaspml

#endif
