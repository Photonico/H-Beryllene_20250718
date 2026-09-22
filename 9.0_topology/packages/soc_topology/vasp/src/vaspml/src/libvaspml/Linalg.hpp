#ifndef LINALG_HPP
#define LINALG_HPP

#include "ShmemArray.hpp"
#include "Tutor.hpp"
#include "types.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <type_traits>

#ifdef _OPENMP
#include <omp.h>
#endif

// First, decide whether to use CBLAS wrapper for BLAS.
#ifdef VASPML_USE_CBLAS

// Intel MKL uses its own non-standard header for CBLAS.
#if defined(VASPML_USE_MKL)
#include <mkl_cblas.h>
#include <mkl_lapacke.h>
// NVIDIA HPC SDK uses a sub-directory for the LP64 CBLAS header.
#elif defined(__NVCOMPILER)
//#ifdef VASPML_USE_CUBLAS
//// using CUBLAS with the nvidia compiler
//#include <cublas_v2.h>
//#else
#include <lp64/cblas.h>
#include <lp64/lapacke.h>
//#endif
// Others use the standard header, e.g.:
// * OpenBLAS
// * NEC Numeric Library Collection
#else
#include <cblas.h>
// NEC does not provide a C LAPACK interface, we have to use Fortran functions instead.
#ifndef __NEC__
#include <lapacke.h>
#endif
#endif

#elif defined(VASPML_USE_CUBLAS)
#include <cublas_v2.h>
#endif // VASPML_USE_CBLAS

namespace vaspml::linalg
{

#ifdef VASPML_USE_CBLAS

// Intel MKL provides its own integer type.
#if defined(VASPML_USE_MKL)
using linalgInt = MKL_INT;
using transpose = enum CBLAS_TRANSPOSE;
// NEC uses its own integer type.
#elif defined(__NEC__)
using linalgInt = cblas_int_t;
using transpose = CBLAS_TRANSPOSE;
// OpenBLAS commes with its own integer type.
#elif defined(OPENBLAS_VERSION)
using linalgInt = blasint;
using transpose = enum CBLAS_TRANSPOSE;
// Others may use the reference implementation type.
#else
using linalgInt = int;
using transpose = CBLAS_TRANSPOSE;
#endif // VASPML_USE_MKL
       //
inline transpose NoTrans = CblasNoTrans;
inline transpose Trans = CblasTrans;
class LinalgContext
{};
#elif VASPML_USE_CUBLAS
using transpose = cublasOperation_t;
inline transpose NoTrans = CUBLAS_OP_N;
inline transpose Trans = CUBLAS_OP_T;
using linalgInt = Int;

class LinalgContext
{
  public:
    LinalgContext(void) { cublasCreate(&handle); }
    ~LinalgContext(void) { cublasDestroy(handle); }

    const cublasHandle_t& get_handle(void) const { return handle; }

  private:
    cublasHandle_t handle;
};

#else

using linalgInt = Int;

#endif

constexpr linalgInt one = 1;

template<typename T>
const Real* get_data_pointer_const(const T& data)
{
    const Real* ptr;
    if constexpr (std::is_same<Vec1Real, T>::value)
    {
        const Real& first_element = data[0];
        ptr = &first_element;
    }
    else if constexpr (std::is_same<Real*, T>::value) { ptr = data; }
    else if constexpr (std::is_same<const Real*, T>::value) { ptr = data; }
    else if constexpr (std::is_same<ShmemArray<Real>, T>::value) { ptr = data.get_dataPointer(); }
    else if constexpr (std::is_same<ShmemArray2D<Real>, T>::value) { ptr = data.get_dataPointer(); }
    else
    {
        const String type = typeid(T).name();
        //global_scope::tutor.bug( "const Real* ERROR:get_data_pointer_const( const T& data )
        //Supplied type not supported\n"
        //                                            "Supplied type " + type );
        return 0;
    }
    return ptr;
}

template<typename T>
Real* get_data_pointer(T& data)
{
    Real* ptr;
    if constexpr (std::is_same<Vec1Real, T>::value)
    {
        Real& first_element = data[0];
        ptr = &first_element;
    }
    else if constexpr (std::is_same<Real*, T>::value) { ptr = data; }
    else if constexpr (std::is_same<ShmemArray<Real>, T>::value) { ptr = data.get_dataPointer(); }
    else if constexpr (std::is_same<ShmemArray2D<Real>, T>::value) { ptr = data.get_dataPointer(); }
    else
    {
        const String type = typeid(T).name();
        //global_scope::tutor.bug( "Real* get_data_pointer( T& data ) Supplied type not supported\n"
        //                         "Supplied type " + type );
        return 0;
    }
    return ptr;
}

/*******************************************************************************************
 * compute L2 norm of a Real vector
 *
 * @param vector input vector of which norm is computed
 * @param size   number of elements of type T which are stored in the supplied vector
 *
 * @note depending on how the code is compiled the function
 * will use openblas or simple self implemented routines
 *******************************************************************************************/
template<typename T>
Real l2Norm(const T& vector, const Int size, [[maybe_unused]] const LinalgContext& linalgContext)
{
    const Real* ptr = get_data_pointer_const(vector);
#ifdef VASPML_USE_CBLAS
    linalgInt n = size;
    return cblas_dnrm2(n, ptr, one);
#elif VASPML_USE_CUBLAS
    linalgInt n = size;
    Real      result = (Real)0.0;
    cublasDnrm2(linalgContext.get_handle(), n, ptr, one, &result);
    return result;
#else
    Real norm = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : norm)
#endif
    for (std::size_t i = 0; i < (std::size_t)size; i++) { norm += ptr[i] * ptr[i]; }
    return std::sqrt(norm);
#endif
}
/*******************************************************************************************
 * rescale vector by a scalar
 *
 * @param vector  input/output vector which will be rescaled
 * @param factor  Real parameter by which the vector will be rescaled
 *
 * @note depending on how the code is compiled the function
 * will use openblas or simple self implemented routines
 *******************************************************************************************/
template<typename T>
void scaleVector(T&                                    vector,
                 const Real                            factor,
                 const Int                             size,
                 [[maybe_unused]] const LinalgContext& linalgContext)
{
    Real* ptr = get_data_pointer(vector);
#ifdef VASPML_USE_CBLAS
    linalgInt n = size;
    cblas_dscal(n, factor, ptr, one);
#elif VASPML_USE_CUBLAS
    linalgInt n = size;
    cublasDscal(linalgContext.get_handle(), n, &factor, ptr, one);
#else
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (std::size_t i = 0; i < (std::size_t)size; i++) { ptr[i] *= factor; }
#endif
}

/*******************************************************************************************
 * matrix matrix multiplication
 *
 * @param m   leading dimension of A and C
 * @param k   second dimension of A leading dimension of B
 * @param n   second dimension of B and C
 * @param A   @f[ \mathbf{A}\in \mathbb{R}^{ m\times k} @f]
 * @param B   @f[ \mathbf{B}\in \mathbb{R}^{ k\times n} @f]
 * @param C   @f[ \mathbf{C}\in \mathbb{R}^{ m\times n} @f]
 *
 * in detail the function is doing the following computation
 * @f[
 C _{ij} = \sum_{k\_indx=1}^{k}A_{i,k\_indx}B_{k\_indx,j} @f]
 *
 * @note depending on how the code is compiled the function
 * will use opencblas, intel mkl cblas or simple self implemented routines
 *
 * @warning when using the opencblas or intel mkl clbas library the matrix C has to be set to zero
 *  before calling matMul
 *******************************************************************************************/
template<typename T, typename U, typename V>
void matMul(const Int                             m,
            const Int                             k,
            const Int                             n,
            const T&                              A,
            const U&                              B,
            V&                                    C,
            [[maybe_unused]] const LinalgContext& linalgContext)
{
    const Real* A_ptr = get_data_pointer_const(A);
    const Real* B_ptr = get_data_pointer_const(B);
    Real*       C_ptr = get_data_pointer(C);
#ifdef VASPML_USE_CBLAS
    linalgInt      nn = n;
    linalgInt      kk = k;
    linalgInt      mm = m;
    constexpr Real one = (Real)1;
    cblas_dgemm(CblasRowMajor,
                NoTrans,
                NoTrans,
                mm,
                nn,
                kk,
                one,
                A_ptr,
                kk,
                B_ptr,
                nn,
                one,
                C_ptr,
                nn);
#elif VASPML_USE_CUBLAS
    linalgInt      nn = n;
    linalgInt      kk = k;
    linalgInt      mm = m;
    constexpr Real one = (Real)1;
    cublasDgemm(linalgContext.get_handle(),
                NoTrans,
                NoTrans,
                nn,
                mm,
                kk,
                &one,
                B_ptr,
                kk,
                A_ptr,
                nn,
                &one,
                C_ptr,
                nn);
#else
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (std::size_t m_indx = 0; m_indx < (std::size_t)m; m_indx++)
    {
        for (std::size_t n_indx = 0; n_indx < (std::size_t)n; n_indx++)
        {
            Real tmp = 0;
            for (std::size_t k_indx = 0; k_indx < (std::size_t)k; k_indx++)
            {
                tmp += A_ptr[m_indx * k + k_indx] * B_ptr[k_indx * n + n_indx];
            }
            C_ptr[m_indx * n + n_indx] += tmp;
        }
    }
#endif
}

/*******************************************************************************************
 * matrix-matrix multiplication C := alpha*(op(A) + A_offset)*(op(B) + B_offset) + beta*C + C_offset
 *
 * op(X) is either op(X) = X or op(X) = XT,
 * @param transA if transA=CblasNoTrans, then op(A) = A; if transA=CblasTrans, then op(A) = AT.
 * @param transB if transb=CblasNoTrans, then op(B) = B; if transb=CblasTrans, then op(B) = BT.
 * @param m Specifies the number of rows of the matrix op(A) and of the matrix C. The value of m
 *must be at least zero.
 * @param n Specifies the number of columns of the matrix op(B) and the number of columns of the
 * matrix C.
 *          The value of n must be at least zero.
 * @param k Specifies the number of columns of the matrix op(A) and the number of rows of the matrix
 * op(B).
 *          The value of k must be at least zero.
 * @param alpha Specifies the scalar alpha.
 * @param A
 * <table>
 * <caption id="description matrix A">Complex table</caption>
 * <tr><th>transA=CblasNoTrans <th>transA=CblasTrans
 * <tr><td rowspan="2">Array, size lda*m; Before entry, the leading k-by-m part of the array a must
 * contain the matrix A.
 *     <td>Array, size lda*k; Before entry, the leading m-by-k part of the array a must contain the
 * matrix A.
 * </table>
 * @param lda
 * <table>
 * <caption id="description lda">Complex table</caption>
 * <tr><th>transA=CblasNoTrans <th>transA=CblasTrans
 * <tr><td rowspan="2">lda must be at least max(1, k). <td>lda must be at least max(1, m).
 * </table>
 * @param B
 * <table>
 * <caption id="description matrix B">Complex table</caption>
 * <tr><th>transB=CblasNoTrans <th>transB=CblasTrans
 * <tr><td rowspan="2">Array, size ldb by k; Before entry, the leading n-by-k part of the array b
 * must contain the matrix B
 *     <td>Array, size ldb by n; Before entry, the leading k-by-n part of the array b must contain
 * the matrix B
 * </table>
 * @param ldb
 * <table>
 * <caption id="description ldb">Complex table</caption>
 * <tr><th>transA=CblasNoTrans <th>transA=CblasTrans
 * <tr><td rowspan="2">ldb must be at least max(1, n).<td>ldb must be at least max(1, k).
 * </table>
 * @param beta Specifies the scalar beta. When beta is equal to zero, then c need not be set on
 * input.
 * @param C Array, size ldc by m. Before entry, the leading n-by-m part of the array C must contain
 * the matrix C, *except when beta is equal to zero, in which case C need not be set on entry.
 * @param ldc must be at least max(1, n).
 *******************************************************************************************/
template<typename T, typename U, typename V>
void matMul(transpose                             transA,
            transpose                             transB,
            const Int                             m,
            const Int                             n,
            const Int                             k,
            const Real                            alpha,
            const T&                              A,
            Int                                   lda,
            const U&                              B,
            Int                                   ldb,
            Real                                  beta,
            V&                                    C,
            Int                                   ldc,
            [[maybe_unused]] const LinalgContext& linalgContext)
{
#ifdef VASPML_USE_CBLAS
    const Real* A_ptr = get_data_pointer_const(A);
    const Real* B_ptr = get_data_pointer_const(B);
    Real*       C_ptr = get_data_pointer(C);
    linalgInt   nn = n;
    linalgInt   kk = k;
    linalgInt   mm = m;
    linalgInt   lda_blas = lda;
    linalgInt   ldb_blas = ldb;
    linalgInt   ldc_blas = ldc;
    cblas_dgemm(CblasRowMajor,
                transA,
                transB,
                mm,
                nn,
                kk,
                alpha,
                A_ptr,
                lda_blas,
                B_ptr,
                ldb_blas,
                beta,
                C_ptr,
                ldc_blas);
#elif VASPML_USE_CUBLAS
    const Real* A_ptr = get_data_pointer_const(A);
    const Real* B_ptr = get_data_pointer_const(B);
    Real*       C_ptr = get_data_pointer(C);
    linalgInt   nn = n;
    linalgInt   kk = k;
    linalgInt   mm = m;
    linalgInt   lda_blas = lda;
    linalgInt   ldb_blas = ldb;
    linalgInt   ldc_blas = ldc;
    // here the parameters have to be jumbled slightly to get the column major
    // format to which is needed for the cublas calls
    cublasDgemm(linalgContext.get_handle(),
                transB,
                transA,
                nn,
                mm,
                kk,
                &alpha,
                B_ptr,
                ldb_blas,
                A_ptr,
                lda_blas,
                &beta,
                C_ptr,
                ldc_blas);
#else
    global_scope::tutor::warning("You compiled VASPML such that no DGEMM call is available!\n"
                                 "Please select either VASPML_USE_CBLAS or VASPML_USE_CUBLAS\n"
                                 "in your makefile include");
#endif
}

/*******************************************************************************************
 *
 * rescale vector by factor and add this to another vector
 *
 * @param factor rescaling factor of first vector
 * @param scaleVector vector to be rescaled by factor
 * @param addVector vector to which the rescaled vector is added. Contains result
 *
 * @f[addVector[i]= factor * scaleVector[i] + addVector[i] @f]
 *
 *******************************************************************************************/
template<typename T, typename U>
void scaleVectorPlusVector(const Real                            factor,
                           const U&                              scaleVector,
                           T&                                    addVector,
                           const Int                             size,
                           [[maybe_unused]] const LinalgContext& linalgContext)
{
    const Real* scaleVector_ptr = get_data_pointer_const(scaleVector);
    Real*       addVector_ptr = get_data_pointer(addVector);
    linalgInt   n = size;
#ifdef VASPML_USE_CBLAS
    cblas_daxpy(n, factor, scaleVector_ptr, one, addVector_ptr, one);
#elif VASPML_USE_CUBLAS
    cublasDaxpy(linalgContext.get_handle(), n, &factor, scaleVector_ptr, one, addVector_ptr, one);
#else
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (std::size_t i = 0; i < (std::size_t)size; i++)
    {
        addVector_ptr[i] += factor * scaleVector_ptr[i];
    }
#endif
}

/*******************************************************************************************
 * compute dot product of two input vectors
 *
 * @param vectorA first input vector
 * @param vectorB second input vector has to match size of vectorA
 * @param size the number of elements of type T stored in vectorA and VectorB
 *
 * @f[dotProduct=\sum_{i=1}^{N} vectorA[i] * vectorB[i] @f]
 *
 * @note vectorA and vectorB have to be of the same size
 *******************************************************************************************/
template<typename T, typename U>
Real dotProduct(const T&                              vectorA,
                const U&                              vectorB,
                const Int                             size,
                [[maybe_unused]] const LinalgContext& linalgContext)
{

    const Real* vectorA_ptr = get_data_pointer_const(vectorA);
    const Real* vectorB_ptr = get_data_pointer_const(vectorB);
#ifdef VASPML_USE_CBLAS
    linalgInt n = size;
    return cblas_ddot(n, vectorA_ptr, one, vectorB_ptr, one);
#elif VASPML_USE_CUBLAS
    linalgInt n = size;
    Real      norm = 0;
    cublasDdot(linalgContext.get_handle(), n, vectorA_ptr, one, vectorB_ptr, one, &norm);
    return norm;
#else
    Real norm = (Real)0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : norm)
#endif
    for (std::size_t i = 0; i < (std::size_t)size; i++) { norm += vectorA_ptr[i] * vectorB_ptr[i]; }
    return norm;
#endif
}

/*******************************************************************************************
 * matrix-Vector multiplication y := alpha*op(A)*x + beta*y
 *
 * op(X) is either op(X) = X or op(X) = XT,
 * @param trans if trans=CblasNoTrans, then op(A) = A; if transA=CblasTrans, then op(A) = AT.
 * @param m Specifies the number of rows of the matrix op(A) and of the matrix C. The value of m
 *must be at least zero.
 * @param n Specifies the number of columns of the matrix op(A) and the number of rows of the vector
 *x
 * @param alpha Specifies the scalar alpha.
 * @param A
 * <table>
 * <caption id="description matrix A">Complex table</caption>
 * <tr><th>transA=CblasNoTrans <th>transA=CblasTrans
 * <tr><td rowspan="2">Array, size lda*m; Before entry, the leading n-by-m part of the array a must
 * contain the matrix A.
 *     <td>Array, size lda*n; Before entry, the leading m-by-n part of the array a must contain the
 * matrix A.
 * </table>
 * @param lda
 * <table>
 * <caption id="description lda">Complex table</caption>
 * <tr><th>transA=CblasNoTrans <th>transA=CblasTrans
 * <tr><td rowspan="2">lda must be at least max(1, n). <td>lda must be at least max(1, n).
 * </table>
 * @param x vector which is multiplied by A
 * @param incx On entry, incx specifies the increment for the elements of x. incx must not be zero.
 * @param beta Specifies the scalar beta. When beta is equal to zero, then y need not be set on
 *input.
 * @param incy On entry, incy specifies the increment for the elements of y. incx must not be zero.
 *******************************************************************************************/
template<typename T, typename U, typename V>
void matVec(transpose                             trans,
            const Int                             m,
            const Int                             n,
            const Real                            alpha,
            const T&                              A,
            const Int                             lda,
            const U&                              x,
            const Int                             incx,
            const Real                            beta,
            V&                                    y,
            const Int                             incy,
            [[maybe_unused]] const LinalgContext& linalgContext)
{

    const Real* A_ptr = get_data_pointer_const(A);
    const Real* x_ptr = get_data_pointer_const(x);
    Real*       y_ptr = get_data_pointer(y);
#ifdef VASPML_USE_CBLAS
    linalgInt nn = n;
    linalgInt mm = m;
    linalgInt ldaInt = lda;
    linalgInt incxInt = incx;
    linalgInt incyInt = incy;

    cblas_dgemv(CblasRowMajor,
                trans,
                mm,
                nn,
                alpha,
                A_ptr,
                ldaInt,
                x_ptr,
                incxInt,
                beta,
                y_ptr,
                incyInt);
#elif VASPML_USE_CUBLAS
    linalgInt nn = n;
    linalgInt mm = m;
    linalgInt ldaInt = lda;
    linalgInt incxInt = incx;
    linalgInt incyInt = incy;
    // This is fucked up. But since nvidia only allows
    // for column major format which does not exist in this c++ code
    // the transpose statement has to be inverted
    transpose transLoc;
    if (trans == NoTrans) transLoc = Trans;
    else transLoc = NoTrans;
    cublasDgemv(linalgContext.get_handle(),
                transLoc,
                nn,
                mm,
                &alpha,
                A_ptr,
                ldaInt,
                x_ptr,
                incxInt,
                &beta,
                y_ptr,
                incyInt);
#endif
}

/*******************************************************************************************
 * compute inverse of a given square matrix by the use of LU factorization
 *
 * @param matrix matrix that will be inverted. Matrix is overwritten on output
 * @param n0 number of rows in supplied matrix
 * @param n1 number of columns in supplied matrix.
 *
 * @note the supplied matrix has to be in row major format. The only useful format for
 * matrices
 *
 * Step 1: \n
 * The LAPACKE_dgetrf function in the LAPACKE library performs an LU decomposition of a given
 *matrix. LU decomposition factors a matrix AA into two matrices: LL (lower triangular matrix) and
 *UU (upper triangular matrix), such that A=L×UA=L×U. This decomposition is useful for solving
 *linear systems, inverting matrices, and calculating determinants. Step 2: \n The LAPACKE_dgetri
 *function computes the inverse of a matrix using the LU decomposition that was previously obtained
 *from LAPACKE_dgetrf. Here's a step-by-step explanation of what the algorithm does: The algorithm
 *first computes the inverse of the upper triangular matrix U. It then solves for the inverse of the
 *original matrix by considering the lower triangular matrix L. Specifically, it performs forward
 *and backward substitution to find the solution.
 *******************************************************************************************/
#ifdef VASPML_USE_CBLAS
template<typename T>
void computeInverseLU(T& matrix, const Int n0, const Int n1)
{
#ifdef __NEC__
    throw std::runtime_error("ERROR: Matrix inversion not implemented for NEC compiler.\n");
#else
    if (n0 != n1)
    {
        std::cout << "Error " << std::endl;
        //global_scope::tutor.bug( "ERROR: computeInverseLU( T& matrix, const Int n0, const Int n1
        //)\n"
        //                         "Dimension n0 does not match n1 \n" );
    }
    Real*      matrix_ptr = get_data_pointer(matrix);
    linalgInt* ipiv = new int[n0];
    linalgInt  lda = n0;
    linalgInt  info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n0, n1, matrix_ptr, lda, ipiv);
    if (info > 0)
    {
        //global_scope::tutor.error( "ERROR: computeInverseLU( T& matrix, const Int n0, const Int n1
        //)\n"
        //                           "LAPACKE_dgetrf( LAPACK_ROW_MAJOR, n0, n1, matrix_ptr, lda,
        //                           ipiv )\n" "returned error code " + std::to_string( info )
        //                         );
    }
    info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n0, matrix_ptr, lda, ipiv);
    if (info > 0)
    {
        //global_scope::tutor.error( "ERROR: computeInverseLU( T& matrix, const Int n0, const Int n1
        //)\n"
        //                           "LAPACKE_dgetri( LAPACK_ROW_MAJOR, n0, matrix_ptr, lda, ipiv
        //                           )\n" "returned error code " + std::to_string( info )
        //                         );
    }
#endif
}

#elif defined(VASPML_USE_CUBLAS)
template<typename T>
void computeInverseLU(T&, const Int, const Int)
{}
#endif

#ifdef VASPML_USE_CUBLAS
/*******************************************************************************************
 * This routine can be used to copy data from the GPU back to the CPU
 *
 * This function copies n elements from a vector x in GPU memory space to a vector
 * y in host memory space. Elements in both vectors are assumed to have a size of elemSize
 * bytes. The storage spacing between consecutive elements is given by incx for the source
 * vector and incy for the destination vector y.
 *
 * @param[in] data1 data which is currently located on the GPU
 * @param[inout] data2 data to which the GPU data is copied to be on CPU
 * @param[in] n is the number of elements in data1 and data2
 *
 * @note data1 and data2 can be the same variable
 *******************************************************************************************/
template<typename T>
void getVectorGPU(const T& data1, T& data2, const Int n)
{
    const Real* data1Ptr = get_data_pointer_const(data1);
    Real*       data2Ptr = get_data_pointer(data2);
    cublasGetVector(n, sizeof(T), data1Ptr, one, data2Ptr, one);
}

/*******************************************************************************************
 * This routine can be used to copy data from the CPU to the GPU
 *
 * This function copies n elements from a vector x in host memory space to a vector y in GPU
 * memory space. Elements in both vectors are assumed to have a size of elemSize bytes.
 * The storage spacing between consecutive elements is given by incx for the source vector x
 * and by incy for the destination vector y.
 *
 * @param[in] data1 data which is currently located on the CPU
 * @param[inout] data2 data to which the CPU data is copied to be on GPU
 * @param[in] n is the number of elements in data1 and data2
 *
 * @note data1 and data2 can be the same variable
 *******************************************************************************************/
template<typename T>
void setVectorGPU(const T& data1, T& data2, const Int n)
{
    const Real* data1Ptr = get_data_pointer_const(data1);
    Real*       data2Ptr = get_data_pointer(data2);
    cublasGetVector(n, sizeof(T), data1Ptr, one, data2Ptr, one);
}
#endif

} // namespace vaspml::linalg
#endif
