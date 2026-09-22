#include "Linalg.hpp"
#include "types.hpp"

#include "Timer.hpp"

#include <iostream>

using namespace vaspml;

Vec1Real makeMatrix(const Int n1, const Int n2)
{

    Vec1Real matrix(n1 * n2);
    Real     value = 0;
    for (Int i = 0; i < n1; i++)
    {
        for (Int j = 0; j < n2; j++)
        {
            //         if ( i == j ) matrix[ i*n2 + j ]  = ( Real ) 2.0;
            matrix[i * n2 + j] = value;
            value++;
        }
    }
    return matrix;
}

void writeMatrix(const Vec1Real& data, const Int n1, const Int n2)
{
    for (Int i = 0; i < n1; i++)
    {
        for (Int j = 0; j < n2; j++) { std::cout << data[i * n2 + j] << "   "; }
        std::cout << std::endl;
    }
}

int main()
{

    Timer timer;
    timer.start("l2Norm");

#ifdef VASPML_USE_CUBLAS
    linalg::LinalgContext linalgContext;
    // vector operations
    Vec1Real test_vec(100000, 1);
    Real     norm = 0;
    for (Int i = 0; i < 100; i++) norm += linalg::l2Norm(test_vec, test_vec.size(), linalgContext);
    timer.stop("l2Norm");

    std::cout << norm << std::endl;

    Real alpha = (Real)2;

    timer.start("scaleVector");
    for (Int i = 0; i < 100; i++)
        linalg::scaleVector(test_vec, alpha, test_vec.size(), linalgContext);
    timer.stop("scaleVector");

    for (Int i = 0; i < 5; i++) std::cout << test_vec[i] << "  ";
    std::cout << std::endl;

    Vec1Real vector1(100000, 1);
    Vec1Real vector2(100000, 1);
    alpha = (Real)3.0;

    timer.start("scaleVectorPlusVector");
    //for ( Int i = 0; i < 100; i++ )
    linalg::scaleVectorPlusVector(alpha, vector1, vector2, vector1.size(), linalgContext);
    //linalg::getVectorGPU( vector2, vector2, vector2.size(), linalgContext );
    timer.stop("scaleVectorPlusVector");

    for (Int i = 0; i < 5; i++) std::cout << vector1[i] << "  " << vector2[i] << "  ";
    std::cout << std::endl;

    Real dotP = (Real)0;
    timer.start("dotProduct");
    for (Int i = 0; i < 100; i++)
        dotP += linalg::dotProduct(vector1, vector2, vector1.size(), linalgContext);
    timer.stop("dotProduct");
    std::cout << "Dot product " << dotP << std::endl;

    // matrix times vector multiplication
    Int      n1 = 5;
    Int      n2 = 8;
    Vec1Real matrixA = makeMatrix(n1, n2);
    Vec1Real vecA(n2, 1);
    Vec1Real vecB(n1);
    Int      one = 1;
    alpha = 1.0;

    linalg::matVec(linalg::NoTrans,
                   n1,
                   n2,
                   alpha,
                   matrixA,
                   n2,
                   vecA,
                   one,
                   alpha,
                   vecB,
                   one,
                   linalgContext);
    linalg::getVectorGPU(vector2, vector2, vector2.size());

    for (Int i = 0; i < n1; i++) std::cout << vecB[i] << "  ";
    std::cout << std::endl;

    Int n3;
    n1 = 5; // m
    n2 = 8; // k
    n3 = 4; // n
    alpha = 1.0;
    Real beta = 0.0;
    matrixA = makeMatrix(n1, n2);
    Vec1Real matrixB = makeMatrix(n2, n3);
    std::cout << "<<<<<<<<<<<<<<<<<< matrixA >>>>>>>>>>>>>>>>>>>" << std::endl;
    writeMatrix(matrixA, n1, n2);
    std::cout << "<<<<<<<<<<<<<<<<<< matrixB >>>>>>>>>>>>>>>>>>>" << std::endl;
    writeMatrix(matrixB, n2, n3);
    Vec1Real matrixC(n1 * n3);

    //   linalg::matMul( linalg::NoTrans, linalg::NoTrans,
    //                   n2, n3, n1, alpha,
    //                   matrixA, n1,
    //                   matrixB, n2,
    //                   0.0,
    //                   matrixC,
    //                   n2, linalgContext
    //                 );

    cublasHandle_t handle;
    cublasCreate(&handle);

    //   this properly multiplies the transpose matrices correctly
    //   checked with python3 B.T * A.T
    //   cublasDgemm( handle,
    //                linalg::Trans,
    //                linalg::Trans,
    //                n1,n3,n2,
    //                &alpha,
    //                matrixA.data(),
    //                n2,
    //                matrixB.data(),
    //                n3,
    //                &beta,
    //                matrixC.data(),
    //                n1, linalgContext
    //         );
    // this works. This gives A * B

    //   cublasDgemm( handle,
    //                linalg::NoTrans,
    //                linalg::NoTrans,
    //                n3,n1,n2,
    //                &alpha,
    //                matrixB.data(),
    //                n3,
    //                matrixA.data(),
    //                n2,
    //                &beta,
    //                matrixC.data(),
    //                n3, linalgContext
    //              );
    linalg::matMul(linalg::NoTrans,
                   linalg::NoTrans,
                   n1,
                   n3,
                   n2,
                   alpha,
                   matrixA,
                   n2,
                   matrixB,
                   n3,
                   beta,
                   matrixC,
                   n3,
                   linalgContext);

    linalg::getVectorGPU(matrixC, matrixC, matrixC.size());
    std::cout << "<<<<<<<<<<<<<<<<<< matrixC >>>>>>>>>>>>>>>>>>>" << std::endl;
    writeMatrix(matrixC, n1, n3);
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    // test 2 A B.T
    n1 = 5; // m
    n2 = 8; // k
    n3 = 4; // n
    alpha = 1.0;
    beta = 0.0;
    matrixA = makeMatrix(n1, n2);
    matrixB = makeMatrix(n3, n2);
    std::cout << "<<<<<<<<<<<<<<<<<< matrixA >>>>>>>>>>>>>>>>>>>" << std::endl;
    writeMatrix(matrixA, n1, n2);
    std::cout << "<<<<<<<<<<<<<<<<<< matrixB >>>>>>>>>>>>>>>>>>>" << std::endl;
    writeMatrix(matrixB, n3, n2);
    matrixC.resize(n1 * n3, 0);
    std::fill(matrixC.begin(), matrixC.end(), 0);

    //   cublasDgemm( handle,
    //                linalg::Trans,
    //                linalg::NoTrans,
    //                n3,n1,n2,
    //                &alpha,
    //                matrixB.data(),
    //                n2,
    //                matrixA.data(),
    //                n2,
    //                &beta,
    //                matrixC.data(),
    //                n3
    //              );
    linalg::matMul(linalg::NoTrans,
                   linalg::Trans,
                   n1,
                   n3,
                   n2,
                   alpha,
                   matrixA,
                   n2,
                   matrixB,
                   n2,
                   beta,
                   matrixC,
                   n3,
                   linalgContext);

    linalg::getVectorGPU(matrixC, matrixC, matrixC.size());
    std::cout << "<<<<<<<<<<<<<<<<<< matrixC >>>>>>>>>>>>>>>>>>>" << std::endl;
    writeMatrix(matrixC, n1, n3);

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    // test 3 A.T B
    n1 = 5; // m
    n2 = 8; // k
    n3 = 4; // n
    alpha = 1.0;
    beta = 0.0;
    matrixA = makeMatrix(n2, n1);
    matrixB = makeMatrix(n2, n3);
    std::cout << "<<<<<<<<<<<<<<<<<< matrixA >>>>>>>>>>>>>>>>>>>" << std::endl;
    writeMatrix(matrixA, n2, n1);
    std::cout << "<<<<<<<<<<<<<<<<<< matrixB >>>>>>>>>>>>>>>>>>>" << std::endl;
    writeMatrix(matrixB, n3, n2);
    matrixC.resize(n1 * n3, 0);
    std::fill(matrixC.begin(), matrixC.end(), 0);

    //cublasDgemm( handle,
    //             linalg::NoTrans,
    //             linalg::Trans,
    //             n3,n1,n2,
    //             &alpha,
    //             matrixB.data(),
    //             n3,
    //             matrixA.data(),
    //             n1,
    //             &beta,
    //             matrixC.data(),
    //             n3, linalgContext
    //           );
    linalg::matMul(linalg::Trans,
                   linalg::NoTrans,
                   n1,
                   n3,
                   n2,
                   alpha,
                   matrixA,
                   n1,
                   matrixB,
                   n3,
                   beta,
                   matrixC,
                   n3,
                   linalgContext);

    linalg::getVectorGPU(matrixC, matrixC, matrixC.size());
    std::cout << "<<<<<<<<<<<<<<<<<< matrixC >>>>>>>>>>>>>>>>>>>" << std::endl;
    writeMatrix(matrixC, n1, n3);

    // test $ A.T B.T
    n1 = 5; // m
    n2 = 8; // k
    n3 = 4; // n
    alpha = 1.0;
    beta = 0.0;
    matrixA = makeMatrix(n2, n1);
    matrixB = makeMatrix(n3, n1);
    std::cout << "<<<<<<<<<<<<<<<<<< matrixA >>>>>>>>>>>>>>>>>>>" << std::endl;
    writeMatrix(matrixA, n2, n1);
    std::cout << "<<<<<<<<<<<<<<<<<< matrixB >>>>>>>>>>>>>>>>>>>" << std::endl;
    writeMatrix(matrixB, n3, n2);
    matrixC.resize(n1 * n3, 0);
    std::fill(matrixC.begin(), matrixC.end(), 0);

    //cublasDgemm( handle,
    //             linalg::Trans,
    //             linalg::Trans,
    //             n3,n1,n2,
    //             &alpha,
    //             matrixB.data(),
    //             n2,
    //             matrixA.data(),
    //             n1,
    //             &beta,
    //             matrixC.data(),
    //             n3, linalgContext
    //           );
    linalg::matMul(linalg::Trans,
                   linalg::Trans,
                   n1,
                   n3,
                   n2,
                   alpha,
                   matrixA,
                   n1,
                   matrixB,
                   n2,
                   beta,
                   matrixC,
                   n3,
                   linalgContext);

    linalg::getVectorGPU(matrixC, matrixC, matrixC.size());
    std::cout << "<<<<<<<<<<<<<<<<<< matrixC >>>>>>>>>>>>>>>>>>>" << std::endl;
    writeMatrix(matrixC, n1, n3);
    cublasDestroy(handle);

#elif defined(VASPML_USE_CBLAS)
    linalg::LinalgContext linalgContext;
    // now I have to redo the whole exercise with the interface function

    // test 1 A * B
    Int      n1 = 5; // m
    Int      n2 = 8; // k
    Int      n3 = 4; // n
    Real     alpha = 1.0;
    Real     beta = 0.0;
    Vec1Real matrixA = makeMatrix(n1, n2);
    Vec1Real matrixB = makeMatrix(n2, n3);
    Vec1Real matrixC(n1 * n3);
    std::cout << "<<<<<<<<<<<<<<<<<< matrixA >>>>>>>>>>>>>>>>>>>" << std::endl;
    writeMatrix(matrixA, n1, n2);
    std::cout << "<<<<<<<<<<<<<<<<<< matrixB >>>>>>>>>>>>>>>>>>>" << std::endl;
    writeMatrix(matrixB, n2, n3);
    linalg::matMul(linalg::NoTrans,
                   linalg::NoTrans,
                   n1,
                   n3,
                   n2,
                   alpha,
                   matrixA,
                   n2,
                   matrixB,
                   n3,
                   beta,
                   matrixC,
                   n3,
                   linalgContext);

    std::cout << "<<<<<<<<<<<<<<<<<< matrixC >>>>>>>>>>>>>>>>>>>" << std::endl;
    writeMatrix(matrixC, n1, n3);

    // test 3 A.T * B.T
    n1 = 5; // m
    n2 = 8; // k
    n3 = 4; // n
    alpha = 1.0;
    beta = 0.0;
    matrixA = makeMatrix(n1, n2);
    matrixB = makeMatrix(n3, n2);
    matrixC.resize(n3 * n1);
    std::cout << "<<<<<<<<<<<<<<<<<< matrixA >>>>>>>>>>>>>>>>>>>" << std::endl;
    writeMatrix(matrixA, n1, n2);
    std::cout << "<<<<<<<<<<<<<<<<<< matrixB >>>>>>>>>>>>>>>>>>>" << std::endl;
    writeMatrix(matrixB, n3, n2);
    linalg::matMul(linalg::NoTrans,
                   linalg::Trans,
                   n1,
                   n3,
                   n2,
                   alpha,
                   matrixA,
                   n2,
                   matrixB,
                   n2,
                   beta,
                   matrixC,
                   n3,
                   linalgContext);

    std::cout << "<<<<<<<<<<<<<<<<<< matrixC >>>>>>>>>>>>>>>>>>>" << std::endl;
    writeMatrix(matrixC, n1, n3);

    // test 3 A.T * B
    n1 = 5; // m
    n2 = 8; // k
    n3 = 4; // n
    alpha = 1.0;
    beta = 0.0;
    matrixA = makeMatrix(n2, n1);
    matrixB = makeMatrix(n2, n3);
    matrixC.resize(n3 * n1);
    std::cout << "<<<<<<<<<<<<<<<<<< matrixA >>>>>>>>>>>>>>>>>>>" << std::endl;
    writeMatrix(matrixA, n2, n1);
    std::cout << "<<<<<<<<<<<<<<<<<< matrixB >>>>>>>>>>>>>>>>>>>" << std::endl;
    writeMatrix(matrixB, n2, n3);
    linalg::matMul(linalg::Trans,
                   linalg::NoTrans,
                   n1,
                   n3,
                   n2,
                   alpha,
                   matrixA,
                   n1,
                   matrixB,
                   n3,
                   beta,
                   matrixC,
                   n3,
                   linalgContext);

    std::cout << "<<<<<<<<<<<<<<<<<< matrixC >>>>>>>>>>>>>>>>>>>" << std::endl;
    writeMatrix(matrixC, n1, n3);

    // test 4 A.T * B.T
    n1 = 5; // m
    n2 = 8; // k
    n3 = 4; // n
    alpha = 1.0;
    beta = 0.0;
    matrixA = makeMatrix(n2, n1);
    matrixB = makeMatrix(n3, n2);
    matrixC.resize(n3 * n1);
    std::cout << "<<<<<<<<<<<<<<<<<< matrixA >>>>>>>>>>>>>>>>>>>" << std::endl;
    writeMatrix(matrixA, n2, n1);
    std::cout << "<<<<<<<<<<<<<<<<<< matrixB >>>>>>>>>>>>>>>>>>>" << std::endl;
    writeMatrix(matrixB, n2, n3);
    linalg::matMul(linalg::Trans,
                   linalg::Trans,
                   n1,
                   n3,
                   n2,
                   alpha,
                   matrixA,
                   n1,
                   matrixB,
                   n2,
                   beta,
                   matrixC,
                   n3,
                   linalgContext);

    std::cout << "<<<<<<<<<<<<<<<<<< matrixC >>>>>>>>>>>>>>>>>>>" << std::endl;
    writeMatrix(matrixC, n1, n3);

#endif

    timer.writeToScreen();
}
