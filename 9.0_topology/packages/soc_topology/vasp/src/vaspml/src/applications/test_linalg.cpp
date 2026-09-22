#include "Linalg.hpp"
#include "ShmemArray.hpp"
#include "types.hpp"

#include <algorithm>
#include <iostream>

using namespace vaspml;

int main()
{

    Vec1Real              test_vec(100, 1);
    ShmemArray<Real>      test_vec2(100);
    linalg::LinalgContext linalgContext;
    for (std::size_t i = 0; i < test_vec2.get_size(); i++) test_vec2.set_value(i, 1);

    Real norm = linalg::l2Norm(test_vec, test_vec.size(), linalgContext);
    std::cout << norm << std::endl;
    norm = linalg::l2Norm(test_vec2, test_vec2.get_size(), linalgContext);
    //norm = linalg::l2Norm( test_vec2, test_vec2.get_size() );
    std::cout << "Shmem array " << norm << std::endl;

    Real factor = 2;
    linalg::scaleVector(test_vec, factor, test_vec.size(), linalgContext);
    linalg::scaleVector(test_vec2.get_dataPointer(), factor, test_vec2.get_size(), linalgContext);

    norm = linalg::l2Norm(test_vec, test_vec.size(), linalgContext);
    std::cout << norm << std::endl;

    norm = linalg::l2Norm(test_vec2.get_dataPointer(), test_vec2.get_size(), linalgContext);
    std::cout << "Shared pointer " << norm << std::endl;

    // matmul test
    Vec1Real A(6);
    A[0] = 1;
    A[1] = 2;
    A[2] = 3;
    A[3] = 4;
    A[4] = 5;
    A[5] = 6;

    Vec1Real B(6);
    B[0] = 1;
    B[1] = 3;
    B[2] = 5;
    B[3] = 2;
    B[4] = 4;
    B[5] = 6;

    Vec1Real C(9);

    linalg::matMul(3, 2, 3, A, B, C, linalgContext);

    std::cout << "Matrix A " << std::endl;
    size_t col = 0;
    for (size_t i = 0; i < 3; i++)
    {
        for (size_t j = 0; j < 2; j++)
        {
            std::cout << A[col] << "   ";
            col++;
        }
        std::cout << std::endl;
    }

    std::cout << "Matrix B " << std::endl;
    col = 0;
    for (size_t i = 0; i < 2; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            std::cout << B[col] << "   ";
            col++;
        }
        std::cout << std::endl;
    }

    std::cout << "Matrix C " << std::endl;
    col = 0;
    for (size_t i = 0; i < 3; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            std::cout << C[col] << "   ";
            col++;
        }
        std::cout << std::endl;
    }

    std::cout << "Shared memoery array matmul test " << std::endl;

    ShmemArray<Real> Ashmem(6);
    Ashmem.set_value(0, 1);
    Ashmem.set_value(1, 2);
    Ashmem.set_value(2, 3);
    Ashmem.set_value(3, 4);
    Ashmem.set_value(4, 5);
    Ashmem.set_value(5, 6);

    ShmemArray<Real> Bshmem(6);
    Bshmem.set_value(0, 1);
    Bshmem.set_value(1, 3);
    Bshmem.set_value(2, 5);
    Bshmem.set_value(3, 2);
    Bshmem.set_value(4, 4);
    Bshmem.set_value(5, 6);

    ShmemArray<Real> Cshmem(9);

    linalg::matMul(3,
                   2,
                   3,
                   Ashmem.get_dataPointer(),
                   Bshmem.get_dataPointer(),
                   Cshmem.get_dataPointer(),
                   linalgContext);

    std::cout << "Shmem Matrix A " << std::endl;
    col = 0;
    for (size_t i = 0; i < 3; i++)
    {
        for (size_t j = 0; j < 2; j++)
        {
            std::cout << Ashmem.get_value(col) << "   ";
            col++;
        }
        std::cout << std::endl;
    }

    std::cout << "Shmem Matrix B " << std::endl;
    col = 0;
    for (size_t i = 0; i < 2; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            std::cout << Bshmem.get_value(col) << "   ";
            col++;
        }
        std::cout << std::endl;
    }

    std::cout << "Shmem Matrix C " << std::endl;
    col = 0;
    for (size_t i = 0; i < 3; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            std::cout << Cshmem.get_value(col) << "   ";
            col++;
        }
        std::cout << std::endl;
    }

    std::fill(C.begin(), C.end(), 0);
    linalg::matMul(3, 2, 3, Ashmem.get_dataPointer(), B, C, linalgContext);
    std::cout << "Shmem Matrix C from mixed multiply" << std::endl;
    col = 0;
    for (size_t i = 0; i < 3; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            std::cout << C[col] << "   ";
            col++;
        }
        std::cout << std::endl;
    }

    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "Next test " << std::endl;
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

    A.resize(12);
    B.resize(12);

    A[0] = 1;
    A[1] = 2;
    A[2] = 3;
    A[3] = 4;
    A[4] = 5;
    A[5] = 6;
    A[6] = 7;
    A[7] = 8;
    A[8] = 9;
    A[9] = 10;
    A[10] = 11;
    A[11] = 12;

    B[0] = 1;
    B[1] = 4;
    B[2] = 7;
    B[3] = 10;
    B[4] = 2;
    B[5] = 5;
    B[6] = 8;
    B[7] = 11;
    B[8] = 3;
    B[9] = 6;
    B[10] = 9;
    B[11] = 12;

    C.resize(12);
    std::fill(C.begin(), C.end(), 0);

    linalg::matMul(4, 3, 4, A, B, C, linalgContext);

    std::cout << "Matrix A " << std::endl;
    col = 0;
    for (size_t i = 0; i < 4; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            std::cout << A[col] << "   ";
            col++;
        }
        std::cout << std::endl;
    }

    std::cout << "Matrix B " << std::endl;
    col = 0;
    for (size_t i = 0; i < 3; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            std::cout << B[col] << "   ";
            col++;
        }
        std::cout << std::endl;
    }

    std::cout << "Matrix C " << std::endl;
    col = 0;
    for (size_t i = 0; i < 4; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            std::cout << C[col] << "   ";
            col++;
        }
        std::cout << std::endl;
    }

    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "Next test " << std::endl;
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

    A.resize(15);
    B.resize(12);

    A[0] = 1;
    A[1] = 2;
    A[2] = 3;
    A[3] = 4;
    A[4] = 5;
    A[5] = 6;
    A[6] = 7;
    A[7] = 8;
    A[8] = 9;
    A[9] = 10;
    A[10] = 11;
    A[11] = 12;
    A[12] = 13;
    A[13] = 14;
    A[14] = 15;

    B[0] = 1;
    B[1] = 4;
    B[2] = 7;
    B[3] = 10;
    B[4] = 2;
    B[5] = 5;
    B[6] = 8;
    B[7] = 11;
    B[8] = 3;
    B[9] = 6;
    B[10] = 9;
    B[11] = 12;

    C.resize(20);
    std::fill(C.begin(), C.end(), 0);

    linalg::matMul(5, 3, 4, A, B, C, linalgContext);

    std::cout << "Matrix A " << std::endl;
    col = 0;
    for (size_t i = 0; i < 5; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            std::cout << A[col] << "   ";
            col++;
        }
        std::cout << std::endl;
    }

    std::cout << "Matrix B " << std::endl;
    col = 0;
    for (size_t i = 0; i < 3; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            std::cout << B[col] << "   ";
            col++;
        }
        std::cout << std::endl;
    }

    std::cout << "Matrix C " << std::endl;
    col = 0;
    for (size_t i = 0; i < 5; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            std::cout << C[col] << "   ";
            col++;
        }
        std::cout << std::endl;
    }
}
