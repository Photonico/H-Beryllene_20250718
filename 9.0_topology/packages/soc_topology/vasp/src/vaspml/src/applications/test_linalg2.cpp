#include "Linalg.hpp"
#include "ShmemArray.hpp"
#include "types.hpp"

#include <algorithm>
#include <cstddef>
#include <iostream>

using namespace vaspml;

#ifndef NVHPC
class Matrix
{

  public:
    Matrix(std::size_t n, std::size_t m)
    {
        this->n = n;
        this->m = m;
        this->len = n * m;
        data.resize(this->len);
    }

    void setValue(std::size_t indx0, std::size_t indx1, Real value)
    {
        std::size_t indx = m * indx0 + indx1;
        data[indx] = value;
    }

    void printMatrix(void)
    {
        std::size_t col = 0;
        for (std::size_t i = 0; i < n; i++)
        {
            for (std::size_t j = 0; j < m; j++)
            {
                std::cout << data[col] << "   ";
                col++;
            }
            std::cout << std::endl;
        }
    }

    Real* getPointer(void)
    {
        Real& first_element = data[0];
        Real* ptr = &first_element;
        return ptr;
    }

  private:
    std::size_t n;
    std::size_t m;
    std::size_t len;
    Vec1Real    data;
};

int main()
{

#ifndef VASPML_USE_CUBLAS
    Matrix A(3, 2);
    Matrix B(2, 4);
    Matrix C(3, 4);
    A.setValue(0, 0, 1);
    A.setValue(0, 1, 2);
    A.setValue(1, 0, 3);
    A.setValue(1, 1, 4);
    A.setValue(2, 0, 5);
    A.setValue(2, 1, 6);
    A.printMatrix();

    B.setValue(0, 0, 2);
    B.setValue(0, 1, 4);
    B.setValue(0, 2, 6);
    B.setValue(0, 3, 8);
    B.setValue(1, 0, 3);
    B.setValue(1, 1, 5);
    B.setValue(1, 2, 7);
    B.setValue(1, 3, 9);
    B.printMatrix();

    cblas_dgemm(CblasRowMajor,
                CblasNoTrans,
                CblasNoTrans,
                3,
                4,
                2,
                1.0,
                A.getPointer(),
                2,
                B.getPointer(),
                4,
                1.0,
                C.getPointer(),
                4);

    std::cout << "~~~~~~~~~~~~C~~~~~~~~~~~~" << std::endl;
    C.printMatrix();
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

    // one matrix transpose exercise
    Matrix AT(2, 3);
    AT.setValue(0, 0, 1);
    AT.setValue(0, 1, 3);
    AT.setValue(0, 2, 5);
    AT.setValue(1, 0, 2);
    AT.setValue(1, 1, 4);
    AT.setValue(1, 2, 6);
    AT.printMatrix();
    Matrix CT(3, 4);

    cblas_dgemm(CblasRowMajor,
                CblasTrans,
                CblasNoTrans,
                3,
                4,
                2,
                1.0,
                AT.getPointer(),
                3,
                B.getPointer(),
                4,
                1.0,
                CT.getPointer(),
                4);
    //   cublasTgemm(CblasTrans, CblasNoTrans, m, n, k, 1, A, k, B, k, 0, C, m);
    // m, n, k                       lda                ldb                      ldc
    //m Specifies the number of rows of the matrix op(A) and of the matrix C
    //n Specifies the number of columns of the matrix op(B) and the number of columns of the matrix
    //C. k Specifies the number of columns of the matrix op(A) and the number of rows of the matrix
    //op(B) lda Layout = CblasRowMajor, transa=CblasTrans lda must be at least max(1, m). ldb Layout
    //= CblasRowMajor, transb=CblasNoTrans ldb must be at least max(1, n).
    //
    //ldc Layout = CblasRowMajor, ldc must be at least max(1, n).

    B.printMatrix();
    std::cout << "~~~~~~~~~~~~CT~~~~~~~~~~~~" << std::endl;
    CT.printMatrix();
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

    Matrix BT(4, 2);
    BT.setValue(0, 0, 2);
    BT.setValue(1, 0, 4);
    BT.setValue(2, 0, 6);
    BT.setValue(3, 0, 8);
    BT.setValue(0, 1, 3);
    BT.setValue(1, 1, 5);
    BT.setValue(2, 1, 7);
    BT.setValue(3, 1, 9);
    BT.printMatrix();
    Matrix CTT(3, 4);

    cblas_dgemm(CblasRowMajor,
                CblasNoTrans,
                CblasTrans,
                3,
                4,
                2,
                1.0,
                A.getPointer(),
                2,
                BT.getPointer(),
                2,
                1.0,
                CTT.getPointer(),
                4);
    std::cout << "~~~~~~~~~~~~CTT~~~~~~~~~~~~" << std::endl;
    CTT.printMatrix();
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
#endif
}
#else
int main()
{
    std::cout << "Not implemented for NVIDIA compiler." << std::endl;
}
#endif
