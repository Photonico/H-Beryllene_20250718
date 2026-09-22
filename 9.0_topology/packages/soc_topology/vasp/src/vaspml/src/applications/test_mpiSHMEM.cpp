#include "Linalg.hpp"
#include "Matrix.hpp"
#include "MlMPI.hpp"
#include "ShmemArray.hpp"
#include "types.hpp"
#include "utils.hpp"

#include <algorithm>
#include <cstddef>
#include <ctime>
#include <iostream>
#include <memory>
#include <random>

using namespace vaspml;

ShVec1Real make_MatrixDataA(void)
{

    ShVec1Real A = std::make_shared<Vec1Real>(15);
    (*A)[0] = 1;
    (*A)[1] = 2;
    (*A)[2] = 3;
    (*A)[3] = 4;
    (*A)[4] = 5;
    (*A)[5] = 6;
    (*A)[6] = 7;
    (*A)[7] = 8;
    (*A)[8] = 9;
    (*A)[9] = 10;
    (*A)[10] = 11;
    (*A)[11] = 12;
    (*A)[12] = 13;
    (*A)[13] = 14;
    (*A)[14] = 15;
    return A;
}

ShVec1Real make_MatrixDataB(void)
{

    ShVec1Real B = std::make_shared<Vec1Real>(12);
    (*B)[0] = 1;
    (*B)[1] = 4;
    (*B)[2] = 7;
    (*B)[3] = 10;
    (*B)[4] = 2;
    (*B)[5] = 5;
    (*B)[6] = 8;
    (*B)[7] = 11;
    (*B)[8] = 3;
    (*B)[9] = 6;
    (*B)[10] = 9;
    (*B)[11] = 12;
    return B;
}

template<class T>
class A
{

  public:
    A(const std::shared_ptr<MlMPI>& mpiIn, Int in)
    {
        testMPI = mpiIn->make_splitShared();
        size = in;
        tvalue = new T[size];
    }
    ~A(void)
    {
        std::cout << "Descrtuct A" << std::endl;
        delete[] tvalue;
    }

  private:
    MlMPI testMPI;
    Int   size;
    T*    tvalue;
};

int main()
{

    linalg::LinalgContext  linalgContext;
    std::shared_ptr<MlMPI> mpi = std::make_shared<MlMPI>();
    int                    a = 1;
    int                    b = 0;
    mpi->make_WorldComm();

    mpi->allReduceSum(a, b);

    std::shared_ptr<MlMPI> mpiShared = std::make_shared<MlMPI>(mpi->make_splitShared());

    //std::cout << a << " a = b " << b << std::endl;
    Real* dataShmem = nullptr;
    Int   n = 4;
    mpiShared->allocateShmemSpace(dataShmem, n);

    //if ( mpiShared -> get_rank() == 0 ){
    //   double sum = 0;
    //   for ( std::size_t i = 0; i < 1000000000; i++ ){
    //      sum+= ( double ) i;
    //   }
    //   std::cout << "sum = " << sum << std::endl;
    //}
    //array->barrier();

    ShmemArray<Real> array(n, mpi);
    //if ( array[ "mpiInter" ].get_rank() == 0 && array[ "mpiShmem" ].get_rank() == 0 ){
    for (std::size_t i = 0; i < (std::size_t)n; i++) { array.set_value(i, (Real)i, 0); }
    //}
    array.distributeInternode();
    array.print();

    array["mpiInter"].barrier();
    mpi->barrier();
    auto [start, end] = mpi->blockLoop(10);
    std::cout << "Start = " << start << "  end = " << end << std::endl;

    ShmemArray2D<Real> A(5, 3, mpi);

    //if ( array[ "mpiInter" ].get_rank() == 0 && array[ "mpiShmem" ].get_rank() == 0 ){
    A.set_value(0, 0, 1);
    A.set_value(0, 1, 2);
    A.set_value(0, 2, 3);
    A.set_value(1, 0, 4);
    A.set_value(1, 1, 5);
    A.set_value(1, 2, 6);
    A.set_value(2, 0, 7);
    A.set_value(2, 1, 8);
    A.set_value(2, 2, 9);
    A.set_value(3, 0, 10);
    A.set_value(3, 1, 11);
    A.set_value(3, 2, 12);
    A.set_value(4, 0, 13);
    A.set_value(4, 1, 14);
    A.set_value(4, 2, 15);
    //}
    A.distributeInternode();
    A["mpiInter"].barrier();
    ShmemArray2D<Real> B(3, 4, mpi);
    //if ( array[ "mpiInter" ].get_rank() == 0 && array[ "mpiShmem" ].get_rank() == 0 ){
    B.set_value(0, 0, 1);
    B.set_value(0, 1, 4);
    B.set_value(0, 2, 7);
    B.set_value(0, 3, 10);
    B.set_value(1, 0, 2);
    B.set_value(1, 1, 5);
    B.set_value(1, 2, 8);
    B.set_value(1, 3, 11);
    B.set_value(2, 0, 3);
    B.set_value(2, 1, 6);
    B.set_value(2, 2, 9);
    B.set_value(2, 3, 12);
    //}
    B.distributeInternode();
    B["mpiInter"].barrier();

    //std::cout << "OKAY " << std::endl;

    ShmemArray2D<Real> C(5, 4, mpi);

    Vec1Real aa(A.get_dataPointer(), A.get_dataPointer() + 15);
    Vec1Real bb(B.get_dataPointer(), B.get_dataPointer() + 12);
    Vec1Real cc(C.get_dataPointer(), C.get_dataPointer() + 20);

    linalg::matMul(5, 3, 4, aa, bb, cc, linalgContext);

    if (array["mpiInter"].get_rank() == 0 && array["mpiShmem"].get_rank() == 5)
    {
        std::size_t col = 0;
        for (std::size_t i = 0; i < 5; i++)
        {
            for (std::size_t j = 0; j < 4; j++)
            {
                std::cout << cc[col] << "   ";
                col++;
            }
            std::cout << std::endl;
        }
    }

    Vec1Int                       lengths = {5, 3, 6};
    ShmemArray2DVariableLen<Real> varLen(lengths, mpi);
    for (std::size_t i = 0; i < (std::size_t)3; i++)
    {
        for (std::size_t j = 0; j < (std::size_t)lengths[i]; j++)
        {
            varLen.set_value(i, j, (Real)i * j);
        }
    }

    varLen["mpiInter"].barrier();
    for (std::size_t i = 0; i < (std::size_t)3; i++)
    {
        for (std::size_t j = 0; j < (std::size_t)lengths[i]; j++)
        {
            std::cout << varLen.get_value(i, j) << "   ";
        }
        std::cout << std::endl;
    }

    Real* testSlice = varLen.get_slice(2);
    for (std::size_t i = 0; (std::size_t)i < (std::size_t)lengths[2]; i++)
    {
        testSlice[i] = (Real)10;
    }

    varLen["mpiInter"].barrier();
    varLen["mpiShmem"].barrier();

    for (std::size_t i = 0; i < (std::size_t)3; i++)
    {
        std::cout << "Index i " << i << " altered ";
        for (std::size_t j = 0; j < (std::size_t)lengths[i]; j++)
        {
            std::cout << varLen.get_value(i, j) << "   ";
        }
        std::cout << std::endl;
    }

    Real summer = 0;
    mpi->barrier();
    for (std::size_t i = 0; i < 1000000; i++) { summer++; }
    mpi->barrier();
    std::cout << mpi->get_rank() << "  " << summer << std::endl;

    std::size_t      nnn = 3859;
    ShmemArray<Real> random_shmem(nnn, mpi);
    if (mpi->get_rank() == 0)
    {
        std::uniform_real_distribution<Real> unif(-10.0, 10.0);
        std::default_random_engine           re;
        std::vector<Real>                    random(nnn);
        std::for_each(random.begin(), random.end(), [&](Real& x) { x = unif(re); });
        for (std::size_t i = 0; i < nnn; i++) { random_shmem.set_value(i, random[i]); }
    }
    random_shmem["mpiShmem"].barrier();
    random_shmem.distributeInternode();
    random_shmem["mpiInter"].barrier();
    auto file = file_io::openFileO(
        "RandVec" + str("%d", array["mpiShmem"].get_rank() + (array["mpiInter"].get_rank() * 8))
        + ".dat");
    for (std::size_t ndesc = 0; ndesc < random_shmem.get_size(); ndesc++)
    {
        file << str("%24.16E", random_shmem.get_value(ndesc)) << std::endl;
    }
    file.close();

    std::cout << "SHmem array " << random_shmem["mpiShmem"].get_rank() << "   "
              << random_shmem["mpiInter"].get_rank() << "  " << mpi->get_rank() << std::endl;
}
