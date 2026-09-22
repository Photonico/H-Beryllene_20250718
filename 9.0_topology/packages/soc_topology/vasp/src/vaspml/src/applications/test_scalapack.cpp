#include "MlMPI.hpp"
#include "ScaLapack.hpp"
#include "Timer.hpp"
#include "math.hpp"
#include "types.hpp"

#include <iostream>
#include <string>
#include <vector>

#include <chrono>
#include <thread>

#ifdef VASPML_SCALAPACK
using namespace vaspml;
using namespace vaspml::ScaLapack;

void unit_test(MlMPI& mpiMain)
{

    ScaLapackArray<Real> array1x1(7, 7, &mpiMain, 1, 1);
    //   ScaLapackArray<Real> array2x3( 7, 7, &mpiMain, 2, 3 );
}

void printArray(const ScaLapackArray<Real>& array, const Int rank)
{
    if (array.get_rank() == rank)
    {
        for (Int indx0 = 0; indx0 < array.get_nRowsLoc(); indx0++)
        {
            for (Int indx1 = 0; indx1 < array.get_nColsLoc(); indx1++)
            {
                std::cout << array.getValue(indx0, indx1) << "   ";
            }
            std::cout << std::endl;
        }
    }
}

void timerFunction(const Int size0, const Int size1, MlMPI& mpiMain)
{

    Timer                timer;
    Int                  timerSteps = 500;
    Int                  nRowTot = size0;
    Int                  nColTot = size1;
    ScaLapackArray<Real> stripes;
    ScaLapackArray<Real> checkBoard;
    stripes.init(nRowTot, nColTot, &mpiMain, mpiMain.get_numberRanks(), 1);
    checkBoard.init(nRowTot, nColTot, &mpiMain);
    // filling stripey
    stripes.allocateArray();
    Real value = 0;
    for (Int indx0 = 0; indx0 < stripes.get_nRowsLoc(); indx0++)
    {
        for (Int indx1 = 0; indx1 < stripes.get_nColsLoc(); indx1++)
        {
            stripes.setValue(indx0, indx1, value);
            value++;
        }
    }
    checkBoard.allocateArray();

    for (Int i = 0; i < timerSteps; i++)
    {
        if (i > 10) timer.start("vaspml_pdgemr2d " + std::to_string(nRowTot * nColTot));
        vaspml_pdgemr2d(stripes, checkBoard);
        if (i > 10) timer.stop("vaspml_pdgemr2d " + std::to_string(nRowTot * nColTot));
    }
    if (checkBoard.get_rank() == 0) timer.writeToScreen();
}

void matrixMultiplication(const Int N1, const Int N2, const Int N3, MlMPI& mpiMain)
{

    ScaLapackArray<Real> matrixA(N1, N2, &mpiMain);
    matrixA.allocateArray();
    ScaLapackArray<Real> matrixB(N2, N3, &mpiMain);
    matrixB.allocateArray();
    ScaLapackArray<Real> matrixC(N1, N3, &mpiMain);
    matrixC.allocateArray();
    // fill matrix A, the block cyclic distribution is not visible by the current
    // initialization
    Real value = 0;
    std::cout << "A mtrix " << matrixA.get_nRowsLoc() << "   " << matrixB.get_nColsLoc()
              << std::endl;
    for (Int i = 0; i < matrixA.get_nRowsLoc(); i++)
    {
        for (Int j = 0; j < matrixA.get_nColsLoc(); j++)
        {
            matrixA.setValue(i, j, value);
            value++;
        }
    }

    //char all = 'A';
    //for ( Int i = 0; i < matrixA.get_nprocs(); i++ ){
    //   matrixA.writeToScreen( i, "<<<<<<<<<<< Matrix A "+ std::to_string( i ) + ">>>>>>>>>>>" );
    //   std::this_thread::sleep_for(std::chrono::milliseconds( 1000 ));
    //}

    Int size;
    if (matrixB.get_nColsLoc() < matrixB.get_nRowsLoc()) size = matrixB.get_nColsLoc();
    else size = matrixB.get_nRowsLoc();
    for (Int j = 0; j < size; j++)
    {
        matrixB.setValue(j, j, 1.0);
        value++;
    }
    //for ( Int i = 0; i < matrixB.get_nprocs(); i++ ){
    //   matrixB.writeToScreen( i, "<<<<<<<<<<< Matrix B "+ std::to_string( i ) + ">>>>>>>>>>>" );
    //   std::this_thread::sleep_for(std::chrono::milliseconds( 1000 ));
    //}

    vaspml_pdgemm(transpose::NoTrans, transpose::NoTrans, matrixA, matrixB, matrixC, 1.0, 0.0);

    for (Int i = 0; i < matrixC.get_nprocs(); i++)
    {
        matrixC.writeToScreen(i, "<<<<<<<<<<< Matrix C " + std::to_string(i) + ">>>>>>>>>>>");
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    }
}
#endif

int main(void)
{
#ifdef VASPML_SCALAPACK
    MlMPI mpiMain;
    mpiMain.make_WorldComm();

    //   unit_test( mpiMain );

    //
    //   std::vector<Int> sizes = { 64, 128, 256, 512, 1024, 2048, 4096, 8192 };
    //        //*64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144,
    //        524288 };
    //   for ( auto size : sizes ){
    //      timerFunction( size, size, mpiMain );
    //   }

    matrixMultiplication(5, 8, 4, mpiMain);

    MPI_Finalize();
#else
    std::cout << "I can't do anything. Please add scalapack" << std::endl;
#endif
}
