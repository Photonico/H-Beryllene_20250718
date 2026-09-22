#include "Linalg.hpp"
#include "Timer.hpp"
#include "types.hpp"

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <random>
#include <string>

#ifdef OPENBLAS
#ifdef USEMKL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif
#endif

#ifdef OPENMP
#include <omp.h>
#endif

using namespace vaspml;

Vec1Real generateRandomMatrix(const size_t n, const size_t m)
{

    Vec1Real matrix(n * m);

    std::random_device               seed;
    std::mt19937                     rng(seed());
    std::uniform_real_distribution<> dist(0.0, 1.0);

    for (auto& x : matrix)
    {
        Real y = dist(rng);
        y = y * (Real)2 - (Real)1;
        x = y;
    }

    return matrix;
}

int main(void)
{

    Timer                 timer;
    std::size_t           timer_steps = 500;
    std::size_t           max_size = 1000;
    linalg::LinalgContext linalgContext;

    //#ifdef OPENBLAS
    //   openblas_set_num_threads( 8 );
    //#endif

#ifdef __INTEL_MKL__
    std::cout << "INTEL MKL is used" << std::endl;
#endif

    for (std::size_t mat_size = 100; mat_size < max_size; mat_size += 100)
    {
        Vec1Real A = generateRandomMatrix(mat_size, mat_size);
        Vec1Real B = generateRandomMatrix(mat_size, mat_size);
        Vec1Real C(mat_size * mat_size);
        for (std::size_t j = 0; j < timer_steps; j++)
        {
            if (j > 10) timer.start("matMul " + std::to_string(mat_size));
            linalg::matMul(mat_size, mat_size, mat_size, A, B, C, linalgContext);
            if (j > 10) timer.stop("matMul " + std::to_string(mat_size));
        }
        std::cout << "mat_size " << mat_size << " done of " << max_size << std::endl;
    }

    timer.writeToScreen();
}
