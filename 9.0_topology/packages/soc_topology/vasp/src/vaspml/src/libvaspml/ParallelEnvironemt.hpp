#ifndef PARALLELENVIRONEMT_HPP
#define PARALLELENVIRONEMT_HPP

#ifdef USE_NVIDIA_GPU
#define VASPML_NV_HOST_DEVICE __host__ __device__
#else
#define VASPML_NV_HOST_DEVICE
#endif // USE_NVIDIA_GPU

#ifdef VASPML_USE_PSTL
#ifdef __INTEL_LLVM_COMPILER
#include <oneapi/dpl/execution>
#else
#include <execution>
#endif // __INTEL_LLVM_COMPILER
#endif // VASPML_USE_PSTL

namespace vaspml
{

#if VASPML_USE_PSTL
#define VASPML_PARALLEL(X) do { X } while( 0 )
#define VASPML_SEQ seq,
#define VASPML_PAR par,
#define VASPML_PAR_UNSEQ par_unseq,
#ifdef __INTEL_LLVM_COMPILER
constexpr auto seq = pstl::execution::seq;
constexpr auto par = pstl::execution::par;
constexpr auto par_unseq = pstl::execution::par_unseq;
#else
constexpr auto seq = std::execution::seq;
constexpr auto par = std::execution::par;
constexpr auto par_unseq = std::execution::par_unseq;
#endif
#else // VASPML_USE_PSTL
#define VASPML_PARALLEL(X)
#define VASPML_SEQ
#define VASPML_PAR
#define VASPML_PAR_UNSEQ
#endif // VASPML_USE_PSTL

enum class ExecutionPolicy
{
    cpuSingleCore,
    gpuStdLib,
};

} //namespace vaspml

#endif // PARALLELENVIRONEMT_HPP
