# 1 "build_info.F"
# 1 "./symbol.inc" 1 

!-------- to be costumized by user (usually done in the makefile)-------
!#define vector              compile for vector machine
!#define essl                use ESSL instead of LAPACK
!#define single_BLAS         use single prec. BLAS

!#define wNGXhalf            gamma only wavefunctions (X-red)
!#define wNGZhalf            gamma only wavefunctions (Z-red)

!#define NGXhalf             charge stored in REAL array (X-red)
!#define NGZhalf             charge stored in REAL array (Z-red)
!#define NOZTRMM             replace ZTRMM by ZGEMM
!#define 1                 compile for parallel machine with 1
!------------- end of user part --------------------------------
# 17

# 62

# 91

!
!   charge density: full grid mode
!










# 113

!
!   charge density complex
!






# 133

!
!   wavefunctions: full grid mode
!

# 182

!
!   wavefunctions complex
!








































!
!   common definitions
!








!
!   mpi parallel macros
!







# 268

# 274

# 279

# 286

# 293



# 302







# 317











# 339









!
! OpenMP macros
!
# 354

# 357

# 366









!
! profiling macros
!
# 381




!
! shmem macros
!
# 390

# 393

# 396

!
! quadruple precision
!
# 414











!
! for the 1 interface
!
# 430

!
! CUDA includes
!
# 438

!
! SIMD related definitions
!
# 1 "./simd.inc" 1 
!!#if   defined(__MIC__) || defined(__AVX512F__)
!!#define SIMD512
!!#undef  SIMD256
!!#elif defined(__AVX__) || defined(__AVX2__)
!!#define SIMD256
!!#undef  SIMD512
!!#endif

# 17





# 25





# 35















# 54





# 443 "./symbol.inc" 2 
!
! memalign macros
!
# 456




!
! Macros for PGI/NV HPC compilers version specific code
!
# 466

# 473



# 485










# 497


# 501


!
! OpenACC macros
!
# 527



































# 568












































!
! combined OpenMP and OpenACC macros
!
# 617



!
! routines replaced in LAPACK >=3.6
!
# 625


!
! Macros for HDF5 error check
!


# 634


!
! Macros for GNU version specific code
!
# 641


!
! For machine learning
!




!
! Extra safe initializations (+ overflow protections)
!
# 655




!
! Macros for memory estimation
!




!
! Offloading related macros
!
# 671


# 695

! Line is included when  !OFFLOADING

! Line is a comment when !OFFLOADING



!replace blas and lapack wrapper calls with
!normal calls if no offloading is used:





































# 747









# 759

! Line is a comment when !_OPENACC or   ACC_OFFLOAD



# 769

! Line is a comment when !ACC_OFFLOAD



# 779

! Line is included when  !_OPENACC

! Line is a comment when !_OPENACC



# 789

! Line is a comment when !_OPENMP or   OMP_OFFLOAD



# 799

! Line is a comment when !OMP_OFFLOAD



# 809

! Line is included when  !_OPENMP

! Line is a comment when !_OPENMP


# 2 "build_info.F" 2 

module build_info

   implicit none
   private

   public :: cpp_options, link_line

# 1 "./build_info.inc" 1 
    character(len=*), parameter :: cpp_options = '&
&-DHOST="LinuxIFC" &
&-DMPI &
&-DMPI_BLOCK=8000 &
&-Duse_collective &
&-DscaLAPACK &
&-DCACHE_SIZE=4000 &
&-Davoidalloc &
&-Dvasp6 &
&-Dtbdyn &
&-Dfock_dblbuf &
&-DVASP_HDF5 &
&-DVASP2WANNIER90'
    character(len=*), parameter :: link_line   = '&
&mpiifort &
&-fc=ifx &
&-mkl=sequential &
&-Llib &
&-ldmy &
&-Lparser &
&-lparser &
&-lstdc++ &
&-L/usr/physics/oneapi/2024.2.1/mkl/2024.2/lib/intel64 &
&-lmkl_scalapack_lp64 &
&-lmkl_blacs_intelmpi_lp64 &
&-L/usr/physics/hdf/5/1.14.1-2_intel2021/lib &
&-lhdf5_fortran &
&/cmt2/lniu6305/Packages/soc_topology_20260922_1435/wannier90-3.1.0/libwannier.a'
    character(len=*), parameter :: fc     = '&
&mpiifort &
&-fc=ifx'
    character(len=*), parameter :: fcl    = '&
&mpiifort &
&-fc=ifx &
&-mkl=sequential'
    character(len=*), parameter :: fflags = '&
&-assume &
&byterecl &
&-w &
&-axCORE-AVX512'
    character(len=*), parameter :: llibs  = '&
&-lstdc++ &
&-L/usr/physics/oneapi/2024.2.1/mkl/2024.2/lib/intel64 &
&-lmkl_scalapack_lp64 &
&-lmkl_blacs_intelmpi_lp64 &
&-L/usr/physics/hdf/5/1.14.1-2_intel2021/lib &
&-lhdf5_fortran &
&/cmt2/lniu6305/Packages/soc_topology_20260922_1435/wannier90-3.1.0/libwannier.a'
    character(len=*), parameter :: incs   = '&
&-I/usr/physics/oneapi/2024.2.1/mkl/2024.2/include/fftw &
&-I/usr/physics/hdf/5/1.14.1-2_intel2021/include'
# 11 "build_info.F" 2 

# 171

end module build_info
