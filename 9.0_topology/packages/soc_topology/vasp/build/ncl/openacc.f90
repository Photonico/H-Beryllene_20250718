# 1 "openacc.F"
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


# 2 "openacc.F" 2 
# 1666


!***********************************************************************
!
! OpenACC related Doxygen documentation page.
!
!***********************************************************************

!> @page openacc OpenACC
!! @tableofcontents
!! @section general All changes
!! All datatypes, subroutines, and functions that have been changed
!! to port VASP to GPU under OpenACC.
!!
!! \li      subrot::add_gamma_from_file
!! \li       hamil::eccp
!! \li       hamil::eccp_tau
!! \li       hamil::eccp_vec
!! \li       david::eddav
!! \li david_inner::eddav_inner
!! \li      subrot::eddiag
!! \li      subrot::eddiag_exact
!! \li    rmm_diis::eddrmm
!! \li       steep::edstep
!! \li            ::elmin
!! \li       hamil::hamilt_local
!! \li       hamil::hamilt_local_tau
!! \li       hamil::hamiltmu
!! \li       hamil::hamiltmu_tau
!! \li       hamil::hamiltmu_vec
!! \li            ::kinhamil
!! \li            ::kinhamil_tau
!! \li            ::kinhamil_vec
!! \li     locproj::lprj_proall
!! \li            ::potlok
!! \li morbitalmag::set_dd_magatom
!! \li        pawm::set_dd_paw
!! \li        meta::set_kineden
!! \li            ::setdij_
!! \li            ::setdij_avec_
!! \li     msphpro::sphpro_fast
!! \li morbitalmag::vectorpot
!! \li        wave::wvreal
!!
!! @section active Ported
!! \li       hamil::eccp
!! \li       hamil::eccp_tau
!! \li       hamil::eccp_vec
!! \li      subrot::eddiag
!! \li    rmm_diis::eddrmm
!! \li            ::elmin
!! \li       hamil::hamilt_local
!! \li       hamil::hamilt_local_tau
!! \li       hamil::hamiltmu
!! \li       hamil::hamiltmu_tau
!! \li       hamil::hamiltmu_vec
!! \li            ::kinhamil
!! \li            ::kinhamil_tau
!! \li            ::kinhamil_vec
!! \li        wave::wvreal
!!
!! @section inactive OpenACC deactivated
!! \li      subrot::add_gamma_from_file
!! \li       david::eddav
!! \li david_inner::eddav_inner
!! \li      subrot::eddiag_exact
!! \li       steep::edstep
!! \li     locproj::lprj_proall
!! \li         pot::potlok
!! \li morbitalmag::set_dd_magatom
!! \li        pawm::set_dd_paw
!! \li        meta::set_kineden
!! \li            ::setdij_
!! \li            ::setdij_avec_
!! \li     msphpro::sphpro_fast
!! \li morbitalmag::vectorpot
