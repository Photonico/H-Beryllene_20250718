# 1 "afqmc.F"
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


# 2 "afqmc.F" 2 
module afqmc

   use afqmc_struct
   use base, only: q, in_struct, symmetry, info_struct
   use hamil_struct_def, only: ham_handle
   use lattice, only: latt
   use tau_mu, only: tau_handle
   use mgrid, only: grid_3d, transit
   use mkpoints, only: kpoints_struct
   use nonl, only: nonl_struct
   use nonlr, only: nonlr_struct
   use poscar, only: type_info
   use pseudo, only: potcar
   use wave, only: wavedes, wavespin

   implicit none

   private
   public afqmc_settings, afqmc_reader, afqmc_propagation

contains

   subroutine afqmc_reader(incar, num_images, settings)
      use incar_reader, only: incar_file, process_incar
      use string, only: lowercase
      use tutor, only: vtutor
      type(incar_file), intent(inout) :: incar
      integer, intent(in) :: num_images
      type(afqmc_settings), intent(out) :: settings
      character(len=:), allocatable :: algo
!
      algo = ''
      call process_incar(incar, 'ALGO', algo)
      if (lowercase(trim(adjustl(algo))) == 'afqmc') then
          call vtutor%error("AFQMC is not implemented yet.")
      end if
!
   end subroutine afqmc_reader

   subroutine afqmc_propagation(settings, kineden, hamiltonian, p, w, wdes, nonlr_s, nonl_s, &
         latt_cur, t_info, info, io, grid, grid_soft, gridc, gridus, c_to_us, &
         soft_to_c, symm, chtot, cvtot, cstrf, chden, dencor, sv, cdij, cqij, &
         crhode, rholm, lmdim, irdmax, n_mix_paw, kpoints)
! parameters
      type(afqmc_settings), intent(in) :: settings
      type(tau_handle) kineden
      type(ham_handle) hamiltonian
      type(potcar) p(:)
      type(wavespin) w
      type(wavedes) wdes
      type(nonlr_struct) nonlr_s
      type(nonl_struct) nonl_s
      type(latt) latt_cur
      type(type_info) t_info
      type(info_struct) info
      type(in_struct) io
      type(grid_3d) grid, grid_soft, gridc, gridus
      type(transit) c_to_us, soft_to_c
      type(symmetry) symm
      type(kpoints_struct) kpoints
      complex(q) chtot(:,:), cvtot(:,:), cstrf(:,:), chden(:,:)
      COMPLEX(q) dencor(:), sv(:,:)
      COMPLEX(q) cdij(:,:,:,:), cqij(:,:,:,:), crhode(:,:,:,:)
      real(q) rholm(:,:)
      integer lmdim, irdmax, n_mix_paw
!
   end subroutine afqmc_propagation

end module afqmc
