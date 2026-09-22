# 1 "coulomb_cutoff_gradients.F"
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


# 2 "coulomb_cutoff_gradients.F" 2 

module coulomb_cutoff_gradients
!! use moffload_struct_def
!! use moffload
      use prec, only: q
      use mgrid, only: grid_3d, inilgrd, gen_rc_sub_grid, gen_rc_grid
      use lattice, only: latt, lattic, dirkar
      use pot_electrostat, only: potion
      use pot_struct_def
      use pseudo_struct_def, only: potcar
      use poscar_struct_def, only: type_info
      use mgrid_struct_def, only: transit
      use base, only: info_struct
      use constant, only: EDEPS, PI, TPI, FELECT, CITPI

      implicit none

contains
    subroutine truncate_forloc(gridc, p, t_info, latt_cur, &
                               coulomb_pot, charge_density, &
                               structure_factor, forces)
       type(grid_3d), intent(in) :: gridc
       type(latt), intent(in) :: latt_cur
       type(type_info), intent(in) :: t_info
       type(potcar), intent(in) :: p(t_info%ntyp)
       type(coulomb_potential), intent(in) :: coulomb_pot
       complex(q) :: charge_density(:)
       complex(q) :: structure_factor(:,:)
       real(q), intent(out) :: forces(:,:)
       integer :: i, nis, nt, n1, n2, n3, nc, ni, nion, ni_add
       complex(q) :: for
       real(q) :: for1, for2, for3
       real(q) :: dummy, factm
       complex(q) :: charge
       complex(q), allocatable :: ion_pot(:)

       
!
       allocate(ion_pot(gridc%MPLWV))
!$acc enter data copyin(ion_pot) 
!
       nis = 1
!
       do nt = 1, t_info%ntyp
!
       ni_add = t_info%NITYP(nt)
!
       do nion = nis, ni_add + nis - 1
!acc kernels present(ion_pot) 
           ion_pot = 0.0_q
!acc end kernels
!
           call potion(gridc, p, latt_cur, t_info, coulomb_pot, &
                       ion_pot, structure_factor, dummy, [nion])
!
           for1 = 0.0_q; for2 = 0.0_q; for3 = 0.0_q
!
!$acc parallel loop reduction(+:for1,for2,for3) &
!$acc private(n2,n3,ni,factm,for) &
!$acc present(gridc,charge_density,ion_pot)
           do nc = 1, gridc%RC%NCOL
      n2 = gridc%RC%I2(nc)
      n3 = gridc%RC%I3(nc)
           do n1 = 1, gridc%RC%NROW
!!          n2 = gridc%RC%I2(nc)
!!          n3 = gridc%RC%I3(nc)
!
               ni = (nc-1) * gridc%RC%NROW + n1
!
               FACTM=1
               
!
               for = cmplx(0.0_q, -1.0_q) * ion_pot(ni) &
                   *  conjg(charge_density(ni))
!
               for1 = for1 - gridc%LPCTX_(n1) * real(for,q)
               for2 = for2 - gridc%LPCTY_(n2) * real(for,q)
               for3 = for3 - gridc%LPCTZ_(n3) * real(for,q)
           end do
           end do
!
           forces(1,nion) = for1 * TPI
           forces(2,nion) = for2 * TPI
           forces(3,nion) = for3 * TPI
           enddo
!
           nis = nis + ni_add
!
       enddo
!
       CALL M_sum_d(gridc%comm, forces(1,1), t_info%nions*3)
!
       call dirkar(t_info%nions, forces, latt_cur%b)
!$acc exit data delete(ion_pot) 
!
       
    end subroutine
end module
