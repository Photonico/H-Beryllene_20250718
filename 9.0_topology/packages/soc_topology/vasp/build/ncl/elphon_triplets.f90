# 1 "elphon_triplets.F"
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


# 2 "elphon_triplets.F" 2 

!> Module to deal with scattering events and their mapping
module triplets
  use prec, only: q
  use lattice, only: latt
  use mkpoints_struct_def, only: kpoints_struct
  implicit none

!> Contains all possible internal representations of (1._q,0._q) k-point
  type kpoint_type
    integer :: kindex !< integer index as computed by the hash table
    integer :: ibz !< integer index in the list of k-points
    integer :: kint(3) !< three intergers defining the k-point
!TODO: integer :: g(3) !< integer part of the k-point
  end type kpoint_type

!> Contains all the representations of a scattering event between 3 states
  type triplet_type
    type(kpoint_type) :: k
    type(kpoint_type) :: kp
    type(kpoint_type) :: q
  end type triplet_type
  contains

!> initialize the elph_kpoint structure from kindex
  pure subroutine kpoint_init_kindex(self, kpoints, kindex)
    use mkpoints, only: kindex_vkpt, kindex_kint, kint_vkpt, ikibz_kindex
    type(kpoint_type),intent(out) :: self
    type(kpoints_struct),intent(in) :: kpoints
    integer,intent(in) :: kindex
!
    self%kindex = kindex
    self%kint = kindex_kint(kpoints, self%kindex)
    self%ibz = kpoints%fbz2ibz(1, self%kindex)
!
  end subroutine kpoint_init_kindex

!> initialize the elph_kpoint structure from kibz
  pure subroutine kpoint_init_kibz(self, kpoints, latt_cur, kibz)
    use mkpoints, only: kindex_vkpt, kindex_kint, kint_vkpt, kint_kindex, kpoint_integers
    type(kpoint_type),intent(out) :: self
    type(kpoints_struct),intent(in) :: kpoints
    type(latt),intent(in) :: latt_cur
    integer,intent(in) :: kibz
!
    self%ibz = kibz
    call kpoint_integers(kpoints,latt_cur,kpoints%vkpt(:,self%ibz),self%kint)
    self%kindex = kint_kindex(kpoints,self%kint)
!
  end subroutine kpoint_init_kibz

!> initialize the elph_kpoint structure from vkpt
  pure subroutine kpoint_init_kint(self, kpoints, kint)
    use mkpoints, only: kpoint_integers, kindex_kint, kint_kindex, ikibz_kindex, kint_mod
    type(kpoint_type),intent(out) :: self
    type(kpoints_struct),intent(in) :: kpoints
    integer,intent(in) :: kint(3)
!
    self%kint = kint_mod(kpoints,kint)
    self%kindex = kint_kindex(kpoints,self%kint)
    self%ibz = kpoints%fbz2ibz(1,self%kindex)
!
  end subroutine kpoint_init_kint

!> initialize the elph_kpoint structure from vkpt
  pure subroutine kpoint_init_vkpt(self, kpoints, latt_cur, vkpt)
    use mkpoints, only: kpoint_integers, kindex_kint, kint_kindex, ikibz_kindex
    type(kpoint_type),intent(out) :: self
    type(kpoints_struct),intent(in) :: kpoints
    type(latt),intent(in) :: latt_cur
    real(q),intent(in) :: vkpt(3)
!
    call kpoint_integers(kpoints,latt_cur,vkpt,self%kint)
    self%kindex = kint_kindex(kpoints,self%kint)
    self%ibz = kpoints%fbz2ibz(1,self%kindex)
!
  end subroutine kpoint_init_vkpt

!> intialize the triplet_type from two integer indexes
  subroutine triplet_init_k_and_kp(self, kpoints, latt_cur, ik_fbz, ikp_fbz)
    use mkpoints, only: kindex_from_k_and_kp
    type(triplet_type),intent(out) :: self
    type(kpoints_struct),intent(in) :: kpoints
    type(latt),intent(in) :: latt_cur
    integer,intent(in) :: ik_fbz !< integer index in the hash table of kp
    integer,intent(in) :: ikp_fbz !< integer index in the hash table of iq
! local variables
    integer :: iq_fbz
    
    call kpoint_init_kindex(self%k, kpoints, ik_fbz)
    call kpoint_init_kindex(self%kp, kpoints, ikp_fbz)
    iq_fbz = kindex_from_k_and_kp(kpoints, latt_cur, ik_fbz, ikp_fbz)
    call kpoint_init_kindex(self%q, kpoints, iq_fbz)
    
  end subroutine triplet_init_k_and_kp

!> intialize the triplet_type from two integer indexes
  subroutine triplet_init_kp_and_q(self, kpoints, latt_cur, ikp_fbz, iq_fbz)
    use mkpoints, only: kindex_from_kp_and_q
    type(triplet_type),intent(out) :: self
    type(kpoints_struct),intent(in) :: kpoints
    type(latt),intent(in) :: latt_cur
    integer,intent(in) :: ikp_fbz !< integer index in the hash table of kp
    integer,intent(in) :: iq_fbz !< integer index in the hash table of iq
! local variables
    integer :: ik_fbz
    
    call kpoint_init_kindex(self%kp, kpoints, ikp_fbz)
    call kpoint_init_kindex(self%q, kpoints, iq_fbz)
    ik_fbz = kindex_from_kp_and_q(kpoints, latt_cur, ikp_fbz, iq_fbz)
    call kpoint_init_kindex(self%k, kpoints, ik_fbz)
    
  end subroutine triplet_init_kp_and_q

end module triplets
