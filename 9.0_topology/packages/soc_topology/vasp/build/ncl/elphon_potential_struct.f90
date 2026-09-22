# 1 "elphon_potential_struct.F"
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


# 2 "elphon_potential_struct.F" 2 

module elphon_potential_struct
    use prec
    use poscar_struct_def
    use mgrid_struct_def
    use pseudo, only: potcar
# 10

   implicit none

!> Contains a 3x3 S integer matrix and the Smith Normal decomposition Q D P:
!>     S = QM D PM
!> where D is a diagonal matrix
!> well as the inverse of the QM and PM matrices
    type snf
        integer :: s(3,3) !< Original matrix to be decomposed
        integer :: qm(3,3) !< Q matrix
        integer :: pm(3,3) !< P matrix
        integer :: d(3)    !< diagonal of the D matrix
        integer :: qinv(3,3) !< inverse of the Q matrix
        integer :: pinv(3,3) !< inverse of the P matrix
    end type snf

!> Data structure containing the force constants
    type phonon_ifc
        integer :: natoms_pc
        integer :: natoms_sc
        integer :: det_p2s

! must be set in init
        real(q), allocatable :: primitive_positions(:,:)
        real(q), allocatable :: supercell_positions(:,:)
        real(q), allocatable :: lattice_point(:,:)
        real(q) :: lattice_sp(3,3)
        real(q) :: lattice_pc(3,3)
        real(q) :: primitive_matrix(3,3)
        integer :: supercell_matrix(3,3)
        integer, allocatable :: atom_indices_in_derivatives(:)
        integer, allocatable :: ityp_pc(:)
        integer, allocatable :: ityp_sc(:)
! related to the long-range part
        real(q) :: dielectric_tensor(3,3)
        real(q), allocatable :: born_eff_charges(:,:,:)
        real(q), allocatable :: quadrupoles(:,:,:,:)
! related to the force constants
        real(q),allocatable :: phonon_masses(:)
        real(q),allocatable :: force_constants(:,:,:,:)
        real(q),allocatable :: shortest_vectors(:,:,:,:)
        real(q),allocatable :: shortest_vector_multiplicities(:,:)


! computed
        integer :: max_points_atoms
        real(q) :: prim2sp(3,3)
        real(q) :: sp2prim(3,3)
        real(q), allocatable :: centers(:,:)
        real(q), allocatable :: supercell_atoms(:,:)
        integer, allocatable :: supercell_atoms_map_to_phelel(:)
        integer, allocatable :: supercell_atoms_map_from_phelel(:)
! points for interpolation of the dynamical matrix
        real(q), allocatable :: sp_lattice_atoms(:,:,:)
        integer, allocatable :: sp_lattice_atoms_size(:,:,:)
        type (latt) :: latt_pc
        type (latt) :: latt_sp

! long range force constants
        logical :: has_lr !< true if we have information about the long-range part
        logical :: has_lr_active !< true if the long-range treatment is activated
        logical :: make_hermitian !< enforce the dynamical matrix to be hermitian
        integer :: ngvecs_pc !< number of g-vectors for the treatment of the long-range part in the unitcell
        integer :: ngvecs_pc_local !< number of g-vectors for the treatment of the long-range part in the unitcell
        integer :: ngvecs_sc !< number of g-vectors for the treatment of the long-range part in the supercell
!> Ewald parameter determining the separation between the reciprocal and real space
!> to evaluate the long-range part of the potential.
!> The evaluation of the real space part is not implemented since it is likely unecessary.
!> This parameter enters the supression argument which has the form:
!> exp( -(G.epsilon.G)/(4*ewald_param**2) )
        real(q) :: ewald_param
!> G vector cutoff for the evaluation of the long-range part
        real(q) :: encutlr
!> Long-range part of the force constants evaluated given a certain Ewald parameter
        real(q),allocatable :: force_constants_lr(:,:,:,:)
        integer,allocatable :: gvecs(:,:) !< coordinates of the G-vectors
        integer,allocatable :: nindpw(:) !< indexes to perform the FFT
    end type phonon_ifc

!> Data structure containing the electron-phonon potential
    type elphon_pot
        integer :: natoms_pc !< total number of atoms in the primitive cell
        integer :: natoms_pc_local !< local number of atoms in the primitive cell
        integer :: natoms_sc !< number of atoms in the supercell
        integer :: det_p2s !< number of primitive lattice cells that build the supercell
        integer :: ncdij !< number of spinor components of the potential
        integer :: lmdim !< max dimenion for the lm channels
        integer :: fft_mesh(3) !< Dimensions of the FFT grid
        integer :: npw !< total number of points in the FFT grid of the primitive cell
        integer :: npr_local !< local number of points in the FFT grid of the primitive cell (depends on ncore)
        integer :: ngridpoints !< total number of grid points in the supercell
        integer :: ngridpoints_local !< local number of grid points in the supercell

! settings
        logical :: skip !< allow to skip the computation of the potential
        logical :: has_shmem !< use shared memory to allocate the large arrays of the potential

! must be set in init
        real(q), allocatable :: primitive_positions(:,:)
        real(q), allocatable :: supercell_positions(:,:)
        real(q), allocatable :: lattice_point(:,:)
        real(q) :: lattice_sp(3,3)
        real(q) :: lattice_pc(3,3) ! xxx lattice vector of primitive cell as columns = [a1 a2 a2 ]
        real(q) :: primitive_matrix(3,3)
        integer :: supercell_matrix(3,3)
        integer, allocatable :: atom_indices_in_derivatives(:)
        integer, allocatable :: ityp_pc(:)
        integer, allocatable :: ityp_sc(:)
! related to the potential
        real(q), allocatable :: grid_point(:,:)
        complex(q), allocatable :: dvdu(:,:,:,:)
        complex(q), allocatable :: dvdu_longrange(:,:,:)
        complex(q), allocatable :: ccdij(:,:,:,:)
        complex(q), allocatable :: ccqij(:,:,:,:)
        complex(q), allocatable :: dqijdu(:,:,:,:,:,:)
        complex(q), allocatable :: ddijdu(:,:,:,:,:,:)
        complex(q), allocatable :: ccdrij_radial(:,:,:,:,:)
        complex(q), allocatable :: ccdrij(:,:,:,:,:)
! related to the long-range part
        real(q) :: dielectric_tensor(3,3)
        real(q), allocatable :: born_eff_charges(:,:,:)
        real(q), allocatable :: quadrupoles(:,:,:,:)
! related to the foce constants
        real(q),allocatable :: phonon_masses(:)

! computed
        integer :: max_points_atoms
        integer :: max_points_lattice
        real(q) :: prim2sp(3,3) ! xxx Transform the coordinates from primitive/current to supercell: [X1,X2,X3]=prim2sp*[x1,x2,x3]
        real(q) :: sp2prim(3,3) ! xxx sp2prim=prim2sp^-1; [A1 A2 A3] =  [a1 a2 a3] sp2prim
        real(q), allocatable :: centers(:,:)
        real(q), allocatable :: supercell_atoms(:,:)
        integer, allocatable :: supercell_atoms_map_to_phelel(:)
        integer, allocatable :: supercell_atoms_map_from_phelel(:)
! points for interpolation of dijdu
        real(q), allocatable :: sp_lattice_atoms(:,:,:)
        integer, allocatable :: sp_lattice_atoms_size(:,:,:)
! points for interpolation of the potential dvdu
        real(q), pointer :: sp_lattice_points(:,:,:)
        integer, allocatable :: sp_lattice_points_size(:,:,:)
        type (latt) :: latt_pc
        type (latt) :: latt_sp
!> grid of the primitive cell
        type (grid_3d) :: grid_pc
        type (grid_3d) :: grid_sc
!> comunicator for the grids
        type(communic) :: comm
!> comunicator for the atoms
        type(communic) :: comm_inb
!> comunicator between k-points
        type(communic) :: comm_kbatch
!> shared memory communicator
        type(communic) :: comm_shmem
!> shared memory segments
# 167

!long-range potential
        logical :: has_lr !< true if we have information about the long-range part
        logical :: has_lr_active !< true if the long-range treatment is activated
        integer :: ngvecs_pc !< number of g-vectors for the treatment of the long-range part in the unitcell
        integer :: ngvecs_pc_local !< number of g-vectors for the treatment of the long-range part in the unitcell
        integer :: ngvecs_sc !< number of g-vectors for the treatment of the long-range part in the supercell
!> Ewald parameter determining the separation between the reciprocal and real space
!> to evaluate the long-range part of the potential.
!> The evaluation of the real space part is not implemented since it is likely unecessary.
!> This parameter enters the supression argument which has the form:
!> exp( -(G.epsilon.G)/(4*ewald_param**2) )
        real(q) :: ewald_param
!> Radius of Coulomb truncation
        real(q) :: rcut
!> G vector cutoff for the evaluation of the long-range part
        real(q) :: encutlr
!> pointer to the pseudopotentials (needed for the LR correction due to the potential)
        type(potcar),pointer :: p(:)
!> Use slow fourier transform to compute the LR part in real space
        integer :: use_lr_sft
!> FFT mesh of the supercell
        integer :: fft_mesh_sc(3)
!> FFT mesh of the supercell prime lattice_spp = matmul(lattice_sp,snf_lat%qm)
        integer :: fft_mesh_scp(3)
!> Smith Normal Form decomposition of the matrix matrix relating
!> the FFT grid points in the primitive cell with the ones in the supercell
        type(snf) :: snf_fft
!> Smith Normal Form decomposition of the matrix matrix relating
!> the the primitive cell and the supercell
        type(snf) :: snf_lat
!> caches for the computation of the long-range
        logical,allocatable :: gmask(:) !< keep track of which g-vectors were set
        integer,allocatable :: gvecs(:,:) !< coordinates of the G-vectors
        integer,allocatable :: gvecs_local(:,:) !< coordinates of the G-vectors
        integer,allocatable :: nindpw_local(:) !< indexes to perform the FFT
    end type elphon_pot


end module elphon_potential_struct
