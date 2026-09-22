# 1 "elphon_kgrid.F"
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


# 2 "elphon_kgrid.F" 2 
!> This module implements symmetry handling using electron-phonon matrix elements
!> Formatter is applied:
!> fprettify -i 4 src/elphon_kgrid.F
module elphon_kgrid
    use prec
    use iso_c_binding, only: c_int, c_long, c_double, c_char
    use wave_struct_def, ONLY: wavespin, wavedes
    implicit none
    private
    integer, parameter :: bz_grid_points_estimation = 1000
    public elph_kgrid, elph_triplets
# 24


!>
!> @brief Data structure to handle reciprocal generalized regular grid.
!>
!> Smith normal form is given by D = P grid_matrix Q.
!>
!> In the following documentation, GR grid and BZ gird means generalized
!> regular grid and Brillouin (1._q,0._q) grid, respectively.
!>
!> One of BZ grid point index can be retrieved by
!>
!>     bz_grid_index = bz_grid_map(gr_grid_index + 1)
!>
!> Ir-grid point index is obtained from GR-grid point index by
!>
!>     ir_grid_index = ir_grid_indices(grg2irgrg(gr_grid_index + 1))
!>
!>
!> D_diag(3) : c_long
!>     Diagonal elements of D in Smith normal form, which correspond to mesh
!>     numbers in Q-transformed reciprocal primitive lattice.
!> P(3, 3) : c_long
!>     P matrix in Smith normal form. Transposed in fortran.
!> Q(3, 3) : c_long
!>     Q matrix in Smith normal form. Transposed in fortran.
!> QDinv(3, 3) : c_double
!>     QD^{-1}. Transposed in fortran.
!> grid_matrix(3, 3) : c_long
!>     Grid generating matrix with respect to reciprocal primitive basis
!>     vectors. Transposed in fortran.
!> regular_rotations(3, 3, 48) : c_long
!>     Rotation matrices with respect to reciprocal primitive basis
!>     vectors. q'=qR in fortran.
!> rotations(3, 3, 48) : c_long
!>     QD^-1 transformed rotation matrices with respect to transformed
!>     reciprocal primitive basis vectors. q'=qR in fortran.
!> PS(3) : c_long
!>     Shift in transformed double grid.
!> sampling_length : c_double
!>     Like that of KPOINTS auto R_k (Angst).
!> reciprocal_lattice(3, 3) : c_double
!>     Basis vectors of reciprocal primitive cell (Angst^-1). Without 2pi,
!>     row vectors in fortran.
!> microzone_lattice(3, 3) : c_double
!>     Basis vectors of microzone cell (Angst^-1). Without 2pi,
!>     row vectors in fortran.
!> n_rotations : c_long
!>     Number of point group operations. May or may not with time reversal
!>     symmetry.
!> n_ir_grid_points : c_long
!>     Number of irreducible grid points.
!> n_bz_grid_points : c_long
!>     Number of grid points in first BZ including its surface. This can be
!>     larger than product(D_diag).
!> grg2irgrg(product(D_diag)) : c_long
!>     Mapping table of GR-grid point index to index in ir_grid_indices.
!>     This array is created in this fortran module.
!>     This array is used only with ir_grid_indices.
!>     Index of ir_grid_indices stored starts with (1._q,0._q), i.e., fortran style.
!> ir_grid_indices(size(unique(grg2irgrg))) : c_long
!>     List of unique GR-grid indices.
!>     This array is created in this fortran module.
!>     Note that grid point index stored starts with (0._q,0._q).
!> ir_grid_weights(size(unique(grg2irgrg))) :  c_long
!>     Counts of unique grid indices in grg2irgrg. Sum = product(D_diag).
!>     This array is created in this fortran module.
!>     This array is used only with ir_grid_indices.
!> bz_grid_addresses(3, n_bz_grid_points) :  c_long
!>     Shortest grid points from Gamma point among translationally equivalent
!>     points in reciprocal lattice. All grid points having distances are
!>     stored.
!> bz_grid_map(product(D_diag) + 1) :  c_long
!>     List of first BZ-grid indices corresponding to GZ-grid points among
!>     translationally equivalent points.
!>     Note that grid point index stored starts with (0._q,0._q).
!> bzg2grg(n_bz_grid_points) :  c_long
!>     Mapping table of BZ-grid point index to GR-grid point index.
!>     Note that grid point index stored starts with (0._q,0._q).
!> relative_grid_addresses(3, 4, 24) :  c_long
!>     Relative vectors of neighboring grid points defined to use tetrahedron
!>     method on generalized regular grid.
!> transformation_matrix(3, 3) :  c_long
!>     Transpose of transformation matrix (T) from standardized conventional
!>     unit cell basis vectors to input cell basis vectors.
!>     It is assumed that the input cell is a primitive cell.
!>     (a_c, b_c, c_c) T = (a_p, b_p, c_p)
!> spacegroup_number :  c_long
!>     Space group number (from 1 to 230)
!> direct_rotations(3, 3, :) :  c_long
!>     Rotation matrices in direct space with respect to input cell basis
!>     vectors. The input cell is expected to be a primitive cell and
!>     number of rotations is equal to or smaller than 48.
!> n_direct_rotations :  c_long
!>     Number of direct rotations.
!>
!> by Atsushi Togo
!>
    type elph_kgrid
        integer(c_long) :: D_diag(3)
        integer(c_long) :: P(3, 3)
        integer(c_long) :: Q(3, 3)
        real(c_double) :: QDinv(3, 3)
        integer(c_long) :: grid_matrix(3, 3)
        integer(c_long) :: regular_rotations(3, 3, 48)
        integer(c_long) :: rotations(3, 3, 48)
        integer(c_long) :: PS(3)
        real(c_double) :: sampling_length
        real(c_double) :: reciprocal_lattice(3, 3)
        real(c_double) :: microzone_lattice(3, 3)
        integer(c_long) :: n_rotations
        integer(c_long) :: n_ir_grid_points
        integer(c_long) :: n_bz_grid_points
        integer(c_long), allocatable :: grg2irgrg(:)
        integer(c_long), allocatable :: ir_grid_indices(:)
        integer(c_long), allocatable :: ir_grid_weights(:)
        integer(c_long), allocatable :: bz_grid_addresses(:, :)
        integer(c_long), allocatable :: bz_grid_map(:)
        integer(c_long), allocatable :: bzg2grg(:)
        integer(c_long) :: relative_grid_addresses(3, 4, 24)
        real(c_double) :: transformation_matrix(3, 3)
        integer(c_long) :: spacegroup_number
        integer(c_long), allocatable :: direct_rotations(:, :, :)
        integer(c_long) :: n_direct_rotations
    end type elph_kgrid

!>
!> @brief Data structure for (q, k, k') triplets
!>
!> ir_triplets(3, num_ir_triplets) :
!>     Irreducible triplets at specified q in BZ-grid.
!> ir_weights(num_ir_triplets) :
!>     Integer weights of irreducible triplets at specified q.
!> num_ir_triplets
!>     Number of irreducible triplets at specified q.
!> rotation_indices
!>     Indices of rotations that map k to irreducible k
!>     The values are in {1..kgrid%n_rotations}.
!> rotated_triplets(3, num_ir_triplets) :
!>     {(Rq, Rk, Rk')} where Rk is in the list of kgrid%ir_grid_indices.
    type elph_triplets
        integer(c_long), allocatable :: ir_triplets(:, :)
        integer(c_long), allocatable :: ir_weights(:)
        integer(c_long) :: num_ir_triplets
        integer(c_long), allocatable :: rotation_indices(:)
        integer(c_long), allocatable :: rotated_triplets(:, :)
    end type elph_triplets

# 1468


end module elphon_kgrid
