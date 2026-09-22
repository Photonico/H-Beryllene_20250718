# 1 "vhdf5.F"
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


# 2 "vhdf5.F" 2 


!> High level subroutines for writting data to the hdf5 files
module vhdf5
  use prec
  use poscar, only : type_info, dynamics, nh_chains
  use base, only : info_struct
  use wave, only: wavedes, wavespin
  use lattice, only : latt, dirkar
  use mkpoints, only : kpoints_struct
  use mkpoints_struct_def, only : skpoints_full
  use locproj_struct, only : LPRJ_function
  use fileio, only : closewav, write_from_buf
  use vhdf5_base

  implicit none

  integer,parameter :: MAX_LEN_GROUP=50

  public :: start_h5
  public :: stop_h5
  public :: open_vaspout
  public :: vh5_write_dosheader
  public :: vh5_write_eigenvalheader
  public :: vh5_write_eigenval
  public :: vh5_write_dos
  public :: vh5_write_efermi
  public :: vh5_write_lattice_ion
  public :: vh5_write_lattice_ion_hist_start
  public :: vh5_write_lattice_ion_hist_step
  public :: vh5_write_energies_hist_start
  public :: vh5_write_energies_hist_step
  public :: vh5_write_stress_forces_hist_start
  public :: vh5_write_stress_forces_hist_step
  public :: vh5_write_pair_correlation_start
  public :: vh5_write_pair_correlation_step
  public :: vh5_write_bandgap
  public :: vh5_write_positions_trailer
  public :: vh5_write_charge
  public :: vh5_write_paw_occupancies
  public :: vh5_write_paw_dij_qij
  public :: vh5_write_atommom
  public :: vh5_write_wavefunctions
  public :: vh5_read_wavefunction_header
  public :: vh5_read_wavefunctions
  public :: vh5_read_force_constants
  public :: vh5_write_projectors
  public :: vh5_write_locproj
  public :: vh5_dump_original_to_h5
  public :: vh5_write_oszicar_start
  public :: vh5_write_oszicar_hist_step

  interface vh5_write_dielectric_energies
    module procedure vh5_write_dielectric_energies_grid, vh5_write_dielectric_energies_mesh
  end interface vh5_write_dielectric_energies

contains

!> Initialize the HDF5 library and open vaspin.h5
  subroutine start_h5(input_file)
    character(len=*), intent(in), optional :: input_file  !< overwrite the default filename vaspin.h5
    character(len=:), allocatable :: filename
! vh5_start() initializes the hdf5 library and needs to be called before using any
! hdf5 related functionality
    if (present(input_file)) then
      filename = input_file
    else
      filename = VASPIN
    end if
    call vh5_error(vh5_start(),"vhdf5.F",71)
    if (hdf5_found) then
      call vh5_error(vh5_file_open_read(filename, ih5infileid),"vhdf5.F",73)
      call vh5_error(vh5_group_open(ih5infileid, grp_input, ih5ininputgroup_id),"vhdf5.F",74)
    end if
  end subroutine start_h5

!> Create the vaspout.h5 file for the HDF5 output
  subroutine open_vaspout(lsynch5_loc, subdir, output_file)
    logical, intent(in) :: lsynch5_loc  !< synchronize the HDF5 output to allow access with other processes
    character(len=*), intent(in), optional :: subdir  !< subdirectory in which the file is created
    character(len=*), intent(in), optional :: output_file  !< overwrite the default filename vaspout.h5
    character(len=:), allocatable :: filename
    if (present(output_file)) then
      filename = output_file
    else
      filename = VASPOUT
    end if
    if (present(subdir)) &
      filename = subdir // filename
    lsynch5 = lsynch5_loc
    call vh5_error(vh5_file_create_or_overwrite(filename, ih5outfileid, lsynch5),"vhdf5.F",92)
    call vh5_error(vh5_group_open_or_create(ih5outfileid, grp_input, ih5outinputgroup_id),"vhdf5.F",93)
    call vh5_error(vh5_group_open_or_create(ih5outfileid, grp_intermediate, ih5intermediategroup_id),"vhdf5.F",94)
    call vh5_write_version(ih5outfileid)
  end subroutine open_vaspout

  subroutine stop_h5(subdir)
    character(len=*), intent(in), optional :: subdir
    character(len=:), allocatable :: subdir_
    logical image_incar
    logical image_kpoints
!
! if no INCAR, try to close vaspin.h5
!
    if (.not. incar_found .and. hdf5_found) then
      call vh5_error(vh5_group_close(ih5ininputgroup_id),"vhdf5.F",107)
      call vh5_error(vh5_file_close(ih5infileid),"vhdf5.F",108)
    endif
!
! write original  input data to vaspout.h5
!
    if (present(subdir)) then
      subdir_ = subdir
    else
      subdir_ = ''
    end if
    if (poscar_found.or.incar_found.or.kpoints_found) then
      call vh5_error(vh5_group_open_or_create(ih5outfileid, grp_original, ih5outoriginalgroup_id),"vhdf5.F",119)
    endif
    if (poscar_found) then
      call vh5_error(vh5_dump_original_to_h5(ih5outoriginalgroup_id, subgrp_poscar, subdir_ // "POSCAR"),"vhdf5.F",122)
    endif
    if (kpoints_found) then
      inquire(file=subdir_ // 'KPOINTS', exist=image_kpoints)
      if (image_kpoints) then
        call vh5_error(vh5_dump_original_to_h5(ih5outoriginalgroup_id, subgrp_kpoints, subdir_ // "KPOINTS"),"vhdf5.F",127)
      else
        call vh5_error(vh5_dump_original_to_h5(ih5outoriginalgroup_id, subgrp_kpoints, "KPOINTS"),"vhdf5.F",129)
      end if
    endif
    if (incar_found) then
      inquire(file=subdir_ // 'INCAR', exist=image_incar)
      if (image_incar) then
        call vh5_error(vh5_dump_original_to_h5(ih5outoriginalgroup_id, subgrp_incar, subdir_ // "INCAR"),"vhdf5.F",135)
      else
        call vh5_error(vh5_dump_original_to_h5(ih5outoriginalgroup_id, subgrp_incar, "INCAR"),"vhdf5.F",137)
      end if
    endif
    if (poscar_found.or.incar_found.or.kpoints_found) then
      call vh5_error(vh5_group_close_writing(ih5outoriginalgroup_id),"vhdf5.F",141)
    endif
!
! close vaspout.h5
!
    call vh5_error(vh5_group_close_writing(ih5outinputgroup_id),"vhdf5.F",146)
    call vh5_error(vh5_group_close_writing(ih5intermediategroup_id),"vhdf5.F",147)
    call vh5_error(vh5_file_close_writing(ih5outfileid),"vhdf5.F",148)
! vh5_end() frees memory allocated by vh5_start() and must be called at the end
! of any program using hdf5
    call vh5_error(vh5_end(),"vhdf5.F",151)
  end subroutine stop_h5

!> @brief write Vasp version into hdf5 file
  subroutine vh5_write_version(fileid)
    use version, only: major, minor, patch
    integer(HID_T), intent(in) :: fileid
    integer(HID_T) :: groupid
    call vh5_error(vh5_group_open_or_create(fileid, grp_version, groupid),"vhdf5.F",159)
    call vh5_error(vh5_write(groupid, "major", major),"vhdf5.F",160)
    call vh5_error(vh5_write(groupid, "minor", minor),"vhdf5.F",161)
    call vh5_error(vh5_write(groupid, "patch", patch),"vhdf5.F",162)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",163)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",164)
  end subroutine vh5_write_version

!> @brief write header of doscar files.
  subroutine vh5_write_dosheader(fileid, typeinfo, jobpar, mdynamics, wave, lattice, info, kpoints, group, subgroup)
    integer, intent(in)              :: jobpar !< if io%lorbit > 10 then jobpar is 1
    integer(HID_T), intent(in)       :: fileid !< hdf5 file handle
    type(type_info), intent(in)      :: typeinfo !< typeinfo data structure (atom types, number of ions, etc.)
    type(dynamics), intent(in)       :: mdynamics !< MD data structure (atom positions, velocities, etc.)
    type(wavedes), intent(in)        :: wave !< wave function description data
    type(latt), intent(in)           :: lattice !< lattice description
    type(info_struct), intent(in)    :: info !< info structure
    type(kpoints_struct), intent(in) :: kpoints !< kpoints data
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_DOS)
!-------------------------------------------------------------------------------------
    real(q) :: aomega
    integer(HID_T) :: groupid, subgroupid
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
!-------------------------------------------------------------------------------------
    my_group    = GRP_RESULTS;  if (present(group))    my_group    = trim(group)
    my_subgroup = GRP_DOS;      if (present(subgroup)) my_subgroup = trim(subgroup)
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",186)
    call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",187)

    aomega = lattice%omega/typeinfo%nions
    call vh5_error(vh5_write(subgroupid, "nions", typeinfo%nions),"vhdf5.F",190)
    call vh5_error(vh5_write(subgroupid, "jobpar", jobpar),"vhdf5.F",191)
    call vh5_error(vh5_write(subgroupid, "aomega", aomega),"vhdf5.F",192)
    call vh5_error(vh5_write(subgroupid, "anorm", lattice%anorm*1E-10_q),"vhdf5.F",193)
    call vh5_error(vh5_write(subgroupid, "potim", mdynamics%potim*1e-15_q),"vhdf5.F",194)
    call vh5_error(vh5_write(subgroupid, "temperature", mdynamics%temp),"vhdf5.F",195)
    call vh5_error(vh5_write(subgroupid, "system", info%sznam1),"vhdf5.F",196)

    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",198)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",199)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",200)
  end subroutine vh5_write_dosheader

!> @brief write header of eigenval files.
  subroutine vh5_write_eigenvalheader(fileid, typeinfo, mdynamics, wave, lattice, info, kpoints, group, subgroup)
    integer(HID_T), intent(in)       :: fileid !< hdf5 file handle
    type(type_info), intent(in)      :: typeinfo !< typeinfo data structure (atom types, number of ions, etc.)
    type(dynamics), intent(in)       :: mdynamics !< MD data structure (atom positions, velocities, etc.)
    type(wavedes), intent(in)        :: wave !< wave function description data
    type(latt), intent(in)           :: lattice !< lattice description
    type(info_struct), intent(in)    :: info !< info structure
    type(kpoints_struct), intent(in) :: kpoints !< kpoints data
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_EIGENVAL)
!-------------------------------------------------------------------------------------
    real(q) :: aomega
    integer(HID_T) :: groupid, subgroupid
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
!-------------------------------------------------------------------------------------

    my_group = GRP_RESULTS;     if (present(group))    my_group = group
    my_subgroup = GRP_EIGENVAL; if (present(subgroup)) my_subgroup = subgroup
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",222)
    call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",223)

    aomega = lattice%omega/typeinfo%nions
    call vh5_error(vh5_write(subgroupid, "nions", typeinfo%nions),"vhdf5.F",226)
    call vh5_error(vh5_write(subgroupid, "nblocks", mdynamics%nblock * mdynamics%kblock),"vhdf5.F",227)
    call vh5_error(vh5_write(subgroupid, "ispin", wave%ispin),"vhdf5.F",228)
    call vh5_error(vh5_write(subgroupid, "aomega", aomega),"vhdf5.F",229)
    call vh5_error(vh5_write(subgroupid, "anorm", lattice%anorm*1E-10_q),"vhdf5.F",230)
    call vh5_error(vh5_write(subgroupid, "potim", mdynamics%potim*1e-15_q),"vhdf5.F",231)
    call vh5_error(vh5_write(subgroupid, "temperature", mdynamics%temp),"vhdf5.F",232)
    call vh5_error(vh5_write(subgroupid, "system", info%sznam1),"vhdf5.F",233)
    call vh5_error(vh5_write(subgroupid, "nelectrons", info%nelect),"vhdf5.F",234)
    call vh5_error(vh5_write(subgroupid, "kpoints", kpoints%nkpts),"vhdf5.F",235)
    call vh5_error(vh5_write(subgroupid, "nb_tot", wave%nb_tot),"vhdf5.F",236)

    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",238)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",239)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",240)
  end subroutine vh5_write_eigenvalheader

!> @brief write electronic energies
  subroutine vh5_write_eigenval(fileid, wave, wavesp, kpoints, kpoints_full, group, subgroup)
    integer(HID_T), intent(in)       :: fileid !< hdf5 file handle
    type(wavedes), intent(in)        :: wave !< wave function description data
    type(wavespin), intent(in)       :: wavesp !< wavefunctions including band index and spin
    type(kpoints_struct), intent(in) :: kpoints !< kpoints data
    type(skpoints_full), pointer, intent(in), optional :: kpoints_full
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_EIGENVAL)
!-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    logical :: write_kpoints_full
!-------------------------------------------------------------------------------------
    write_kpoints_full = .false.
    if (present(kpoints_full)) then
! make sure kpoints_full has been set up
      if (associated(kpoints_full)) write_kpoints_full = .true.
    endif
    my_group    = GRP_RESULTS;  if (present(group))    my_group    = trim(group)
    my_subgroup = GRP_EIGENVAL; if (present(subgroup)) my_subgroup = trim(subgroup)
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",264)
    call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",265)

    call vh5_error(vh5_write(subgroupid, "kpoint_coords", kpoints%vkpt),"vhdf5.F",267)
    call vh5_error(vh5_write(subgroupid, "kpoints_symmetry_weight", kpoints%wtkpt),"vhdf5.F",268)
    call vh5_error(vh5_write(subgroupid, "eigenvalues", real(wavesp%celtot)),"vhdf5.F",269)
    call vh5_error(vh5_write(subgroupid, "fermiweights", wavesp%fertot),"vhdf5.F",270)
! write tetrahedra information if available
    if (kpoints%ltet) then
       call vh5_error(vh5_write(subgroupid, 'num_tetrahedra', kpoints%ntet),"vhdf5.F",273)
       call vh5_error(vh5_write(subgroupid, 'volume_weight_tetrahedra', kpoints%volwgt),"vhdf5.F",274)
       call vh5_error(vh5_write(subgroupid, 'coordinate_id_tetrahedra', kpoints%idtet(:,1:kpoints%ntet)),"vhdf5.F",275)
     endif
! write full kpoint mesh
    if (write_kpoints_full) then
      call vh5_error(vh5_write(subgroupid, "kpoint_coords_full", kpoints_full%vkpt),"vhdf5.F",279)
      call vh5_error(vh5_write(subgroupid, "kpoints_symmetry_weight_full", kpoints_full%wtkpt),"vhdf5.F",280)
      call vh5_error(vh5_write(subgroupid, "kpoints_symmetry_symop", kpoints_full%isymop),"vhdf5.F",281)
      call vh5_error(vh5_write(subgroupid, "kpoints_symmetry_mapping", kpoints_full%nequiv),"vhdf5.F",282)
    endif

    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",285)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",286)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",287)
  end subroutine vh5_write_eigenval

!> @brief write electronic energies
  subroutine vh5_write_band_velocities(fileid, velocities, group, subgroup)
    integer(HID_T), intent(in) :: fileid !< hdf5 file handle
    real(q), intent(in) :: velocities(:,:,:,:) !< band velocities
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_EIGENVAL)
!-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
!-------------------------------------------------------------------------------------
    my_group    = GRP_RESULTS;  if (present(group))    my_group    = trim(group)
    my_subgroup = GRP_EIGENVAL; if (present(subgroup)) my_subgroup = trim(subgroup)
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",302)
    call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",303)

    call vh5_error(vh5_write(subgroupid, "velocities", velocities),"vhdf5.F",305)

    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",307)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",308)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",309)
  end subroutine vh5_write_band_velocities

!> @brief write phonon frequencies and eigenvectors
  subroutine vh5_write_phonons(fileid, prim_lat, prim_pos, prim_nityp, prim_ityp, q_mesh, frequ, evecs, group, subgroup)
    use tutor, ONLY: vtutor
    integer(HID_T), intent(in)        :: fileid !< hdf5 file handle
    real(q), intent(in)               :: prim_lat(3,3) !< lattice vectors of the primitive cell
    real(q), intent(in)               :: prim_pos(:,:) !< primitive positions within the primitive cell
    integer, intent(in)               :: prim_nityp(:) !< number of atoms of each type in the primitive cell
    character(len=*), intent(in)      :: prim_ityp(:)  !< string identifier of each of the atomic types
    type (kpoints_struct), intent(in) :: q_mesh
    real(q), intent(in) :: frequ(:, :)
    complex(q), intent(in) :: evecs(:, :, :)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_EIGENVAL)
!-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    integer :: nmodes, nqpoints, nions
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
!-------------------------------------------------------------------------------------
    my_group    = GRP_RESULTS;  if (present(group))    my_group    = trim(group)
    my_subgroup = GRP_PHONONS; if (present(subgroup)) my_subgroup = trim(subgroup)
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",332)
    call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",333)

    call vh5_error(vh5_write(subgroupid, "primitive/scale", 1.0_q),"vhdf5.F",335)
    call vh5_error(vh5_write(subgroupid, "primitive/lattice_vectors", prim_lat),"vhdf5.F",336)
    call vh5_error(vh5_write(subgroupid, "primitive/position_ions", prim_pos),"vhdf5.F",337)
    call vh5_error(vh5_write(subgroupid, "primitive/number_ion_types", prim_nityp),"vhdf5.F",338)
    call vh5_error(vh5_write(subgroupid, "primitive/ion_types", prim_ityp),"vhdf5.F",339)

    call vh5_error(vh5_write(subgroupid, "qpoint_coords", q_mesh%vkpt),"vhdf5.F",341)
    call vh5_error(vh5_write(subgroupid, "qpoints_symmetry_weight", q_mesh%wtkpt),"vhdf5.F",342)

    if (q_mesh%nkpts/=size(frequ,2)) call vtutor%bug('Inconsistent dimension for frequencies array',"vhdf5.F",344)
    if (q_mesh%nkpts/=size(evecs,3)) call vtutor%bug('Inconsistent dimension for frequencies array',"vhdf5.F",345)
    nmodes   = size(frequ,1)
    nqpoints = size(frequ,2)
    nions = nmodes/3
    call vh5_error(vh5_write(subgroupid, "nions", nions),"vhdf5.F",349)
    call vh5_error(vh5_write(subgroupid, "frequencies", frequ),"vhdf5.F",350)
    call vh5_error(vh5_write(subgroupid, "eigenvectors", reshape(evecs,[3,nions,nmodes,nqpoints])),"vhdf5.F",351)

    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",353)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",354)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",355)
  end subroutine vh5_write_phonons

!> Write phonon density of states
  subroutine vh5_write_phonon_dos(fileid, energies, phonon_dos, group, subgroup)
    integer(HID_T), intent(in)        :: fileid !< hdf5 file handle
    real(q), intent(in) :: energies(:)
    real(q), intent(in) :: phonon_dos(:, :, :)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_EIGENVAL)
!-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    integer :: nmodes
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
!-------------------------------------------------------------------------------------
    my_group    = GRP_RESULTS; if (present(group))    my_group    = trim(group)
    my_subgroup = GRP_PHONONS; if (present(subgroup)) my_subgroup = trim(subgroup)
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",372)
    call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",373)

    call vh5_error(vh5_write(subgroupid, "dos_mesh", energies),"vhdf5.F",375)
    call vh5_error(vh5_write(subgroupid, "dos", sum(sum(phonon_dos,dim=2),dim=2)),"vhdf5.F",376)
    call vh5_error(vh5_write(subgroupid, "dospar", phonon_dos),"vhdf5.F",377)
    call vh5_error(vh5_write(subgroupid, "directions", ["x","y","z"]),"vhdf5.F",378)

    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",380)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",381)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",382)
  end subroutine vh5_write_phonon_dos

!> @brief Same as above but simpler since it does not use the q_data structure
  subroutine vh5_write_phonons_simple(fileid, vkpt, eigenvalues, eigenvectors,k_path_length, group, subgroup)
    integer(HID_T), intent(in)        :: fileid !< hdf5 file handle
    real(q),intent(in) :: vkpt(:,:)
    real(q),intent(in) :: eigenvalues(:,:)
    complex(q),intent(in) :: eigenvectors(:,:,:)
    real(q),intent(in) :: k_path_length(:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_EIGENVAL)
!-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    integer :: i, nmodes, nqpoints
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    integer(SIZE_T) :: start2(2), count2(2), start3(3), count3(3)
!-------------------------------------------------------------------------------------
    my_group    = GRP_RESULTS;  if (present(group))    my_group    = trim(group)
    my_subgroup = GRP_PHONONS; if (present(subgroup)) my_subgroup = trim(subgroup)
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",402)
    call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",403)

    nqpoints = size(vkpt,2)
    nmodes = size(eigenvectors,1)
    call vh5_error(vh5_write(subgroupid, "kpoint_coords", vkpt),"vhdf5.F",407)
    call vh5_error(vh5_write(subgroupid, "k_path_length", k_path_length),"vhdf5.F",408)

    call vh5_error(vh5_write(subgroupid, "nions", nmodes/3),"vhdf5.F",410)
    call vh5_error(vh5_write(subgroupid, "eigenvalues", eigenvalues),"vhdf5.F",411)
    call vh5_error(vh5_write(subgroupid, "eigenvectors", eigenvectors),"vhdf5.F",412)

    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",414)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",415)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",416)
  end subroutine vh5_write_phonons_simple

!> @brief write electronic density of states
  subroutine vh5_write_dos(fileid, wave, kpoints, dos, dosi, dospar, efermi, nionp, lpar, group, subgroup, scale)
    integer(HID_T), intent(in)       :: fileid !< hdf5 file handle
    type(wavedes), intent(in)        :: wave !< wave function description data
    type(kpoints_struct), intent(in) :: kpoints !< kpoints data
    real(q), intent(in) :: dos(:,:) !< density of states
    real(q), intent(in) :: dosi(:,:) !< integrated dos
    real(q), intent(in) :: dospar(:,:,:,:) !< partial density of states
    real(q), intent(in) :: efermi !< fermi energy
    integer, intent(in) :: nionp !< number of ions
    integer, intent(in) :: lpar
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_DOS)
    real(q), intent(in),optional :: scale !< scale dos and idos before writing
!-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    integer :: i, nedos
    real(q) :: my_scale, deltae
    real(q) :: energies(size(dos,1))
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
!-------------------------------------------------------------------------------------

    my_group    = GRP_RESULTS;  if (present(group))    my_group    = trim(group)
    my_subgroup = GRP_DOS;      if (present(subgroup)) my_subgroup = trim(subgroup)
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",443)
    call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",444)

    my_scale    = 1.d0; if (present(scale)) my_scale = scale
    nedos = size(dos,1)

    call vh5_error(vh5_write(subgroupid, "nionp", nionp),"vhdf5.F",449)
    call vh5_error(vh5_write(subgroupid, "ncdij", wave%ncdij),"vhdf5.F",450)
    call vh5_error(vh5_write(subgroupid, "lpar", lpar),"vhdf5.F",451)

    call vh5_error(vh5_write(subgroupid, "emax", kpoints%emax),"vhdf5.F",453)
    call vh5_error(vh5_write(subgroupid, "emin", kpoints%emin),"vhdf5.F",454)
    call vh5_error(vh5_write(subgroupid, "nedos", nedos),"vhdf5.F",455)
    call vh5_error(vh5_write(subgroupid, "efermi", efermi),"vhdf5.F",456)

    deltae=(kpoints%emax-kpoints%emin)/(nedos-1)
    do i=1,nedos
        energies(i)=kpoints%emin+deltae*(i-1)
    end do
    call vh5_error(vh5_write(subgroupid, "energies", energies),"vhdf5.F",462)
    call vh5_error(vh5_write(subgroupid, "dos", dos*my_scale),"vhdf5.F",463)
    call vh5_error(vh5_write(subgroupid, "dosi", dosi*my_scale),"vhdf5.F",464)
    call vh5_error(vh5_write(subgroupid, "dospar", dospar),"vhdf5.F",465)

    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",467)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",468)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",469)
  end subroutine vh5_write_dos

!> @brief: small helper to only update the fermi level in the h5 archive
  subroutine vh5_write_efermi(fileid, efermi, group, subgroup)
    integer(HID_T), intent(in)       :: fileid !< hdf5 file handle
    real(q), intent(in) :: efermi !< fermi energy
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_DOS)
!-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
!-------------------------------------------------------------------------------------

    my_group    = GRP_RESULTS;  if (present(group))    my_group    = trim(group)
    my_subgroup = GRP_DOS;      if (present(subgroup)) my_subgroup = trim(subgroup)
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",485)
    call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",486)

    call vh5_error(vh5_write(subgroupid, "efermi", efermi),"vhdf5.F",488)

    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",490)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",491)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",492)
  end subroutine vh5_write_efermi

!> @brief write lattice parameters and positions
  subroutine vh5_write_lattice_ion(fileid, sznam, typeinfo, scale, amat, selective_dynamics, posions, group, subgroup)
    integer(HID_T), intent(in)    :: fileid !< hdf5 file handle
    character(LEN=*), intent(in)  :: sznam !< string with system name read from poscar
    type(type_info), intent(in)   :: typeinfo !< typeinfo data structure (atom types, number of ions, etc.)
    real(q), intent(in) :: scale !< scale for the lattice vectors amat
    real(q), intent(in) :: amat(3,3) !< lattice vectors
    logical :: selective_dynamics !< logical determining whether we are running selective dynamics
    real(q), dimension(3,typeinfo%nions), intent(in) :: posions !< ion positions
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_POSITIONS)
!-------------------------------------------------------------------------------------
    integer(HID_T) :: posgroupid
    integer(HID_T) :: unitgroupid
    character(LEN=40) :: my_group, my_subgroup
!-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;      if (present(group))    my_group = group
    my_subgroup = GRP_POSITIONS; if (present(subgroup)) my_subgroup = subgroup
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), unitgroupid),"vhdf5.F",513)
    call vh5_error(vh5_group_open_or_create(unitgroupid, trim(my_subgroup), posgroupid),"vhdf5.F",514)

! lattice and system information
    call vh5_error(vh5_write(posgroupid, "system", sznam),"vhdf5.F",517)
    call vh5_error(vh5_write(posgroupid, "scale", scale),"vhdf5.F",518)
    call vh5_error(vh5_write(posgroupid, "lattice_vectors", amat/scale),"vhdf5.F",519)
    call vh5_error(vh5_write(posgroupid, "direct_coordinates", 1),"vhdf5.F",520)
    call vh5_error(vh5_write(posgroupid, "ion_types", typeinfo%type(1:typeinfo%ntyp)),"vhdf5.F",521)
    call vh5_error(vh5_write(posgroupid, "ion_sha256", typeinfo%sha256(1:typeinfo%ntyp)),"vhdf5.F",522)
    call vh5_error(vh5_write(posgroupid, "number_ion_types", typeinfo%nityp(1:typeinfo%ntyp)),"vhdf5.F",523)

!positions of the ions
    call vh5_error(vh5_write(posgroupid, "position_ions", posions),"vhdf5.F",526)

!selective dynamics
    call vh5_error(vh5_write(posgroupid, "selective_dynamics", selective_dynamics),"vhdf5.F",529)
    if (selective_dynamics) then
       call vh5_error(vh5_write(posgroupid, "selective_dynamics_ions", typeinfo%lsfor),"vhdf5.F",531)
    end if

    call vh5_error(vh5_group_close_writing(posgroupid),"vhdf5.F",534)
    call vh5_error(vh5_group_close_writing(unitgroupid),"vhdf5.F",535)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",536)
  end subroutine vh5_write_lattice_ion

!> @brief start writting of lattice vectors and ion positions history to hdf5 file
  subroutine vh5_write_lattice_ion_hist_start(fileid, subgroup, nions, dyn)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    integer,intent(in) :: nions
    type(dynamics), intent(in) :: dyn
!-------------------------------------------------------------
    integer(HID_T) :: groupid
    integer(SIZE_T) :: dims(3), maxdims(3), chunk_size(3)
!-------------------------------------------------------------
    call vh5_error(vh5_group_open_or_create(fileid, trim(subgroup), groupid),"vhdf5.F",549)

    dims = [3,nions,1]
    maxdims = [integer(SIZE_T)::3,nions,H5S_UNLIMITED_F]
    chunk_size = [1,1,min(vh5_default_chunk_size, max(dyn%nsw / dyn%nblock, 1))]
    call vh5_error(vh5_create_double_array_nd(groupid, "position_ions", 3, dims, maxdims, chunk_size),"vhdf5.F",554)

    dims = [3,3,1]; maxdims = [integer(SIZE_T)::3,3,H5S_UNLIMITED_F]
    call vh5_error(vh5_create_double_array_nd(groupid, "lattice_vectors", 3, dims, maxdims),"vhdf5.F",557)

    call vh5_error(vh5_write(groupid, "scale", 1.0_q),"vhdf5.F",559)

    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",561)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",562)
  end subroutine vh5_write_lattice_ion_hist_start

!> @brief write the lattice vectors and ion positions for (1._q,0._q) ionic step to hdf5 file
  subroutine vh5_write_lattice_ion_hist_step(fileid, subgroup, nstep, amat, posions)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    integer, intent(in)          :: nstep
    real(q), intent(in) :: amat(3,3)
    real(q), intent(in) :: posions(:,:)
!-------------------------------------------------------------------------------------
    integer :: nions
    integer(HID_T) :: groupid
    integer(SIZE_T) :: start(3), count(3), extend(3)
!-------------------------------------------------------------------------------------
    call vh5_error(vh5_group_open_or_create(fileid, trim(subgroup), groupid),"vhdf5.F",577)

    nions = size(posions,2)
    start = [1,1,nstep]
    count = [3,nions,1]
    extend = [0,0,1]
    call vh5_error(vh5_write_double_subarray_nd(groupid, "position_ions", 3, start, count, posions, extend),"vhdf5.F",583)

    count = [3,3,1]
    call vh5_error(vh5_write_double_subarray_nd(groupid, "lattice_vectors", 3, start, count, amat, extend),"vhdf5.F",586)

    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",588)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",589)
  end subroutine vh5_write_lattice_ion_hist_step

!> create entry in HDF5 file to store the ion velocities
  subroutine vh5_write_velocity_hist_start(fileid, subgroup, nions, dyn)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    integer,intent(in) :: nions
    type(dynamics), intent(in) :: dyn
!-------------------------------------------------------------
    integer(HID_T) :: groupid
    integer(SIZE_T) :: dims(3), maxdims(3), chunk_size(3)
!-------------------------------------------------------------
    call vh5_error(vh5_group_open_or_create(fileid, trim(subgroup), groupid),"vhdf5.F",602)

    dims = [3,nions,1]
    maxdims = [integer(SIZE_T)::3,nions,H5S_UNLIMITED_F]
    chunk_size = [1,1,min(vh5_default_chunk_size, max(dyn%nsw / dyn%nblock, 1))]
    call vh5_error(vh5_create_double_array_nd(groupid, "ion_velocities", 3, dims, maxdims, chunk_size),"vhdf5.F",607)

    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",609)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",610)
  end subroutine vh5_write_velocity_hist_start

!> write the ion velocities for (1._q,0._q) ionic step to hdf5 file
  subroutine vh5_write_velocity_hist_step(fileid, subgroup, nstep, lattice, dyn)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    integer, intent(in)          :: nstep
    type(latt), intent(in) :: lattice
    type(dynamics), intent(in) :: dyn
!-------------------------------------------------------------------------------------
    integer :: nions
    integer(HID_T) :: groupid
    integer(SIZE_T) :: start(3), count(3), extend(3)
    real(q), allocatable :: velocities(:,:)
!-------------------------------------------------------------------------------------
    call vh5_error(vh5_group_open_or_create(fileid, trim(subgroup), groupid),"vhdf5.F",626)

! convert velocities to A/fs
    allocate(velocities, source=dyn%vel)
    if (dyn%ibrion /= 44 .and. dyn%ibrion /= 40) then
       do nions = 1, size(velocities, 2)
          velocities(:,nions) = velocities(:,nions) / dyn%potim
          call dirkar(1,velocities(:,nions),lattice%a)
       end do
    end if

    nions = size(velocities,2)
    start = [1,1,nstep]
    count = [3,nions,1]
    extend = [0,0,1]
    call vh5_error(vh5_write_double_subarray_nd(groupid, "ion_velocities", 3, start, count, velocities, extend),"vhdf5.F",641)

    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",643)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",644)
  end subroutine vh5_write_velocity_hist_step

!> @brief start writing of enegies history to hdf5 file
  subroutine vh5_write_energies_hist_start(fileid, subgroup, ibrion)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    integer, intent(in) :: ibrion
!-------------------------------------------------------------
    integer :: extent
    integer(HID_T) :: groupid
    integer(SIZE_T) :: dims(2), maxdims(2)
    character(LEN=24) :: energies_tags(7)
!-------------------------------------------------------------
    call vh5_error(vh5_group_open_or_create(fileid, trim(subgroup), groupid),"vhdf5.F",658)

    if (ibrion == 0) then
        energies_tags = [ 'ion-electron   TOTEN   ',&
                          'kinetic energy EKIN    ',&
                          'kin. lattice   EKIN_LAT',&
                          'temperature    TEIN    ',&
                          'nose potential ES      ',&
                          'nose kinetic   EPS     ',&
                          'total energy   ETOTAL  ' ]
        call vh5_error(vh5_write(groupid, "energies_tags", energies_tags),"vhdf5.F",668)
        extent = 7
    else
        energies_tags(1:3) = [ 'free energy    TOTEN   ',&
                               'energy without entropy ',&
                               'energy(sigma->0)       ']
        call vh5_error(vh5_write(groupid, "energies_tags", energies_tags(1:3)),"vhdf5.F",674)
        extent = 3
    endif

! Create an extensible array that will contain all the energies
    dims = [extent,1]
    maxdims = [integer(SIZE_T)::extent,H5S_UNLIMITED_F]
    call vh5_error(vh5_create_double_array_nd(groupid, "energies", 2, dims, maxdims),"vhdf5.F",681)

    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",683)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",684)
  end subroutine vh5_write_energies_hist_start

!> @brief write the energies for (1._q,0._q) ionic step to hdf5 file
  subroutine vh5_write_energies_hist_step(fileid, subgroup, nstep, energies)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    integer, intent(in)          :: nstep
    real(q), intent(in) :: energies(:)
!-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid
    integer(SIZE_T) :: start(2), count(2), extend(2)
!-------------------------------------------------------------------------------------
    call vh5_error(vh5_group_open_or_create(fileid, trim(subgroup), groupid),"vhdf5.F",697)

    start = [1,nstep]
    count = [size(energies),1]
    extend = [0,1]
    call vh5_error(vh5_write_double_subarray_nd(groupid, "energies", 2, start, count, energies, extend),"vhdf5.F",702)

    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",704)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",705)
  end subroutine vh5_write_energies_hist_step

!> @brief write bse index
  subroutine vh5_write_bseindex(fileid, bseindex, group, subgroup)
    integer(HID_T), intent(in)  :: fileid
    integer,intent(in) :: bseindex(:,:,:,:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
!-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
!-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",720)
    call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",721)

    call vh5_error(vh5_write(subgroupid, "bse_index", bseindex),"vhdf5.F",723)

    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",725)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",726)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",727)
  end subroutine vh5_write_bseindex

!> @brief write bands included in bse and the Fermi energy
  subroutine vh5_write_bsebands(fileid, bsebands, efermi, group, subgroup)
    use bse_struct, only : banddesc
    integer(HID_T), intent(in)  :: fileid
    type(banddesc), intent(in) :: bsebands(:)
    real(q), intent(in), optional :: efermi
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
!-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    integer :: isp
    character(LEN=40) :: my_group, my_subgroup
!-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",745)
    call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",746)

    call vh5_error(vh5_write(subgroupid, "bse_vbmax", (/ (bsebands(isp)%vbmax, isp=1,size(bsebands)) /)),"vhdf5.F",748)
    call vh5_error(vh5_write(subgroupid, "bse_vbmin", (/ (bsebands(isp)%vbmin, isp=1,size(bsebands)) /)),"vhdf5.F",749)
    call vh5_error(vh5_write(subgroupid, "bse_cbmin", (/ (bsebands(isp)%cbmin, isp=1,size(bsebands)) /)),"vhdf5.F",750)
    call vh5_error(vh5_write(subgroupid, "bse_cbmax", (/ (bsebands(isp)%cbmax, isp=1,size(bsebands)) /)),"vhdf5.F",751)
    if (present(efermi)) then
       call vh5_error(vh5_write(subgroupid, "efermi", efermi),"vhdf5.F",753)
    end if

    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",756)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",757)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",758)
  end subroutine vh5_write_bsebands

!> @brief write bse fatbands
  subroutine vh5_write_bsefatband(fileid, fatbands, group, subgroup)
    integer(HID_T), intent(in)  :: fileid
    complex(q),intent(in) :: fatbands(:,:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
!-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
!-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",773)
    call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",774)

    call vh5_error(vh5_write(subgroupid, "bse_fatbands", fatbands),"vhdf5.F",776)

    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",778)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",779)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",780)
  end subroutine vh5_write_bsefatband

!> @brief start writting of stress and forces history to hdf5 file
  subroutine vh5_write_stress_forces_hist_start(fileid, subgroup, nions, dyn)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    integer,intent(in) :: nions
    type(dynamics), intent(in) :: dyn
!-------------------------------------------------------------
    integer(HID_T) :: posgroupid
    integer(SIZE_T) :: dims(3), maxdims(3), chunk_size(3)
!-------------------------------------------------------------
    call vh5_error(vh5_group_open_or_create(fileid, trim(subgroup), posgroupid),"vhdf5.F",793)

    dims = [3,nions,1]
    maxdims = [integer(SIZE_T)::3,nions,H5S_UNLIMITED_F]
    chunk_size = [1,1,min(vh5_default_chunk_size, max(dyn%nsw / dyn%nblock, 1))]
    call vh5_error(vh5_create_double_array_nd(posgroupid, "forces", 3, dims, maxdims, chunk_size),"vhdf5.F",798)

    dims = [3,3,1]
    maxdims = [integer(SIZE_T)::3,3,H5S_UNLIMITED_F]
    call vh5_error(vh5_create_double_array_nd(posgroupid, "stress", 3, dims, maxdims),"vhdf5.F",802)

    call vh5_error(vh5_group_close_writing(posgroupid),"vhdf5.F",804)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",805)
  end subroutine vh5_write_stress_forces_hist_start

!> @brief write the stress and the forces for (1._q,0._q) ionic step to hdf5 file
  subroutine vh5_write_stress_forces_hist_step(fileid, subgroup, nstep, stress, forces)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    integer, intent(in)          :: nstep
    real(q), intent(in) :: stress(3,3)
    real(q), intent(in) :: forces(:,:)
!-------------------------------------------------------------------------------------
    integer :: nions
    integer(HID_T) :: posgroupid
    integer(SIZE_T) :: start(3), count(3), extend(3)
!-------------------------------------------------------------------------------------
    call vh5_error(vh5_group_open_or_create(fileid, trim(subgroup), posgroupid),"vhdf5.F",820)

    nions = size(forces,2)
    start = [1,1,nstep]
    count = [3,nions,1]
    extend = [0,0,1]
    call vh5_error(vh5_write_double_subarray_nd(posgroupid, "forces", 3, start, count, forces, extend),"vhdf5.F",826)

    count = [3,3,1]
    call vh5_error(vh5_write_double_subarray_nd(posgroupid, "stress", 3, start, count, stress, extend),"vhdf5.F",829)

    call vh5_error(vh5_group_close_writing(posgroupid),"vhdf5.F",831)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",832)
  end subroutine vh5_write_stress_forces_hist_step

!> @brief start writing output of OSZICAR to hdf5 file and store constants
  subroutine vh5_write_oszicar_start(fileid, subgroup)
      integer(HID_T), intent(in) :: fileid
      character(len=*), intent(in) :: subgroup
!
      integer(HID_T) :: posgroupid
      integer(SIZE_T) :: dims(2), maxdims(2)
      character(len=6), allocatable :: label(:)

      call vh5_error(vh5_group_open_or_create(fileid, trim(subgroup), posgroupid),"vhdf5.F",844)
      dims = [7, 1]; maxdims = [integer(SIZE_T)::7, H5S_UNLIMITED_F]
      call vh5_error(vh5_create_double_array_nd(posgroupid, "oszicar", 2, dims, maxdims),"vhdf5.F",846)
      call vh5_error(vh5_group_open_or_create(fileid, trim(subgroup), posgroupid),"vhdf5.F",847)
! Used to store if the calculation converged electronically
      dims = [1, 1]; maxdims = [integer(SIZE_T)::1, H5S_UNLIMITED_F]
      call vh5_error(vh5_create_double_array_nd(posgroupid, "electronic_step_converged", 2, dims, maxdims),"vhdf5.F",850)
      call vh5_error(vh5_group_open_or_create(fileid, trim(subgroup), posgroupid),"vhdf5.F",851)

      allocate(label(7))
      label = [character(len=6) :: "N", "E", "dE", "deps", "ncg", "rms", "rms(c)"]
      call vh5_error(vh5_write(posgroupid, "oszicar_label", label),"vhdf5.F",855)
      call vh5_error(vh5_write(posgroupid, "electronic_step_converged_dtype", "bool"),"vhdf5.F",856)
  end subroutine vh5_write_oszicar_start

!> @brief write the output data of OSZICAR for (1._q,0._q) electronic step
  subroutine vh5_write_oszicar_hist_step(fileid, subgroup, convergence_data)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    real(q), intent(in) :: convergence_data(:)
!
    integer(HID_T) :: posgroupid
    integer(SIZE_T) :: start(2), count(2), extend(2)
    integer, save :: completed_steps = 1

    call vh5_error(vh5_group_open_or_create(fileid, trim(subgroup), posgroupid),"vhdf5.F",869)
    start = [1,completed_steps]
    count = [7,1]
    extend = [0,1]
    completed_steps = completed_steps + 1
    call vh5_error(vh5_write_double_subarray_nd(posgroupid, "oszicar", 2, start, count, convergence_data, extend),"vhdf5.F",874)
    call vh5_error(vh5_group_close_writing(posgroupid),"vhdf5.F",875)
  end subroutine vh5_write_oszicar_hist_step

!> @brief write the output of the OSZICAR at the end of the electronic step
  subroutine vh5_write_oszicar_end_elmin_step(fileid, subgroup, is_elmin_converged)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    logical, intent(in)          :: is_elmin_converged

    integer(HID_T) :: posgroupid
    integer(SIZE_T) :: start(2), count(2), extend(2)
    integer, save :: completed_steps = 1
    integer :: converged(1)

    call vh5_error(vh5_group_open_or_create(fileid, trim(subgroup), posgroupid),"vhdf5.F",889)
    start = [1, completed_steps]
    count = [1, 1]
    extend = [0, 1]

    completed_steps = completed_steps + 1
    if (is_elmin_converged) then
        converged = 0
    else
        converged = 1
    end if
    call vh5_error(vh5_write_integer_subarray_nd(posgroupid, "electronic_step_converged", 2, start, count, converged, extend),"vhdf5.F",900)
    call vh5_error(vh5_group_close_writing(posgroupid),"vhdf5.F",901)
  end subroutine vh5_write_oszicar_end_elmin_step

!> @brief start writing of charges and magnetic, or orbital-moments history to hdf5 file
  subroutine vh5_write_moments_hist_start(fileid, subgroup, shape_, prefix, dyn)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    integer, intent(in)          :: shape_(3)  !< [n_orbitals,n_ions,n_components]
    character(LEN=*), intent(in) :: prefix     !< select "spin" or "orbital" moments
    type(dynamics), intent(in)   :: dyn
!-------------------------------------------------------------
    integer(HID_T) :: groupid
    integer(SIZE_T) :: dims(4), maxdims(4), chunk_size(4)
    character(len=6), allocatable :: component_tags(:)
    character(len=1), allocatable :: orbital_tags(:)
!-------------------------------------------------------------
    if (prefix == "spin") then
       allocate(component_tags(4),orbital_tags(4))
       component_tags = ['charge', 'x     ', 'y     ', 'z     ']
       orbital_tags = ['s', 'p', 'd', 'f']
    else if (prefix == "orbital") then
       allocate(component_tags(3),orbital_tags(3))
       component_tags = ['x     ', 'y     ', 'z     ']
       orbital_tags = ['p', 'd', 'f']
    end if
    call vh5_error(vh5_group_open_or_create(fileid, trim(subgroup), groupid),"vhdf5.F",926)
    call vh5_error(vh5_write(groupid, "magnetism/" // prefix // "_moments/components", component_tags(:shape_(3))),"vhdf5.F",927)
    call vh5_error(vh5_write(groupid, "magnetism/" // prefix // "_moments/orbitals", orbital_tags(:shape_(1))),"vhdf5.F",928)
! Create an extensible array that will contain all the charges and magnetic moments
    dims = [shape_, 1]
    maxdims = [int(shape_, kind=SIZE_T), H5S_UNLIMITED_F]
    chunk_size = [1,1,1,min(vh5_default_chunk_size, max(dyn%nsw / dyn%nblock, 1))]
    call vh5_error(vh5_create_double_array_nd(groupid, "magnetism/" // prefix // "_moments/values", size(dims), dims, maxdims, chunk_size),"vhdf5.F",933)

    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",935)
    deallocate(component_tags,orbital_tags)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",937)
  end subroutine vh5_write_moments_hist_start

!> @brief write the charges and magnetic, or orbital-moments for (1._q,0._q) ionic step to hdf5 file
  subroutine vh5_write_moments_hist_step(locid, subgroup, nstep, moments, prefix)
    integer(HID_T), intent(in)   :: locid
    character(LEN=*), intent(in) :: subgroup
    integer, intent(in)          :: nstep          !< ionic step
    real(q), intent(in)          :: moments(:,:,:) !< values of the moments in shape [n_orbitals,n_ions,n_components]
    character(LEN=*), intent(in) :: prefix         !< select "spin" or "orbital" moments
!-------------------------------------------------------------------------------------
    integer ierr
    integer(HID_T) :: groupid
    integer(SIZE_T) :: start(4), count_(4), extend(4)
!-------------------------------------------------------------------------------------
    call vh5_error(vh5_group_open_or_create(locid, trim(subgroup), groupid),"vhdf5.F",952)
    start = [1,1,1,nstep]
    count_ = [shape(moments),1]
    extend = [0,0,0,1]
    ierr = vh5_write_double_subarray_nd(groupid, "magnetism/" // prefix // "_moments/values", size(start), &
        start, count_, moments, extend); call vh5_error(ierr,"vhdf5.F",957)

    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",959)
    call vh5_error(vh5_file_flush(locid),"vhdf5.F",960)
  end subroutine vh5_write_moments_hist_step

  subroutine vh5_write_pair_correlation_start(locid, group, typeinfo, dyn, pair_correlation)
    use base, only: paco_struct
    integer(HID_T), intent(in) :: locid
    character(len=*), intent(in) :: group
    type(type_info), intent(in) :: typeinfo
    type(dynamics), intent(in) :: dyn
    type(paco_struct), intent(in) :: pair_correlation
!
    integer point, num_points, num_pairs, max_steps
    integer(HID_T) :: groupid
    integer(SIZE_T) :: dims(3), maxdims(3)
    real(q) step, itype, jtype, index_
    real(q), allocatable :: distances(:)
    character(len=5), allocatable :: labels(:)
!
    call vh5_error(vh5_group_open_or_create(locid, trim(group), groupid),"vhdf5.F",978)
!
    num_points = pair_correlation%npaco + 1
    num_pairs = typeinfo%ntyp * (typeinfo%ntyp + 1) / 2 + 1  ! add 1 for the total pair correlation function
    max_steps = dyn%nsw / (dyn%nblock * dyn%kblock)
    dims = [num_points, num_pairs, 0]
    maxdims = [num_points, num_pairs, max_steps]
    call vh5_error(vh5_create_double_array_nd(groupid, "function", size(dims), dims, maxdims),"vhdf5.F",985)
!
    allocate(distances(num_points))
    step = pair_correlation%apaco / pair_correlation%npaco
    distances = [(step * (point - 0.5), point=1, num_points)]
    call vh5_error(vh5_write(groupid, "distances", distances),"vhdf5.F",990)
!
    allocate(labels(num_pairs))
    index_ = 1
    labels(index_) = "total"
    do itype = 1, typeinfo%ntyp
       do jtype = itype, typeinfo%ntyp
          index_ = index_ + 1
          labels(index_) = trim(typeinfo%type(itype)) // "~" // trim(typeinfo%type(jtype))
       end do
    end do
    call vh5_error(vh5_write(groupid, "labels", labels),"vhdf5.F",1001)
!
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",1003)
    call vh5_error(vh5_file_flush(locid),"vhdf5.F",1004)
!
  end subroutine vh5_write_pair_correlation_start

  subroutine vh5_write_pair_correlation_step(locid, group, step, dyn, pair_correlation)
    use base, only: paco_struct
    integer(HID_T), intent(in) :: locid
    character(len=*), intent(in) :: group
    integer, intent(in) :: step
    type(dynamics), intent(in) :: dyn
    type(paco_struct), intent(in) :: pair_correlation
!
    integer shape_(2), sample
    integer(HID_T) :: groupid
    integer(SIZE_T) :: start(3), count_(3), extend(3)
!
    shape_ = shape(pair_correlation%sipaco)
    sample = step / (dyn%nblock * dyn%kblock)
    start = [1, 1, sample]
    count_ = [shape_, 1]
    extend = [0, 0, 1]
    call vh5_error(vh5_group_open_or_create(locid, trim(group), groupid),"vhdf5.F",1025)
    call vh5_error(vh5_write_double_subarray_nd(groupid, "function", 3, start, count_, pair_correlation%sipaco, extend),"vhdf5.F",1026)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",1027)
    call vh5_error(vh5_file_flush(locid),"vhdf5.F",1028)
!
  end subroutine vh5_write_pair_correlation_step


  subroutine vh5_write_bandgap(locid, group, method, bandgap)
    use bandgap_struct
    use tutor, only: vtutor
    integer(HID_T), intent(in) :: locid
    character(len=*), intent(in) :: group
    integer, intent(in) :: method
    type(bandgap_info), intent(in) :: bandgap(:)
!
    character(len=23), parameter :: labels(14) = [ &
        "valence band maximum   ", &
        "conduction band minimum", &
        "direct gap bottom      ", &
        "direct gap top         ", &
        "Fermi energy           ", &
        "kx (VBM)               ", &
        "ky (VBM)               ", &
        "kz (VBM)               ", &
        "kx (CBM)               ", &
        "ky (CBM)               ", &
        "kz (CBM)               ", &
        "kx (direct)            ", &
        "ky (direct)            ", &
        "kz (direct)            "]
    character(len=:), allocatable :: dataset
    real(q) data_(size(labels))
    integer(HID_T) groupid
    integer(SIZE_T) dims(3), maxdims(3), start(3), extend(3)
    integer component
    integer, save :: step(2) = 1
!
    select case (method)
    case (use_fermi_weights)
        dataset = "band/gap_from_weight"
    case (use_local_fermi_energy)
        dataset = "band/gap_from_kpoint"
    case default
        call vtutor%bug("Bandgap method to write to hdf5 not implemented", "vhdf5.F", 1069)
    end select
!
    dims = [size(labels), 1, 1]
    maxdims = [size(labels, kind=SIZE_T), size(bandgap, kind=SIZE_T), H5S_UNLIMITED_F]
    extend = [0, 1, 1]
    call vh5_error(vh5_group_open_or_create(locid, trim(group) // "/electron", groupid),"vhdf5.F",1075)
    if (all(step == 1)) &
      call vh5_error(vh5_write(groupid, "band/labels", labels),"vhdf5.F",1077)
    if (step(method) == 1) &
      call vh5_error(vh5_create_double_array_nd(groupid, dataset, size(dims), dims, maxdims),"vhdf5.F",1079)
    do component = 1, size(bandgap)
        associate(gap => bandgap(component))
            data_ = [ &
                gap%fundamental_valence%eigenvalue, &
                gap%fundamental_conduction%eigenvalue, &
                gap%direct_valence%eigenvalue, &
                gap%direct_conduction%eigenvalue, &
                gap%fermi_energy, &
                gap%fundamental_valence%kvector, &
                gap%fundamental_conduction%kvector, &
                gap%direct_valence%kvector &
            ]
            start = [1, component, step(method)]
            call vh5_error(vh5_write_double_subarray_nd(groupid, dataset, size(dims), start, dims, data_, extend),"vhdf5.F",1093)
        end associate
    end do
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",1096)
    step(method) = step(method) + 1
    call vh5_error(vh5_file_flush(locid),"vhdf5.F",1098)
!
  end subroutine vh5_write_bandgap


!> @brief write optical transitions
  subroutine vh5_write_optical_transitions(fileid, opticaltransitions, group, subgroup)
    integer(HID_T), intent(in)  :: fileid
    real(q),intent(in) :: opticaltransitions(:,:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
!-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
!-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",1115)
    call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",1116)

    call vh5_error(vh5_write(subgroupid, "opticaltransitions", opticaltransitions),"vhdf5.F",1118)

    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",1120)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",1121)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",1122)
  end subroutine vh5_write_optical_transitions

!> @brief write contcar trailer
!>
!> writes the informations contained in the contcar trailer in the "traditional"
!> VASP output files to the hdf5 file. Writes lattice velocities and velocities
!>
!> @param fileid hdf5 file handle
!> @param lattice lattice description (def. in lattice.F)
!> @param typeinfo typeinfo data structure (atom types, number of ions, etc. defined in
!>        poscar.f)
!> @param mdynamics MD data structure (atom positions, velocities, etc. defined in poscar.f)
  subroutine vh5_write_positions_trailer(fileid, lattice, typeinfo, mdynamics,nhchains, group, subgroup)
    use poscar, only: has_lattice_velocities
    integer(HID_T), intent(in) :: fileid
    type(latt), intent(in):: lattice
    type(type_info), intent(in):: typeinfo
    type(dynamics), intent(in):: mdynamics
    type(nh_chains), intent(in)      :: nhchains !< nose hoover chains thermostat data
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_POSITIONS)
!-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    integer:: ni
    real(q), dimension(3,typeinfo%nions) :: velocities_out
    real(q), dimension(3) :: tmp
    character(LEN=40) :: my_group, my_subgroup
!-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;      if (present(group))    my_group = group
    my_subgroup = GRP_POSITIONS; if (present(subgroup)) my_subgroup = subgroup
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",1153)
    call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",1154)

    if (has_lattice_velocities(mdynamics, lattice)) then
       call vh5_error(vh5_write(subgroupid, "lattice_velocities", lattice%avel/mdynamics%potim),"vhdf5.F",1157)
    end if

    if (typeinfo%nionp > typeinfo%nions) then
       call vh5_error(vh5_write(subgroupid, "number_empty_sphere_types", typeinfo%nityp(typeinfo%ntyp+1:typeinfo%ntypp)),"vhdf5.F",1161)
       call vh5_error(vh5_write(subgroupid, "position_empty_spheres", mdynamics%posion(:,typeinfo%nions+1:typeinfo%nionp)),"vhdf5.F",1162)
    endif

    if (mdynamics%ibrion/=44 .and. mdynamics%ibrion/=40) then
       do ni = 1, typeinfo%nions
          tmp = mdynamics%vel(:,ni)/mdynamics%potim
          call dirkar(1,tmp,lattice%a)
          velocities_out(:,ni) = tmp
       end do
    else
       velocities_out = mdynamics%vel
    end if
    call vh5_error(vh5_write(subgroupid, "ion_velocities", velocities_out),"vhdf5.F",1174)
    call vh5_error(vh5_write(subgroupid, "direct_coordinates_velocities", 0),"vhdf5.F",1175)
    if (mdynamics%init==1) then
       call vh5_error(vh5_write(subgroupid, "dyn_init", mdynamics%init),"vhdf5.F",1177)
       call vh5_error(vh5_write(subgroupid, "potim", mdynamics%potim),"vhdf5.F",1178)
       IF (nhchains%LINIT) THEN
         call vh5_error(vh5_write(subgroupid, "nhc_thermostat_x", nhchains%x(1:nhchains%NCHAINSMAX)),"vhdf5.F",1180)
         call vh5_error(vh5_write(subgroupid, "nhc_thermostat_p", nhchains%p(1:nhchains%NCHAINSMAX)),"vhdf5.F",1181)
       ELSE
         call vh5_error(vh5_write(subgroupid, "nose_thermostat", mdynamics%snose),"vhdf5.F",1183)
       ENDIF
       call vh5_error(vh5_write(subgroupid, "predictor_coordinates", mdynamics%posion(:,1:typeinfo%nions)),"vhdf5.F",1185)
       call vh5_error(vh5_write(subgroupid, "predictor_coordinates_2", mdynamics%d2),"vhdf5.F",1186)
       call vh5_error(vh5_write(subgroupid, "predictor_coordinates_3", mdynamics%d3),"vhdf5.F",1187)
    end if

    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",1190)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",1191)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",1192)
  end subroutine vh5_write_positions_trailer

!> @brief write force constants
  subroutine vh5_write_force_constants(fileid, force_constants, group, subgroup)
    integer(HID_T), intent(in)       :: fileid
    real(q) :: force_constants(:,:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
!-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
!-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;           if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",1207)
    call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",1208)

    call vh5_error(vh5_write(subgroupid, "force_constants", force_constants),"vhdf5.F",1210)

    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",1212)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",1213)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",1214)
  end subroutine vh5_write_force_constants

!> @brief write hessian matrix
  subroutine vh5_write_hessian(fileid, hessian, group, subgroup)
    integer(HID_T), intent(in)       :: fileid
    real(q) :: hessian(:,:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
!-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
!-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",1229)
    call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",1230)

    call vh5_error(vh5_write(subgroupid, "hessian", hessian),"vhdf5.F",1232)

    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",1234)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",1235)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",1236)
  end subroutine vh5_write_hessian

!> @brief write born effective charges
  subroutine vh5_write_elastic_modulus(fileid, prefix, elastic_modulus, group, subgroup)
    integer(HID_T), intent(in)  :: fileid
    character(len=*), intent(in) :: prefix
    real(q),intent(in) :: elastic_modulus(:,:,:,:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
!-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
!-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",1252)
    call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",1253)

    call vh5_error(vh5_write(subgroupid, prefix // "_elastic_modulus", elastic_modulus),"vhdf5.F",1255)

    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",1257)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",1258)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",1259)
  end subroutine vh5_write_elastic_modulus

!> @brief write born effective charges
  subroutine vh5_write_born_charges(fileid, born_charges, group, subgroup)
    integer(HID_T), intent(in)  :: fileid
    real(q),intent(in) :: born_charges(:,:,:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
!-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
!-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",1274)
    call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",1275)

    call vh5_error(vh5_write(subgroupid, "born_charges", born_charges),"vhdf5.F",1277)

    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",1279)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",1280)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",1281)
  end subroutine vh5_write_born_charges

!> @brief write internal strain tensor
  subroutine vh5_write_internal_strain(fileid, internal_strain, group, subgroup)
    integer(HID_T), intent(in)  :: fileid
    real(q),intent(in) :: internal_strain(:,:,:,:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
!-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
!-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",1296)
    call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",1297)

    call vh5_error(vh5_write(subgroupid, "internal_strain", internal_strain),"vhdf5.F",1299)

    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",1301)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",1302)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",1303)
  end subroutine vh5_write_internal_strain

!> @brief write piezoelectric tensor
  subroutine vh5_write_piezoelectric_tensor(fileid, prefix, piezoelectric_tensor, group, subgroup)
    integer(HID_T), intent(in)  :: fileid
    character(len=*), intent(in) :: prefix
    real(q),intent(in) :: piezoelectric_tensor(:,:,:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
!-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
!-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",1319)
    call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",1320)

    call vh5_error(vh5_write(subgroupid, prefix // "_piezoelectric_tensor", piezoelectric_tensor),"vhdf5.F",1322)

    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",1324)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",1325)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",1326)
  end subroutine vh5_write_piezoelectric_tensor

!> @brief write dipole moment
  subroutine vh5_write_dipole_moment(fileid, prefix, dipole_moment, group, subgroup)
    integer(HID_T), intent(in)  :: fileid
    character(len=*), intent(in) :: prefix
    real(q),intent(in) :: dipole_moment(:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
!-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
!-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",1342)
    call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",1343)

    call vh5_error(vh5_write(subgroupid, prefix // "_dipole_moment", dipole_moment),"vhdf5.F",1345)

    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",1347)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",1348)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",1349)
  end subroutine vh5_write_dipole_moment

!> @brief write spin resolved dipole moment
  subroutine vh5_write_spin_resolved_dipole_moment(fileid, dipole_moment, group, subgroup)
    integer(HID_T), intent(in)  :: fileid
    real(q),intent(in) :: dipole_moment(:,:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
!-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
!-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",1364)
    call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",1365)

    call vh5_error(vh5_write(subgroupid, "spin_resolved_dipole_moment", dipole_moment),"vhdf5.F",1367)

    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",1369)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",1370)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",1371)
  end subroutine vh5_write_spin_resolved_dipole_moment

!> @brief write dielectric tensor
  subroutine vh5_write_dielectric_static(fileid, dielectric_tag, dielectric_tensor, method, group, subgroup)
    integer(HID_T), intent(in)  :: fileid
    character(LEN=*),intent(in) :: dielectric_tag
    real(q),intent(in) :: dielectric_tensor(:,:)
    character(LEN=*), intent(in), optional :: method
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_DIELECTRIC)
!-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
!-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",1388)
    call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",1389)

    call vh5_error(vh5_write(subgroupid, dielectric_tag // "_dielectric_tensor", dielectric_tensor),"vhdf5.F",1391)
    if (present(method)) &
       call vh5_error(vh5_write(subgroupid, "method_dielectric_tensor", method),"vhdf5.F",1393)

    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",1395)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",1396)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",1397)
  end subroutine vh5_write_dielectric_static

  subroutine vh5_write_dielectric_energies_grid(fileid, energy_step, num_energies, prefix, group, subgroup)
    integer(HID_T), intent(in)    :: fileid !< hdf5 file handle
    real(q), intent(in) :: energy_step
    integer, intent(in) :: num_energies
    character(len=*), intent(in), optional :: prefix
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
    real(q) :: energies(num_energies)
    integer :: ii
    energies = [(energy_step * (ii - 1), ii = 1, num_energies)]
    call vh5_write_dielectric_energies_mesh(fileid, energies, prefix, group, subgroup)
  end subroutine vh5_write_dielectric_energies_grid

  subroutine vh5_write_dielectric_energies_mesh(fileid, energies, prefix, group, subgroup)
    integer(HID_T), intent(in)    :: fileid !< hdf5 file handle
    real(q), intent(in) :: energies(:)
    character(len=*), intent(in), optional :: prefix
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
!-------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
    character(LEN=:), allocatable :: dataset
!-------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",1426)
    call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",1427)

    dataset = "energies_dielectric_function"
    if (present(prefix)) dataset = prefix // "_" // dataset
    call vh5_error(vh5_write(subgroupid, dataset, energies),"vhdf5.F",1431)

    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",1433)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",1434)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",1435)
  end subroutine vh5_write_dielectric_energies_mesh

!> @brief write energy dependent dielectric function
  subroutine vh5_write_dielectric_dynamic(fileid, dielectric_function, prefix, group, subgroup)
    integer(HID_T), intent(in)    :: fileid !< hdf5 file handle
    complex(q), intent(in) :: dielectric_function(:,:,:)
    character(len=*), intent(in) :: prefix
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
!-------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
!-------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",1451)
    call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",1452)

    call vh5_error(vh5_write(subgroupid, prefix // "_dielectric_function", dielectric_function),"vhdf5.F",1454)

    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",1456)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",1457)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",1458)
  end subroutine vh5_write_dielectric_dynamic

!> @brief write energy dependent conductivity
  subroutine vh5_write_conductivity_dynamic(fileid, nedos, edos, conductivity, efermi, group, subgroup)
    integer(HID_T), intent(in)    :: fileid !< hdf5 file handle
    real(q),intent(in) :: efermi
    real(q),intent(in) :: edos(nedos), conductivity(nedos,3,3)
    integer,intent(in) :: nedos
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_CONDUCTIVITY)
!-------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    real(q) :: energies(nedos)
    integer :: i
    character(LEN=40) :: my_group, my_subgroup
!-------------------------------------------------------------
    my_group = GRP_RESULTS;         if (present(group))    my_group = group
    my_subgroup = GRP_CONDUCTIVITY; if (present(subgroup)) my_subgroup = subgroup
    call vh5_error(vh5_group_open_or_create(fileid, my_group, groupid),"vhdf5.F",1477)
    call vh5_error(vh5_group_open_or_create(groupid, my_subgroup, subgroupid),"vhdf5.F",1478)

! write stuff
    call vh5_error(vh5_write(subgroupid, "energies", edos),"vhdf5.F",1481)
    call vh5_error(vh5_write(subgroupid, "conductivity_dynamic", conductivity),"vhdf5.F",1482)
    call vh5_error(vh5_write(subgroupid, "efermi", efermi),"vhdf5.F",1483)

! close/deallocate
    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",1486)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",1487)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",1488)
  end subroutine vh5_write_conductivity_dynamic

!> @brief write contents of waveder file
  subroutine vh5_write_waveder(fileid, cder_between_states, kpoints, group, subgroup)
    integer(HID_T), intent(in)    :: fileid !< hdf5 file handle
    COMPLEX(q),intent(in) :: cder_between_states(:,:,:,:,:)
    real(q),intent(in) :: kpoints(:,:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_CONDUCTIVITY)
!-------------------------------------------------------------
    real(qs) :: test
    integer(HID_T) :: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
!-------------------------------------------------------------
    my_group = GRP_RESULTS;         if (present(group))    my_group = group
    my_subgroup = GRP_CONDUCTIVITY; if (present(subgroup)) my_subgroup = subgroup
    call vh5_error(vh5_group_open_or_create(fileid, my_group, groupid),"vhdf5.F",1505)
    call vh5_error(vh5_group_open_or_create(groupid, my_subgroup, subgroupid),"vhdf5.F",1506)

! write stuff
    call vh5_error(vh5_write(subgroupid, "cder_between_states", cder_between_states),"vhdf5.F",1509)
    call vh5_error(vh5_write(subgroupid, "kpoint_coords", kpoints),"vhdf5.F",1510)

! close/deallocate
    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",1513)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",1514)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",1515)
  end subroutine vh5_write_waveder

!> @brief write charge (single spin or all spins)
!>
!> @param fileid hdf5 file handle
!> @param subgroup subgroup to store the position data in
!> @param grid grid_3d grid description type (from mgrid.F)
!> @param chargetot charge array. shape is (ispin, nodes)
  subroutine vh5_write_charge(fileid, subgroup, grid, chargetot)
    use mgrid_struct_def, only: grid_3d
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    type(grid_3d) :: grid
    complex(q) :: chargetot(:,:)
!-------------------------------------------------------------------------------------
    complex(q), allocatable :: work(:,:)
    integer spin
!-------------------------------------------------------------------------------------
! bring density to real space
    allocate(work, source=chargetot)
    do spin = 1, size(work, 2)
      call FFT3D(work(:,spin), grid, 1)
    end do
!
    call vh5_write_grid_quantity(fileid, subgroup, "charge", grid, work)
!
  end subroutine vh5_write_charge

!> ---------------------------- SUBROUTINE VH5 WRITE POTENTIAL---------------------------
!> @brief write local potential on the plane wave grid
!>  for spin-unpolarized calculations it writes (1._q,0._q) dataset
!>  for spin-polarized calculations (ISPIN=2) it writes spin up and spin down potential
!>  for noncollinear calculations it writes the spinor representation (4 datasets)
!> -_------------------------------------------------------------------------------------
  subroutine vh5_write_potential(fileid, key, gridc, cvtot, group, subgroup)
    use mgrid_struct_def, only: grid_3d
    integer(HID_T), intent(in)             :: fileid    !< id of h5 file
    character(len=*), intent(in) :: key                 !< key used to write the potential to the file
    character(len=*), intent(in), optional :: group     !< group in h5 file, default: GRP_RESULTS
    character(len=*), intent(in), optional :: subgroup  !< subgroup in h5 file, default: GRP_POTENTIAL
    type(grid_3d), intent(in) :: gridc                  !< real space grid
    complex(q), intent(in)    :: cvtot(:,:)             !< potential computed by POTLOK
!-------------------------------------------------------------------------------------
    character(len=:), allocatable :: my_group
!-------------------------------------------------------------------------------------
    if (present(group)) then
      my_group = group
    else
      my_group = GRP_RESULTS
    end if
    if (present(subgroup)) then
      my_group = my_group // "/" // subgroup
    else
      my_group = my_group // "/" // GRP_POTENTIAL
    end if
!
    call vh5_write_grid_quantity(fileid, my_group, key, gridc, cvtot)
!
  end subroutine vh5_write_potential

!*******************************************************************
!> write a quantity on a grid to the hdf5 file
!>
!> this function assumes the following conditions are fulfilled
!> - the quantity is given in real space and is real (imaginary part is discarded
!> - the grid is consistent with other quantities in the same group
!*******************************************************************
  subroutine vh5_write_grid_quantity(fileid, group, key, grid, quantity)
    use prec, only: q
    use mgrid, only: grid_3d
    use string, only: str
    use tutor, only: vtutor
    integer(HID_T), intent(in) ::  fileid   !< id of h5 file
    character(len=*), intent(in) :: group   !< name of the group in the HDF5 file
    character(len=*), intent(in) :: key     !< name of the dataset in the HDF5 file
    type(grid_3d), intent(in) :: grid       !< grid on which the quantity is defined
    complex(q), intent(in) :: quantity(:,:) !< quantity which is written to the HDF5 file
!-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid
    integer(HSIZE_T) :: start(4), count_(4), dimensions(4)
    real(q), allocatable :: work(:)
    integer size_work, status_, spin, plane
!-------------------------------------------------------------------------------------
! check for consistency of grid data
    if (grid%nplwv /= grid%ngx * grid%ngy * grid%ngz) then
       call vtutor%bug('internal ERROR: vh5_write_grid_quantity gridc%nplwv = ' // &
         str(grid%nplwv) // 'is not compatible with NGX = ' // str(grid%ngx) // &
         ', NGY = ' // str(grid%ngy) // ', NGZ = ' // str(grid%ngz), "vhdf5.F", 1603)
    end if
!
! allocate work arrays
    size_work = grid%ngx * grid%ngy
    allocate(work(size_work), stat=status_)
    if (status_ /= 0) then
       write (*,*) "WARNING: insufficient memory to write " // key // " to hdf5 file"
       return
    end if
!
    call vh5_error(vh5_group_open_or_create(fileid, group, groupid),"vhdf5.F",1614)
!
! write grid dimensions
    call vh5_error(vh5_write(groupid, "grid", [grid%ngx, grid%ngy, grid%ngz]),"vhdf5.F",1617)
!
! write potential
    dimensions = [grid%ngx, grid%ngy, grid%ngz, size(quantity, 2)]
    call vh5_error(vh5_create_double_array_nd(groupid, key, 4, dimensions),"vhdf5.F",1621)
    do spin = 1, size(quantity, 2)
! write (1._q,0._q) plane of the grid quantity
       do plane = 1, grid%ngz
          call mrg_grid_rl_plane(grid, work, quantity(:,spin), plane)
          start = [1, 1, plane, spin]
          count_ = [grid%ngx, grid%ngy, 1, 1]
          call vh5_error(vh5_write_double_subarray_nd(groupid, key, 4, start, count_, work),"vhdf5.F",1628)
       end do
    end do
!
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",1632)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",1633)
!
  end subroutine vh5_write_grid_quantity

!*******************************************************************
!> initialize the array in the hdf5 file for writing partial charges
!>
!> this function creates a 6D array where partial charges are written
!> to. Two 1D integer arrays for band and kpoint numbers are also
!> written to the file. If they just contain a 0, this means that
!> there will be an average over all selected bands or kpoints.
!*******************************************************************
  subroutine vh5_init_parcharge_output(fileid, grid, nspins, bands, kpoints, group, subgroup)
    use mgrid, only : grid_3d
    integer(HID_T), intent(in) :: fileid !< hdf5 file handle
    type(grid_3d), intent(in) :: grid      ! grid for potential / charge
    integer, intent(in) :: nspins, bands(:), kpoints(:)  ! nr. of spin channels and list of bands and kpoints
    character(len=*), intent(in), optional :: group !< group name (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_PARCHG)
!-------------------------------------------------------------------------------------
    character(len=20), parameter :: dataset_name = "parchg"
    integer, parameter :: rank = 6
    integer(HID_T) :: groupid, subgroupid
    integer(HSIZE_T) :: dimensions(rank)
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    
!-------------------------------------------------------------------------------------
    my_group    = GRP_RESULTS;  if (present(group))    my_group    = trim(group)
    my_subgroup = SUBGRP_PARCHG; if (present(subgroup)) my_subgroup = trim(subgroup)

! open the groups
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",1664)
    call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",1665)

! write grid data
    call vh5_error(vh5_write(subgroupid, "grid", [grid%ngx, grid%ngy, grid%ngz]),"vhdf5.F",1668)

! write lists of bands and kpoints
    call vh5_error(vh5_write(subgroupid, "bands", bands),"vhdf5.F",1671)
    call vh5_error(vh5_write(subgroupid, "kpoints", kpoints),"vhdf5.F",1672)

! prepare the array dimensions and the array
    dimensions = [grid%ngx, grid%ngy, grid%ngz, nspins, size(bands), size(kpoints)]
    call vh5_error(vh5_create_double_array_nd(subgroupid, dataset_name, rank, dimensions),"vhdf5.F",1676)
  
! close the groups again
    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",1679)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",1680)
 end subroutine vh5_init_parcharge_output

!*******************************************************************
!> write the partial charge density to the hdf5 file
!>
!> This function writes a parchg array to the hdf5 file, eiter for
!> a single spin channel or for two of them. The array is written
!> to the given band and kpoint index.
!*******************************************************************

 subroutine vh5_write_to_parcharge_array(fileid, grid, chgtot, iband, ikpt, group, subgroup)
    use mgrid, only : grid_3d
    use string, only: str
    USE tutor, only: vtutor
    integer(HID_T), intent(in) :: fileid !< hdf5 file handle
    type(grid_3d), intent(in) :: grid      ! grid for potential / charge
    complex(q), intent(in) :: chgtot(:,:) !<The (partial) charge density
    integer, intent(in) :: iband !< counting band index
    integer, intent(in) :: ikpt !< counting kpoint index
    character(len=*), intent(in), optional :: group !< group name (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_PARCHG)
!-------------------------------------------------------------------------------------
    complex(q), allocatable :: work_fft(:,:) !< work array for fft3d
    real(q), allocatable :: work(:) !< work array for real space charge writing
    character(len=20), parameter :: dataset_name = "parchg"
    integer, parameter :: rank = 6
    integer(HSIZE_T) :: start(6), count_(6)
    integer(HID_T) :: groupid, subgroupid, kpts_id, band_id
    integer :: spin, size_work, status_, plane
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
!-------------------------------------------------------------------------------------
    my_group    = GRP_RESULTS;  if (present(group))    my_group    = trim(group)
    my_subgroup = SUBGRP_PARCHG; if (present(subgroup)) my_subgroup = trim(subgroup)

! check for consistency of grid data
    if (grid%nplwv /= grid%ngx * grid%ngy * grid%ngz) then
       call vtutor%bug('internal ERROR: vh5_write_grid_quantity gridc%nplwv = ' // &
       str(grid%nplwv) // 'is not compatible with NGX = ' // str(grid%ngx) // &
       ', NGY = ' // str(grid%ngy) // ', NGZ = ' // str(grid%ngz), "vhdf5.F", 1719)
    end if

    allocate(work_fft, source=chgtot)
    do spin = 1, size(work_fft, 2)
       call FFT3D(work_fft(:,spin), grid, 1)
    end do

!
! allocate work arrays
    size_work = grid%ngx * grid%ngy
    allocate(work(size_work), stat=status_)
    if (status_ /= 0) then
       write (*,*) "WARNING: insufficient memory to write " // dataset_name // " to hdf5 file"
       return
    end if

! open the groups
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",1737)
    call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",1738)
    
! write the partial charge density
    do spin = 1, size(work_fft, 2)
! write (1._q,0._q) plane of the grid quantity
       do plane = 1, grid%ngz
          call mrg_grid_rl_plane(grid, work, work_fft(:,spin), plane)
          start = [1, 1, plane, spin, iband, ikpt]
          count_ = [grid%ngx, grid%ngy, 1, 1, 1, 1]
          call vh5_error(vh5_write_double_subarray_nd(subgroupid, dataset_name, rank, start, count_, work),"vhdf5.F",1747)
       end do
    end do

! close the groups again
    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",1752)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",1753)
 end subroutine vh5_write_to_parcharge_array

  subroutine vh5_write_average_potential(fileid, latt_cur, direction, avg_potential, vacuum_potential, group, subgroup)
    use prec
    use LATTICE
    integer(HID_T), intent(in) :: fileid
    integer :: direction
    type(latt)      :: latt_cur
    real(q), intent(in) ::  avg_potential(:), vacuum_potential(2)
    character(len=*), intent(in), optional :: group     !< group in h5 file, default: GRP_RESULTS
    character(len=*), intent(in), optional :: subgroup  !< subgroup in h5 file, default: GRP_POTENTIAL
    integer(HID_T)       :: groupid, subgroupid
    character(len=20)    :: my_group, my_subgroup
    real(q), allocatable :: distance(:)
    integer :: i

    my_group = GRP_RESULTS;      if (present(group))    my_group = trim(group)
    my_subgroup = GRP_POTENTIAL; if (present(subgroup)) my_subgroup = trim(subgroup)

! Generate grid information for the averaged potential
    allocate(distance(size(avg_potential,1)))
    do i = 1, size(avg_potential,1)
      distance(i) = real(i,q) / size(avg_potential,1,q) * latt_cur%ANORM(direction)
    end do

    call vh5_error(vh5_group_open_or_create(fileid, my_group, groupid),"vhdf5.F",1779)
    call vh5_error(vh5_group_open_or_create(groupid, my_subgroup, subgroupid),"vhdf5.F",1780)
    call vh5_error(vh5_write(subgroupid, "vacuum_potential", vacuum_potential),"vhdf5.F",1781)
    call vh5_error(vh5_write(subgroupid, "average_potential_along_IDIPOL", avg_potential),"vhdf5.F",1782)
    call vh5_error(vh5_write(subgroupid, "distance_along_IDIPOL", distance),"vhdf5.F",1783)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",1784)
  end subroutine vh5_write_average_potential

!*******************************************************************
!> write the PAW occupancies to a hdf5 file
!*******************************************************************
  subroutine vh5_write_paw_occupancies(fileid, subgroup, p, t_info, loverl, rholm_store, comm, ispin)
    use pseudo
    use poscar
    use wave
    use constant
    use paw

    integer(HID_T), intent(in)       :: fileid
    character(LEN=*), intent(in)    :: subgroup
    type (type_info)   t_info
    type (potcar),target::  p(t_info%ntyp)
    logical  loverl          ! overlap matrix used ?
    real(q) rholm_store(:)   ! storage for the channel occupancies
    type(communic) :: comm
    integer :: ispin

! local variables
    integer(HID_T) :: unitgroupid
    type (potcar),pointer :: pp
    integer nt, ni, i
    integer ibase, iadd, lymax, lmmax
    integer, allocatable::  nelements(:)
    integer, external :: maxl_aug
!    integer node_me, ionode
    character(LEN=40) :: datasetname

    real(q), allocatable::  buffer(:,:)

!     node_me=0
!     ionode=0
! #ifdef 1
!     node_me=comm%node_me
!     ionode =comm%ionode
! #endif

    allocate(nelements(t_info%nions))

    call vh5_error(vh5_group_open_or_create(fileid, trim(subgroup), unitgroupid),"vhdf5.F",1827)
!=======================================================================
! quick return if possible
!=======================================================================
    if (.not.loverl .or. mimic_us ) return

    lymax =maxl_aug(t_info%ntyp,p)
    lmmax=(lymax+1)**2

    allocate( buffer(lmmax*lmmax, t_info%nions))

!=======================================================================
! cycle all ions and write the required elements
!=======================================================================
    ibase=1

    ion: do ni=1,t_info%nions
       buffer(:,ni)=0
       iadd  =0

       nt=t_info%ityp(ni)
       pp=> p(nt)

       nelements(ni) = 0
       if (do_local(ni)) then
          call retrieve_rholm( buffer(:,ni), rholm_store(ibase:), &
               metric(ibase:), iadd, pp, .true., nelements(ni))

          if (nelements(ni) > lmmax*lmmax) then
             call vtutor%bug('internal error: wrt_rho_paw running out of buffer', "vhdf5.F", 1856)
          endif
          ibase=ibase+iadd
       endif

       CALL m_sum_i(comm, nelements(ni), 1)
       CALL m_sum_d(comm, buffer(:,ni), nelements(ni))

!  IF (NODE_ME==IONODE) THEN
!  write(iu,'("augmentation occupancies",2i4)') ni, nelements
!  write(iu,'(5e15.7)') (buffer(i),i=1,nelements)
!  ENDIF

    enddo ion
    write(datasetname,'(A,I0.2)') "aug_occupancies_ions_s",ispin
    call vh5_error(vh5_write(unitgroupid, datasetname, nelements),"vhdf5.F",1871)
    write(datasetname,'(A,I0.2)') "aug_occupancies_s",ispin
    call vh5_error(vh5_write(unitgroupid, datasetname, buffer),"vhdf5.F",1873)

# 1911

    deallocate(nelements)
    deallocate(buffer)
    call vh5_error(vh5_group_close_writing(unitgroupid),"vhdf5.F",1914)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",1915)
  end subroutine vh5_write_paw_occupancies

!> write PAW strength (Dij) and overlap (Qij) to the hdf5 file
!> note that qij and dij are distributed over NCORE so we must gather and write
  subroutine vh5_write_paw_dij_qij(fileid,wdes,t_info,p,lmdim,cdij,cqij,group)
    use wave, only: ni_global
    use pseudo, only: potcar
    integer(HID_T), intent(in)   :: fileid
    type (wavedes), intent(in)   :: wdes
    type (type_info), intent(in) :: t_info
    type (potcar), intent(in) :: p(:)
    integer, intent(in) :: lmdim
    COMPLEX(q),intent(in) :: cdij(:,:,:,:)
    COMPLEX(q),intent(in) :: cqij(:,:,:,:)
    character(len=*),optional,intent(in) :: group
!local variables
    integer(HID_T) :: groupid
    character(:),allocatable :: my_group
    integer :: i, j, isp, at, ni_local, ni, lmax
    complex(q),allocatable :: cdij_tot(:,:,:,:)
    complex(q),allocatable :: cqij_tot(:,:,:,:)
    integer,allocatable :: lps(:,:)

! collect the data among the mpi ranks
    allocate(cdij_tot(lmdim, lmdim, t_info%nions, wdes%ncdij))
    allocate(cqij_tot(lmdim, lmdim, t_info%nions, wdes%ncdij))
    cdij_tot = cmplx(0.0_q,0.0_q,q)
    cqij_tot = cmplx(0.0_q,0.0_q,q)
    do isp = 1, wdes%ncdij
       do ni_local = 1, wdes%nions
          ni = ni_global(ni_local, wdes%comm_inb)
          cqij_tot(:, :, ni, isp) = cmplx(cqij(:, :, ni_local, isp),kind=q)
          cdij_tot(:, :, ni, isp) = cmplx(cdij(:, :, ni_local, isp),kind=q)
       enddo
    enddo
    CALL M_sum_z(wdes%comm_inb, cdij_tot, size(cdij_tot))
    CALL M_sum_z(wdes%comm_inb, cqij_tot, size(cqij_tot))

    allocate(lps(maxval(p(:)%lmax),size(p)))
    lps(:,:)=-1
    do i=1,size(p)
        lmax = p(i)%lmax
        lps(:lmax,i) = p(i)%lps(:lmax)
    enddo

! write the data to hdf5 file
    my_group = GRP_RESULTS//'/'//GRP_PAW
    if (present(group)) my_group = trim(my_group)
    call vh5_error(vh5_group_open_or_create(fileid, my_group, groupid),"vhdf5.F",1964)

! write data
    call vh5_error(vh5_write(groupid,'dij',cdij_tot),"vhdf5.F",1967)
    call vh5_error(vh5_write(groupid,'qij',cqij_tot),"vhdf5.F",1968)
    call vh5_error(vh5_write(groupid,'lmdim',lmdim),"vhdf5.F",1969)
    call vh5_error(vh5_write(groupid,'ldim',p(:)%ldim),"vhdf5.F",1970)
    call vh5_error(vh5_write(groupid,'lmax',p(:)%lmax),"vhdf5.F",1971)
    call vh5_error(vh5_write(groupid,'lreal',p(:)%lreal),"vhdf5.F",1972)
    call vh5_error(vh5_write(groupid,'lps',lps),"vhdf5.F",1973)

    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",1975)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",1976)
  end subroutine vh5_write_paw_dij_qij

  subroutine vh5_write_atommom(fileid, subgroup, t_info)
    integer(HID_T), intent(in)       :: fileid
    character(LEN=*), intent(in)    :: subgroup
    type (type_info)   t_info

    integer(HID_T) :: unitgroupid

    call vh5_error(vh5_group_open_or_create(fileid, trim(subgroup), unitgroupid),"vhdf5.F",1986)
    call vh5_error(vh5_write(unitgroupid, "atomic_moments", t_info%atomom),"vhdf5.F",1987)
    call vh5_error(vh5_group_close_writing(unitgroupid),"vhdf5.F",1988)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",1989)

  end subroutine vh5_write_atommom


!> @brief write wavefunctions
!>
!> @param fileid hdf5 file handle
!> @param io io handle
!> @param wave_des wave function descriptor (see wave.f)
!> @param wave_spin wave functions (see wave.f)
!> @param latt_ini lattice definition
!> @param efermi fermi energy
!> @param nbands_dump
  subroutine vh5_write_wavefunctions(fileid, io, wave_des, wave_spin, latt_ini, efermi, nbands_dump)
    USE base
    USE string
    USE wave
    USE lattice
    USE main_mpi
    USE tutor, only: vtutor

    integer(HID_T), intent(in)       :: fileid
    type(in_struct) :: io
    type(wavedes)   :: wave_des
    type(wavespin)  :: wave_spin
    type(latt)      :: latt_ini
    real (q)        :: efermi
    integer, optional :: nbands_dump
!----------------------------------------------------------------------
! local variables
    integer(HID_T) :: groupid, spin_groupid, k_groupid
    integer :: ierr, i, j, k, npl, irec, ireclw_old, isp
    integer(HSIZE_T) :: dimensions(2), start(2), count(2)
!    integer :: node_me, ionode
    integer npl_tot, iu0
    integer nb_tot_old
    real(q) :: rdum, rispin, rnb_tot, rnkpts, rnpl
    logical, pointer :: lconjg(:)=>null()
    integer, pointer :: index(:)=>null()
    character(LEN=256) :: spingroupname, kgroupname
!----------------------------------------------------------------------
! local work arrays
    type(wavedes1) wdes1
! FIXME CPTWFP is replaced by CPTFWP via macro in symbol.inc. still need to find out
! why, and what the exact effects are
    complex(q), allocatable :: CPTWFP(:)
!----------------------------------------------------------------------
!
!     node_me=0
!     ionode=0
! #ifdef 1
!     node_me=wave_des%comm%node_me
!     ionode =wave_des%comm%ionode
! #endif

! backup old number of bands
    nb_tot_old = wave_des%nb_tot
    if ( present( nbands_dump ) ) then
       if( nbands_dump > 1 ) then
          wave_des%nb_tot = nbands_dump
       else
          call vtutor%bug(" silly number of bands passed "//str(nbands_dump),"vhdf5.F", 2051 )
       endif
    endif
!
! report to  stdout
!

    if (io%iu6>=0) write(io%iu6,*)'writing wavefunctions to vaspwave.h5'
    if (io%iu0>=0) write(io%iu0,*)'writing wavefunctions to vaspwave.h5'
!
! set up constants
!
    npl_tot = maxval(wave_des%nplwkp_tot)
# 2067

    io%ireclw=max((npl_tot+1)/2,7)*io%icmplx
    ireclw_old=max(max((npl_tot+1)/2,6),((wave_des%nb_tot*3+1)/2+2))*io%icmplx

!
! allocate work arrays
!
    allocate(CPTWFP(npl_tot))
    rispin=wave_des%ispin
    rdum  =io%ireclw
!
! write constants
!
! IF (NODE_ME==IONODE) THEN
    call vh5_error(vh5_group_open_or_create(fileid, GRP_WAVE, groupid),"vhdf5.F",2081)

    call vh5_error(vh5_write(groupid, "rdum", rdum),"vhdf5.F",2083)
    call vh5_error(vh5_write(groupid, "rispin", rispin),"vhdf5.F",2084)
! ENDIF

    irec=2
! IF (NODE_ME==IONODE) THEN
! in order to increase exchangeability of WAVECAR files across IEEE platforms
! avoid INTEGERS on output, write REAL(q) items instead (same below with RNPL)
    rnkpts =wave_des%nkpts
    rnb_tot=wave_des%nb_tot
!
! write some more parameters
!
    call vh5_error(vh5_write(groupid, "rnkpts", rnkpts),"vhdf5.F",2096)
    call vh5_error(vh5_write(groupid, "rnb_tot", rnb_tot),"vhdf5.F",2097)
    call vh5_error(vh5_write(groupid, "efermi", efermi),"vhdf5.F",2098)
    call vh5_error(vh5_write(groupid, "enmax", wave_des%enmax),"vhdf5.F",2099)
    call vh5_error(vh5_write(groupid, "amat", latt_ini%a),"vhdf5.F",2100)
! ENDIF
!
! loop over spin, kpoints, bands
!
    spin: DO ISP=1,wave_des%ISPIN
       write(spingroupname, '(2A)') "spin_"//trim(str(isp))
       call vh5_error(vh5_group_open_or_create(groupid, trim(spingroupname), spin_groupid),"vhdf5.F",2107)
       kpoints: DO K=1,wave_des%NKPTS
          write(kgroupname, '(2A)') "kpoint_"//trim(str(k))
          call vh5_error(vh5_group_open_or_create(spin_groupid, trim(kgroupname), k_groupid),"vhdf5.F",2110)

          IF (wave_des%COMM_KINTER%NCPU > 1) THEN
             IF (wave_des%COMM_KINTER%NODE_ME /= 1) THEN
                IF (MOD(K-1,wave_des%COMM_KINTER%NCPU).NE.wave_des%COMM_KINTER%NODE_ME-1) THEN
                   CYCLE
                ELSE
                   CALL M_send_z(wave_des%COMM_KINTER, 1, wave_spin%CPTWFP(1,1,K,ISP), SIZE(wave_spin%CPTWFP,1)*SIZE(wave_spin%CPTWFP,2))
                END IF
             ELSE
                IF (MOD(K-1,wave_des%COMM_KINTER%NCPU).NE.wave_des%COMM_KINTER%NODE_ME-1) THEN
                   CALL M_recv_z(wave_des%COMM_KINTER, MOD(K-1,wave_des%COMM_KINTER%NCPU)+1, wave_spin%CPTWFP(1,1,K,ISP), SIZE(wave_spin%CPTWFP,1)*SIZE(wave_spin%CPTWFP,2))
                ENDIF
             ENDIF
          END IF

          CALL SETWDES(wave_des,WDES1,K)
# 2129

          npl=wave_des%nplwkp_tot(k)
          rnpl=npl
! write number of plane waves, k-point coordinates and all eigenvalues and
! occupation numbers for current k-point K
! write eigenvalues in real format
! IF (NODE_ME==IONODE) THEN
          call vh5_error(vh5_write(k_groupid, "num_planewaves", npl),"vhdf5.F",2136)
          call vh5_error(vh5_write(k_groupid, "vkpt", wave_des%vkpt(:,k)),"vhdf5.F",2137)
          call vh5_error(vh5_write(k_groupid, "celtot", wave_spin%celtot(1:wave_des%nb_tot,k,isp)),"vhdf5.F",2138)
          call vh5_error(vh5_write(k_groupid, "fertot", wave_spin%fertot(1:wave_des%nb_tot,k,isp)),"vhdf5.F",2139)
! ENDIF

          dimensions = [npl, wave_des%nb_tot]
# 2145

          call vh5_error(vh5_create_complex_array_nd(k_groupid, "wave", 2, dimensions),"vhdf5.F",2146)


          DO J=1,wave_des%nb_tot
# 2152

             CALL MRG_PW_BAND(WDES1, J, CPTWFP, wave_spin%CPTWFP(1,1,K,ISP))

             irec=irec+1

             start = [1,J]
             count = [npl,1]
             call vh5_error(vh5_write_double_complex_subarray_nd(k_groupid, "wave", 2, start, count, CPTWFP),"vhdf5.F",2159)
          enddo

# 2164

          call vh5_error(vh5_group_close_writing(k_groupid),"vhdf5.F",2165)
       enddo kpoints
       call vh5_error(vh5_group_close_writing(spin_groupid),"vhdf5.F",2167)
    enddo spin

    deallocate(CPTWFP)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",2171)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",2172)
  end subroutine vh5_write_wavefunctions

!> read force constants from hdf5 file
  subroutine vh5_read_force_constants(fileid,force_constants)
    use tutor, only: vtutor
    use string, only: str
    integer(HID_T),intent(in) :: fileid
    real(q),allocatable,intent(inout) :: force_constants(:,:)
! local variables
    integer :: dims(2)
    real(q) :: fc(1,1)
    integer(HID_T) :: ih5inresultsgroup_id
    integer :: ierr

    ierr = get_dimensions(fileid, '/results/linear_response/force_constants', 2, dims)
    if (ierr/=0) then
        call vtutor%error("Could not find '/results/linear_response/force_constants' in " // VASPIN)
    endif
! check the dimensions
    if (any(shape(force_constants)/=dims)) then
        call vtutor%error("The dimensions of the force_constants array in the HDF5 file "//str(dims(1))//&
                          " are incompatible with the current VASP run "//str(size(force_constants,1))//&
                          ". Please check if the POSCAR is the same you used for the supercell calculation")
    endif
    call vh5_error(vh5_read(fileid, '/results/linear_response/force_constants', force_constants),"vhdf5.F",2197)
  end subroutine vh5_read_force_constants

!> @brief read wave function header
!>
!> @param fileid hdf5 file handle
!> @param wave_des wave function descriptor (see wave.f)
!> @param wave_spin wave functions (see wave.f)
!> @param latt_ini lattice definition
!> @param latt_cur lattice definition
!> @param enmaxi
!> @param istart
!> @param iu0
  function vh5_read_wavefunction_header(fileid, wave_des, latt_ini, latt_cur, enmaxi, istart, iu0 ) result(ierr)
    USE base
    USE wave
    USE lattice
    USE main_mpi

    integer(HID_T), intent(in)       :: fileid
    type(wavedes)   :: wave_des
    type(wavespin)  :: wave_spin
    type(latt)      :: latt_ini
    type(latt)      :: latt_cur
    real (q)        :: enmaxi
    integer         :: istart
    integer         :: iu0
!-------------------------------------------------------------
! local variables
    integer(HID_T) :: groupid
    integer :: ierr, i, j, k, n, ionode, irec, ireclw_old, isp
    integer :: node_me, npl
    integer :: npl_tot, nbandf, nkptsf
    real(q) :: rdum, rispin, rnb_tot, rnkpts, rnpl!, enmax
    logical :: ldiff
!-------------------------------------------------------------

    node_me=0
    ionode=0

    node_me=wave_des%comm%node_me
    ionode =wave_des%comm%ionode


! open group
    ierr = vh5_group_open(fileid, GRP_WAVE, groupid)
    if(ierr .EQ. 0) then
       if (iu0 >= 0) then
          IF (NODE_ME==IONODE) write(iu0,*) 'found vaspwave.h5, reading the header'
       endif

       call vh5_error(vh5_read(groupid, "rnkpts", rnkpts),"vhdf5.F",2248)
       call vh5_error(vh5_read(groupid, "rnb_tot", rnb_tot),"vhdf5.F",2249)
       call vh5_error(vh5_read(groupid, "enmax", enmaxi),"vhdf5.F",2250)
       call vh5_error(vh5_read(groupid, "amat", latt_ini%a),"vhdf5.F",2251)


       NKPTSF=NINT(rnkpts)
       NBANDF=NINT(rnb_tot)

       CALL LATTIC(LATT_INI)

       IF (ISTART==2 .AND. ENMAXI /= wave_des%ENMAX) THEN
          CALL vtutor%error('ERROR: ENMAX changed please set ISTART to 1')
       ENDIF

       IF (.NOT. (NBANDF==wave_des%NB_TOT .OR. (wave_des%NRSPINORS==2 .AND.  NBANDF*wave_des%NRSPINORS==wave_des%NB_TOT))) THEN
          IF (IU0 >= 0) WRITE(IU0,'(2X,A,I6,A,I6)') &
               'number of bands has changed, file:',NBANDF,' present:',wave_des%NB_TOT
          IF (IU0 >= 0) WRITE(IU0,'(2X,A,I6,A,I6)') &
               'trying to continue reading WAVECAR, but it might fail'
       ENDIF
       IF (NKPTSF/=wave_des%NKPTS) THEN
          IF (IU0 >= 0) WRITE(IU0,'(2X,A,I6,A,I6)') &
               'number of k-points has changed, file:',NKPTSF,' present:',wave_des%NKPTS
          IF (IU0 >= 0) WRITE(IU0,'(2X,A,I6,A,I6)') &
               'trying to continue reading WAVECAR, but it might fail'
       ENDIF
       IF (ISTART ==1) THEN

          LDIFF=.FALSE.
          DO I=1,3
             DO J=1,3
                IF (ABS(LATT_INI%A(I,J)-LATT_CUR%A(I,J)) > 1E-4) LDIFF=.TRUE.
             ENDDO
          ENDDO
          IF (ENMAXI /= wave_des%ENMAX) LDIFF=.TRUE.
          IF (LDIFF) THEN
             IF (IU0>=0) &
                  WRITE(IU0,*)'WAVECAR: different cutoff or change in lattice found'
          ENDIF

       ENDIF

       RETURN
    else
       latt_ini%a = latt_cur%a
       call lattic(latt_ini)
       if (istart==3) then
          if (iu0 >= 0) then
             write(iu0,*) "ERROR: can't restart (ISTART=3) with wavefunctions on file"
          endif
       else
          istart = 0
       endif
    endif

  end function vh5_read_wavefunction_header


!> @brief read wave functions
!>
!> @param fileid hdf5 file handle
!> @param io io handle
!> @param wave_des wave function descriptor (see wave.f)
!> @param wave_spin wave functions (see wave.f)
!> @param latt_cur lattice definition
!> @param latt_ini lattice definition
!> @param istart
!> @param efermi fermi energy
  function vh5_read_wavefunctions(fileid, io, wave_des, wave_spin, grid, latt_cur, latt_ini, istart, efermi ) result(ierr)
    USE base
    USE string
    USE wave_high
    use mgrid
    USE lattice
    USE main_mpi

    integer(HID_T), intent(in)       :: fileid
    type(in_struct) :: io
    type(wavedes)   :: wave_des
    type(wavespin)  :: wave_spin
    type(grid_3d)   :: grid
    type(latt)      :: latt_ini
    type(latt)      :: latt_cur
    real (q)        :: efermi
    integer         :: istart
!-------------------------------------------------------------
! local variables
    integer(HID_T) :: groupid, spin_groupid, k_groupid
    integer :: ierr, i, j, k, n, ionode, irec, ireclw_old, isp, kpl
    integer :: node_me, npl, nmaxf, ifail, irecl, is, ispinread, iu0
    integer :: npl_tot, nbandf, nkptsf, malloc, ncomm, nplread
    real(q) :: rdum, rispin, rnb_tot, rnkpts, rnpl, enmax, rbandf, rkptsf, rtag, enmaxf
    logical :: ldiff
    logical, pointer :: lconjg(:)=>null()
    integer, pointer :: index(:)=>null()
    character(LEN=256) :: spingroupname, kgroupname
!----------------------------------------------------------------------
! local work arrays
    type(wavedes1) wave_des1
! FIXME CPTWFP is replaced by CPTFWP via macro in symbol.inc. still need to find out
! why, and what the exact effects are
    complex(q), allocatable :: CPTWFP(:)
    REAL(q) VKPT(3)
    INTEGER,ALLOCATABLE    :: IND(:),INDI(:)
    COMPLEX(q),ALLOCATABLE :: CW1(:),CW2(:)
    CHARACTER(:),ALLOCATABLE :: MSG
    LOGICAL ::  SINGLE_PREC, SPAN_RECORDS
    REAL(q), ALLOCATABLE :: INBUF(:)
    integer(HSIZE_T) :: start(2), count(2)
!----------------------------------------------------------------------

    iu0 = io%iu0

    node_me=0
    ionode=0

    node_me=wave_des%comm%node_me
    ionode =wave_des%comm%ionode


    ierr = vh5_group_open(fileid, GRP_WAVE, groupid)
    if(ierr .EQ. 0) then
       if (iu0 >= 0) then
          IF (NODE_ME==IONODE) write(iu0,*) 'found vaspwave.h5, reading file, group = ', GRP_WAVE
       endif

       RTAG=0
       call vh5_error(vh5_read_double_scalar(groupid, "rispin", rispin),"vhdf5.F",2376)

       ISPINREAD=NINT(RISPIN)
       call vh5_error(vh5_read(groupid, "rnkpts", rkptsf),"vhdf5.F",2379)
       call vh5_error(vh5_read(groupid, "rnb_tot", rbandf),"vhdf5.F",2380)
       call vh5_error(vh5_read(groupid, "enmax", enmaxf),"vhdf5.F",2381)
       call vh5_error(vh5_read(groupid, "amat", latt_ini%a),"vhdf5.F",2382)
       call vh5_error(vh5_read(groupid, "efermi", efermi),"vhdf5.F",2383)

       IREC=2
       NKPTSF=NINT(RKPTSF)
       NBANDF=NINT(RBANDF)

       CALL LATTIC(LATT_INI)
!=======================================================================
! read WAVECAR file, number of bands agree
!=======================================================================
       IF (.NOT.(wave_des%NRSPINORS==2 .AND. NBANDF*wave_des%NRSPINORS == wave_des%NB_TOT)) THEN

          spin:    DO ISP=1, MIN(wave_des%ISPIN, ISPINREAD)
             write(spingroupname, '(2A)') "spin_"//trim(str(isp))
             call vh5_error(vh5_group_open(groupid, trim(spingroupname), spin_groupid),"vhdf5.F",2397)
             kpoints: DO K=1,wave_des%NKPTS
                IF ( K> NKPTSF ) THEN
                   IREC=IRECL
                   write(kgroupname, '(2A)') "kpoint_"//trim(str(NKPTSF))
                ELSE
                   IRECL=IREC
                   write(kgroupname, '(2A)') "kpoint_"//trim(str(k))
                ENDIF
                call vh5_error(vh5_group_open(spin_groupid, trim(kgroupname), k_groupid),"vhdf5.F",2406)

                CALL SETWDES(wave_des,wave_des1,K)
# 2411


                IF (NODE_ME==IONODE) THEN
                ifail = 0
                call vh5_error(vh5_read(k_groupid, "num_planewaves", nplread),"vhdf5.F",2415)
                call vh5_error(vh5_read(k_groupid, "vkpt", vkpt),"vhdf5.F",2416)
                call vh5_error(vh5_read(k_groupid, "celtot", wave_spin%CELTOT(:,K,ISP)),"vhdf5.F",2417)
                call vh5_error(vh5_read(k_groupid, "fertot", wave_spin%FERTOT(:,K,ISP)),"vhdf5.F",2418)
                ENDIF

                CALL M_bcast_i( wave_des%COMM, IFAIL, 1)
                IF (IFAIL /=0) GOTO 230

                CALL M_bcast_i( wave_des%COMM, IREC, 1)
                CALL M_bcast_i( wave_des%COMM, NPLREAD, 1)
                NPL=wave_des%NPLWKP_TOT(K)

                MALLOC=MAX(NPL, NPLREAD)
                ALLOCATE(CW1(MALLOC),CW2(MALLOC),IND(MALLOC),INDI(MALLOC))

                IF (NODE_ME==IONODE) THEN
! create an index to allow for change of cutoff or cell size
                IF (ISTART==2) THEN
! gK: fix ISTART=2 in INWAV_FAST
                   CALL REPAD_INDEX_ARRAY(GRID, wave_des%VKPT(:,K), VKPT, LATT_INI%B,  LATT_INI%B, &
                        wave_des%ENMAX, ENMAXF, NPL/wave_des%NRSPINORS, NPLREAD/wave_des%NRSPINORS, IND, INDI, MALLOC, IFAIL )
                ELSE
                   CALL REPAD_INDEX_ARRAY(GRID, wave_des%VKPT(:,K), VKPT, LATT_CUR%B,  LATT_INI%B, &
                        wave_des%ENMAX, ENMAXF, NPL/wave_des%NRSPINORS, NPLREAD/wave_des%NRSPINORS, IND, INDI, MALLOC, IFAIL )
                ENDIF
                ENDIF

                CALL M_bcast_i( wave_des%COMM, IFAIL, 1)
                IF (IFAIL /=0) GOTO 220

                band: DO J=1,NBANDF
                   IREC=IREC+1
                   IF (J>wave_des%NB_TOT) CYCLE
                   start = [1,J]
                   count = [nplread,1]

                   IF (NODE_ME==IONODE) THEN
                   call vh5_error(vh5_read_double_complex_subarray_nd(k_groupid, "wave", 2, start, count, cw2),"vhdf5.F",2453)

                   CW1=0
                   DO IS=1,wave_des%NRSPINORS
! store the wave function coefficients according to new cutoff
                      CALL REPAD_WITH_INDEX_ARRAY( MALLOC, IND, INDI, &
                           CW1((IS-1)*NPL/wave_des%NRSPINORS+1), CW2((IS-1)*NPLREAD/wave_des%NRSPINORS+1))
                   ENDDO

                   ENDIF

                   IF (wave_des%COMM_KINTER%NCPU.GT.1) THEN
                      CALL M_bcast_z(wave_des%COMM_KINTER,CW1,SIZE(CW1))
                   END IF

# 2471

                   CALL DIS_PW_BAND(wave_des1, J, CW1, wave_spin%CPTWFP(1,1,K,ISP))

                ENDDO band

# 2478

                DEALLOCATE(CW1,CW2,IND,INDI)

                call vh5_error(vh5_group_close(k_groupid),"vhdf5.F",2481)
             ENDDO kpoints


             IF (NKPTSF > wave_des%NKPTS) THEN
                IREC=IREC+(NKPTSF-wave_des%NKPTS)*(wave_des%NB_TOT+1)
             ENDIF

! copy eigenvalues and weights to all nodes
             NCOMM=wave_des%NB_TOT*wave_des%NKPTS
             CALL M_bcast_z(wave_des%COMM, wave_spin%CELTOT(1,1,ISP),NCOMM )
             CALL M_bcast_d(wave_des%COMM, wave_spin%FERTOT(1,1,ISP),NCOMM )

             call vh5_error(vh5_group_close(spin_groupid),"vhdf5.F",2494)
          ENDDO spin

          IF (ISPINREAD > wave_des%ISPIN .AND. IU0>=0 ) THEN
             WRITE(IU0,*) 'down-spin wavefunctions not read'
          ENDIF

          IF (NBANDF<wave_des%NB_TOT) THEN
             IF (IU0>=0) WRITE(IU0,*) 'random initialization beyond band ',NBANDF
             CALL WFINIT(wave_des, wave_spin, 1E10_q, NBANDF+1) ! ENINI=1E10 not cutoff restriction
          ENDIF

          IF (IU0 >= 0) WRITE(IU0,*) 'the vaspwave.h5 file was read successfully'

          IF (wave_des%ISPIN<=ISPINREAD) THEN
             RETURN
          ENDIF
!
!  spin down is missing
!
          IF (IU0>=0) &
               WRITE(IU0,*) 'No down-spin wavefunctions found', &
               &             ' --> setting down-spin equal up-spin ...'

          DO K=1,wave_des%NKPTS
             wave_spin%CELTOT(1:wave_des%NB_TOT,K,2)=wave_spin%CELTOT(1:wave_des%NB_TOT,K,1)
             wave_spin%FERTOT(1:wave_des%NB_TOT,K,2)=wave_spin%FERTOT(1:wave_des%NB_TOT,K,1)
             NPL=wave_des%NPLWKP(K)
             wave_spin%CPTWFP(1:NPL,1:wave_des%NBANDS,K,2)=wave_spin%CPTWFP(1:NPL,1:wave_des%NBANDS,K,1)
          ENDDO


          RETURN
!=======================================================================
! read collinear WAVECAR file for a non collinear run
!=======================================================================
       ELSE
          wave_spin%CPTWFP=0

          spin2:    DO ISP=1,ISPINREAD
             IF (IU0>=0.AND. ISP==1) &
                  WRITE(IU0,*) 'reading wavefunctions of collinear run, up'
             IF (IU0>=0.AND. ISP==2) &
                  WRITE(IU0,*) 'reading wavefunctions of collinear run, down'
             write(spingroupname, '(A,I0.2)') "spin_", isp
             call vh5_error(vh5_group_open(groupid, trim(spingroupname), spin_groupid),"vhdf5.F",2539)
             kpoints2: DO K=1,NKPTSF

                write(kgroupname, '(A,I0.4)') "kpoint_", k
                call vh5_error(vh5_group_open(spin_groupid, trim(kgroupname), k_groupid),"vhdf5.F",2543)

                IF ( K> NKPTSF ) THEN
                   IREC=IRECL
                ELSE
                   IRECL=IREC
                ENDIF

                CALL SETWDES(wave_des,wave_des1,K)


                IF (NODE_ME==IONODE) THEN
                ifail = 0
                call vh5_error(vh5_read_integer_scalar(k_groupid, "num_planewaves", nplread),"vhdf5.F",2556)
                call vh5_error(vh5_read_double_array_1d(k_groupid, "vkpt", vkpt),"vhdf5.F",2557)
                call vh5_error(vh5_read_double_complex_array_1d(k_groupid, "celtot", wave_spin%CELTOT(:,K,ISP)),"vhdf5.F",2558)
                call vh5_error(vh5_read_double_array_1d(k_groupid, "fertot", wave_spin%FERTOT(:,K,ISP)),"vhdf5.F",2559)
                ENDIF

                CALL M_bcast_i( wave_des%COMM, IFAIL, 1)
                IF (IFAIL /=0) GOTO 230

                CALL M_bcast_i( wave_des%COMM, IREC, 1)
                CALL M_bcast_i( wave_des%COMM, NPLREAD, 1)
                NPL=wave_des%NPLWKP_TOT(K)/2

                MALLOC=MAX(NPL, NPLREAD)
                ALLOCATE(CW1(2*MALLOC),CW2(2*MALLOC),IND(MALLOC),INDI(MALLOC))

                IF (NODE_ME==IONODE) THEN
                IF (ISTART==2) THEN
! gK: fix ISTART=2 in INWAV_FAST
                   CALL REPAD_INDEX_ARRAY(GRID, wave_des%VKPT(:,K), VKPT, LATT_INI%B,  LATT_INI%B, &
                        wave_des%ENMAX, ENMAXF, NPL, NPLREAD, IND, INDI, MALLOC, IFAIL )
                ELSE
                   CALL REPAD_INDEX_ARRAY(GRID, wave_des%VKPT(:,K), VKPT, LATT_CUR%B,  LATT_INI%B, &
                        wave_des%ENMAX, ENMAXF, NPL, NPLREAD, IND, INDI, MALLOC, IFAIL )
                ENDIF ! istart == 2
                ENDIF

                CALL M_bcast_i( wave_des%COMM, IFAIL, 1)
                IF (IFAIL /=0) GOTO 220

                IF (ISP==1) THEN
                   DO J=1,NBANDF
                      IREC=IREC+1
                      start = [1,J]
                      count = [nplread,1]

                      IF (NODE_ME==IONODE) THEN
                      call vh5_error(vh5_read_double_complex_subarray_nd(k_groupid, "wave", 2, start, count, cw2),"vhdf5.F",2593)

                      CW1=0
                      CALL REPAD_WITH_INDEX_ARRAY( MALLOC, IND, INDI, CW1, CW2)
                      ENDIF


                      IF (wave_des%COMM_KINTER%NCPU.GT.1) THEN
                         CALL M_bcast_z(wave_des%COMM_KINTER,CW1,SIZE(CW1))
                      END IF ! wave_des

                      CALL DIS_PW_BAND(wave_des1, J, CW1, wave_spin%CPTWFP(1,1,K,1))

! copy immediately to second panel
                      wave_spin%CELTOT(J+NBANDF,K,1)=wave_spin%CELTOT(J,K,1)
                      wave_spin%FERTOT(J+NBANDF,K,1)=wave_spin%FERTOT(J,K,1)

                      CW1(NPL+1:2*NPL)=CW1(1:NPL)
                      CW1(1:NPL)=0

                      CALL DIS_PW_BAND(wave_des1, J+NBANDF, CW1, wave_spin%CPTWFP(1,1,K,1))
                   ENDDO
                ELSE ! ISP==1
                   DO J=1,NBANDF
                      IREC=IREC+1
                      start = [1,J]
                      count = [nplread,1]

                      IF (NODE_ME==IONODE) THEN
                      call vh5_error(vh5_read_double_complex_subarray_nd(k_groupid, "wave", 2, start, count, cw2),"vhdf5.F",2622)

                      CW1=0
                      CALL REPAD_WITH_INDEX_ARRAY( MALLOC, IND, INDI, CW1(NPL+1), CW2)
                      ENDIF

                      IF (wave_des%COMM_KINTER%NCPU.GT.1) THEN
                         CALL M_bcast_z(wave_des%COMM_KINTER,CW1,SIZE(CW1))
                      END IF

                      CALL DIS_PW_BAND(wave_des1, J+NBANDF, CW1, wave_spin%CPTWFP(1,1,K,1))
                   ENDDO

                ENDIF  ! ISP == 1

                DEALLOCATE(CW1,CW2,IND,INDI)

                call vh5_error(vh5_group_close(k_groupid),"vhdf5.F",2639)
             ENDDO kpoints2

             IF (NKPTSF > wave_des%NKPTS) THEN
                IREC=IREC+(NKPTSF-wave_des%NKPTS)*(NBANDF+1)
             ENDIF

             call vh5_error(vh5_group_close(spin_groupid),"vhdf5.F",2646)
          ENDDO spin2

! copy eigenvalues and weights to all nodes
          NCOMM=wave_des%NB_TOT*wave_des%NKPTS
          CALL M_bcast_z(wave_des%COMM, wave_spin%CELTOT(1,1,1),NCOMM )
          CALL M_bcast_d(wave_des%COMM, wave_spin%FERTOT(1,1,1),NCOMM )

          IF (IU0 >= 0) WRITE(IU0,*) 'the vaspwave.h5 file was read successfully'
          RETURN
       ENDIF ! num bands agree

       MSG = 'ERROR: while reading vaspwave.h5, file is incompatible'
       IF (wave_des%ENMAX /= ENMAXF) &
          MSG = MSG // ' the energy cutoff has changed (new,old) ' // &
             str(wave_des%ENMAX) // " " // str(ENMAXF)
       IF (NBANDF /= wave_des%NB_TOT) &
          MSG = MSG // ' the number of bands has changed (new,old) ' // &
             str(wave_des%NB_TOT) // " " // str(NBANDF)
       IF (NKPTSF /= wave_des%NKPTS) &
          MSG = MSG // ' the number of k-points has changed (new,old) ' // &
             str(wave_des%NKPTS) // " " // str(NKPTSF)

       CALL vtutor%error(MSG)

    else
! could not open the group WAVE_GRP in file, handle error

       CALL vtutor%error('Error reading hdf5 wave function file: group does not exist')

    endif

! can not do anything with WAVECAR
! hard stop, pull all breaks
!=======================================================================
200 CONTINUE

    CALL vtutor%error('ERROR: while reading vaspwave.h5, header is corrupt')

220 CONTINUE

    CALL vtutor%error('ERROR: while reading vaspwave.h5, plane wave coefficients changed ' &
       // str(NPL) // ' ' // str(NPLREAD))

230 CONTINUE

    CALL vtutor%error('ERROR: while reading eigenvalues from vaspwave.h5 ' &
       // str(K) // ' ' // str(ISP))

240 CONTINUE

    CALL vtutor%error('ERROR: while reading plane wave coef. from vaspwave.h5 ' &
       // str(K) // ' ' // str(ISP) // ' ' // str(J))

  end function vh5_read_wavefunctions

!> @brief write projectors to hdf5 file
  subroutine vh5_write_projectors(fileid, typeinfo, wave_des, wave_spin, lorbit, &
       lpar, par, lchar, phase, group, subgroup)
    integer(HID_T), intent(in)   :: fileid !< hdf5 file handle
    type(type_info), intent(in)  :: typeinfo
    type(wavedes)   :: wave_des
    type(wavespin)  :: wave_spin
    integer :: lorbit, lpar
    real(q), dimension(wave_des%nb_tot, wave_des%nkpts, lpar, typeinfo%nionp, wave_des%ncdij) :: par
    character(len=*), dimension(:) :: lchar
    complex(q), dimension(:,:,:,:,:) :: phase
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_PROJECTORS)
!-------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    character(len=256) :: lcharout
    integer :: nl
    integer(HSIZE_T), dimension(5) :: par_shape;
    integer(HSIZE_T), dimension(5) :: phase_shape;
    character(LEN=40) :: my_group, my_subgroup
!-------------------------------------------------------------
    my_group    = GRP_RESULTS;    if (present(group))    my_group    = trim(group)
    my_subgroup = GRP_PROJECTORS; if (present(subgroup)) my_subgroup = trim(subgroup)
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",2725)
    call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",2726)

! write procar type
    if (lorbit == 11) then
       call vh5_error(vh5_write(subgroupid, "procar_type", 'PROCAR lm decomposed'),"vhdf5.F",2730)
    else if (lorbit == 12) then
       call vh5_error(vh5_write(subgroupid, "procar_type", 'PROCAR lm decomposed + phase'),"vhdf5.F",2732)
    else
       call vh5_error(vh5_write(subgroupid, "procar_type", 'PROCAR new format'),"vhdf5.F",2734)
    end if
!
! write parameters
! FIXME these parameters are also written elsewhere. we need to discuss how to structure
! a file where these things are only written once. ATM, I'm writing them out because
! I'm not sure whether they correspond to the same state of the program as the other
! instances
    call vh5_error(vh5_write(subgroupid, "nb_tot", wave_des%nb_tot),"vhdf5.F",2742)
    call vh5_error(vh5_write(subgroupid, "nionp", typeinfo%nionp),"vhdf5.F",2743)
    call vh5_error(vh5_write(subgroupid, "nkpoints", wave_des%nkpts),"vhdf5.F",2744)
    call vh5_error(vh5_write(subgroupid, "kpoints", wave_des%vkpt),"vhdf5.F",2745)
    call vh5_error(vh5_write(subgroupid, "celtot", wave_spin%celtot),"vhdf5.F",2746)
!
! write lchars
!
    call vh5_error(vh5_write(subgroupid, "lchar", lchar),"vhdf5.F",2750)
!
! write par
!
    par_shape = shape(par)
    call vh5_error(vh5_write_double_array_nd(subgroupid, "par", 5, par_shape, par),"vhdf5.F",2755)
!
! write phase
!
    if (lorbit == 12) then
       phase_shape = shape(phase)
       call vh5_error(vh5_write_double_complex_array_nd(subgroupid, "phase", 5, phase_shape, phase),"vhdf5.F",2761)
    end if
!
! calculate and write parsum and sumion
!

!
! close/deallocate
!
    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",2770)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",2771)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",2772)
  end subroutine vh5_write_projectors

  subroutine vh5_write_locproj(fileid, lprj_functions, lprj_covl, lnoncol, group, subgroup)
    use string, only: str

    integer(HID_T), intent(in)       :: fileid !< hdf5 file handle
    TYPE (LPRJ_function), ALLOCATABLE :: lprj_functions(:)
    COMPLEX(q), ALLOCATABLE :: LPRJ_COVL(:,:,:,:)
    integer(HID_T) :: nproj !nr. of proj.
    integer(HID_T):: ip, groupid, subgroupid, paramgroupid
    logical, intent(in), optional :: lnoncol
    character(len=*), intent(in), optional :: group
    character(len=*), intent(in), optional :: subgroup
    character(len=15), allocatable :: radial_type_list(:)
    character(len=10), allocatable :: ang_type_list(:)
    REAL(q), allocatable :: lprj_coord_list(:,:)           !< position of the localized orbital in reduced coordinates
    CHARACTER(LEN=15) :: PRKIND(3)=(/"PAW projector  ","PS partial wave","Hydrogen-like  "/)
    CHARACTER(LEN=10) :: LABEL(1:7,-5:3)
    DATA LABEL(1:7,-5) &
     &/"  sp3d2-1 ","  sp3d2-2 ","  sp3d2-3 ","  sp3d2-4 ","  sp3d2-5 ","  sp3d2-6 ","          "/
      DATA LABEL(1:7,-4) &
     &/"   sp3d-1 ","   sp3d-2 ","   sp3d-3 ","   sp3d-4 ","   sp3d-5 ","          ","          "/
      DATA LABEL(1:7,-3) &
     &/"   sp3-1  ","   sp3-2  ","   sp3-3  ","   sp3-4  ","          ","          ","          "/
      DATA LABEL(1:7,-2) &
     &/"   sp2-1  ","   sp2-2  ","   sp2-3  ","          ","          ","          ","          "/
      DATA LABEL(1:7,-1) &
     &/"    sp-1  ","    sp-2  ","          ","          ","          ","          ","          "/
      DATA LABEL(1:7, 0) &
     &/"      s   ","          ","          ","          ","          ","          ","          "/
      DATA LABEL(1:7, 1) &
     &/"     py   ","     pz   ","     px   ","          ","          ","          ","          "/
      DATA LABEL(1:7, 2) &
     &/"    dxy   ","    dyz   ","    dz2   ","    dxz   ","   dx2-y2 ","          ","          "/
      DATA LABEL(1:7, 3) &
     &/"fy(3x2-y2)","    fxyz  ","    fyz2  ","    fz3   ","    fxz2  "," fz(x2-y2)","fx(x2-3y2)"/
    character(LEN=40) :: my_group, my_subgroup

    my_group    = GRP_RESULTS;  if (present(group))    my_group    = trim(group)
    my_subgroup = GRP_LOCPROJ; if (present(subgroup)) my_subgroup = subgroup

! create / open groups
    call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",2815)
    call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",2816)
    call vh5_error(vh5_group_open_or_create(subgroupid, trim("parameters"), paramgroupid),"vhdf5.F",2817)

! allocate lists of strings
    nproj=SIZE(LPRJ_COVL,4)
    allocate(radial_type_list(nproj))
    allocate(ang_type_list(nproj))
    allocate(lprj_coord_list(3,nproj))
! write string info in lists first
    do ip=1,nproj
      radial_type_list(ip) = trim(adjustl(prkind(lprj_functions(ip)%radial_type)))
      ang_type_list(ip) = trim(adjustl(LABEL(lprj_functions(ip)%m,lprj_functions(ip)%l)))
      lprj_coord_list(:,ip) = lprj_functions(ip)%r
    enddo

! write locproj parameters
    call vh5_error(vh5_write(paramgroupid, "site", lprj_functions(:)%poscar_site),"vhdf5.F",2832)
    call vh5_error(vh5_write(paramgroupid, "coordinates", lprj_coord_list),"vhdf5.F",2833)
    call vh5_error(vh5_write(paramgroupid, "rad_type", radial_type_list),"vhdf5.F",2834)
    call vh5_error(vh5_write(paramgroupid, "ang_type", ang_type_list),"vhdf5.F",2835)
    call vh5_error(vh5_write(paramgroupid, "orbital", lprj_functions(:)%n),"vhdf5.F",2836)
    call vh5_error(vh5_write(paramgroupid, "alpha", lprj_functions(:)%za),"vhdf5.F",2837)
    call vh5_error(vh5_write(paramgroupid, "species_number", lprj_functions(:)%species),"vhdf5.F",2838)
    if (present(lnoncol)) then
      call vh5_error(vh5_write(paramgroupid, "lnoncollinear", lnoncol),"vhdf5.F",2840)
    endif

! the locproj data itself and a small format header
    call vh5_error(vh5_write(subgroupid, "data", CONJG(lprj_covl)),"vhdf5.F",2844)
    call vh5_error(vh5_write(subgroupid, "format", "[proj_index, spin_index, kpts_index, band_index]"),"vhdf5.F",2845)

    call vh5_error(vh5_group_close_writing(paramgroupid),"vhdf5.F",2847)
    call vh5_error(vh5_group_close_writing(subgroupid),"vhdf5.F",2848)
    call vh5_error(vh5_group_close_writing(groupid),"vhdf5.F",2849)
    call vh5_error(vh5_file_flush(fileid),"vhdf5.F",2850)
  end subroutine

!> write the k-points along path to a file
  subroutine kpoints_along_path_hdf5write(fileid,kpoints_ikibz,kpoints_dists,kpoints_coord,group,subgroup)
     integer(HID_T), intent(in) :: fileid !< hdf5 file handle
     real(q) :: kpoints_dists(:)
     integer :: kpoints_ikibz(:)
     real(q) :: kpoints_coord(:,:)
     character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
     character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_PROJECTORS)
!----------------------------------------------------------------------
     integer(HID_T) :: groupid, subgroupid
     character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
!----------------------------------------------------------------------
     my_group    = GRP_RESULTS;         if (present(group))    my_group    = trim(group)
     my_subgroup = SUBGRP_KPOINTS_LINE; if (present(subgroup)) my_subgroup = trim(subgroup)

     call vh5_error(vh5_group_open_or_create(fileid, trim(my_group), groupid),"vhdf5.F",2868)
     call vh5_error(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid),"vhdf5.F",2869)

     call vh5_error(vh5_write(subgroupid, "kpoints_ikibz", kpoints_ikibz),"vhdf5.F",2871)
     call vh5_error(vh5_write(subgroupid, "kpoints_dists", kpoints_dists),"vhdf5.F",2872)
     call vh5_error(vh5_write(subgroupid, "kpoint_coords", kpoints_coord),"vhdf5.F",2873)
  end subroutine kpoints_along_path_hdf5write

end module vhdf5
# 2888

