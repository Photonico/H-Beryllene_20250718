#include "symbol.inc"

module elphon_derivative
   use prec
   use base, only: in_struct
   use msymmetry, only: prim_cell_t
   use poscar_struct_def, only: type_info
   use wave_struct_def, only: wavedes
   use mgrid_struct_def, only: grid_3d, transit
   use tutor, only: vtutor

#ifdef VASP_HDF5
   use hdf5, only: hid_t
   use vhdf5_base, only: &
      vh5_file_create_or_overwrite, &
      vh5_file_close_writing, &
      vh5_write, &
      vh5_error
#endif

   implicit none
   private

   public elphon_deriv_settings, elphon_deriv_reader
   public elphon_deriv_io

   ! _test_mani temporary
   public calc_elphon_derivatives

   type elphon_deriv_settings
      logical :: generate = .false.
      real(q) :: prim_latt(3, 3)
      logical :: use_prim_latt = .false.
      integer :: fft_mesh(3) = [0, 0, 0]
   end type

   type grid_interpolator_t
      type(grid_3d), public :: source_grid
      type(grid_3d), public :: target_grid

      real(q), private :: recip_target_lattice(3, 3)
      real(q), private :: recip_source_lattice(3, 3)
      complex(q), allocatable, private :: buffer(:)
      integer, private :: target_to_source(3, 3)
      logical, private :: is_aligned
      type(transit), private :: aligned_transit
      integer, allocatable, private :: general_transit(:)
      integer, allocatable, private :: plane_wave_weight(:)
   contains
      procedure, public :: init => grid_interpolator_init
      procedure, public :: interpolate => grid_interpolator_interpolate
      procedure, private :: set_general_transit => grid_interpolator_set_general_transit

      ! _test_mani write finalizer
   end type

   type min_image_mapper_t
      real(q) :: cart_to_source(3, 3)
      real(q) :: standard_to_cart(3, 3)
      real(q) :: source_to_standard(3, 3)
      real(q) :: shift_vectors(3, 8)
   contains
      procedure, public :: init => min_image_mapper_init
      procedure, public :: get_shift => min_image_mapper_get_shift
   end type

contains

   subroutine min_image_mapper_init(this, source_lattice)
      class(min_image_mapper_t), intent(out) :: this
      real(q), intent(in) :: source_lattice(3, 3)
      real(q) :: std_latt(3, 3)
      integer :: shift_idx
      integer :: x, y, z
      integer :: i_dum
      real(q) :: r_dum(6)

      std_latt = source_lattice
      call lattyp(std_latt(:, 1), std_latt(:, 2), std_latt(:, 3), i_dum, r_dum, -1)
      this%standard_to_cart = std_latt
      this%cart_to_source = matrix_inverse_3x3_r(source_lattice)
      this%source_to_standard = matmul(matrix_inverse_3x3_r(std_latt), source_lattice)

      shift_idx = 0
      do z = -1, 0
         do y = -1, 0
            do x = -1, 0
               shift_idx = shift_idx + 1
               this%shift_vectors(:, shift_idx) = matmul(std_latt, [x, y, z])
            enddo
         enddo
      enddo
   end subroutine

   function min_image_mapper_get_shift(this, point, square_dist) result(shift)
      integer :: shift(3)
      class(min_image_mapper_t), intent(in) :: this
      real(q), intent(in) :: point(3)
      real(q), optional, intent(out) :: square_dist

      real(q) :: trafo_point(3)
      real(q) :: shifted_point(3)
      real(q) :: min_dist, cur_dist
      integer :: shift_idx
      integer :: min_idx

      trafo_point = matmul(this%source_to_standard, point)
      trafo_point = modulo(trafo_point, 1._q)
      trafo_point = matmul(this%standard_to_cart, trafo_point)

      min_dist = huge(1._q)
      do shift_idx = 1, 8
         shifted_point = trafo_point + this%shift_vectors(:, shift_idx)
         cur_dist = dot_product(shifted_point, shifted_point)
         if (cur_dist < min_dist) then
            min_dist = cur_dist
            min_idx = shift_idx
         endif
      enddo

      shifted_point = trafo_point + this%shift_vectors(:, min_idx)
      shifted_point = matmul(this%cart_to_source, shifted_point)
      shift = nint(shifted_point - point)

      if (present(square_dist)) square_dist = min_dist
   end function

   subroutine elphon_deriv_reader(incar, settings)
      use incar_reader, only: incar_file, process_incar

      type(incar_file), intent(inout) :: incar
      type(elphon_deriv_settings), intent(out) :: settings

      integer :: ierr

      call process_incar(incar, "elph_pot_generate", settings%generate)
      call process_incar(incar, "elph_pot_lattice", settings%prim_latt, ierr)
      settings%use_prim_latt = ierr == 0
      call process_incar(incar, "elph_pot_fft_mesh", settings%fft_mesh)

#ifndef VASP_HDF5
      if (settings%generate) call vtutor%error( &
         "In order to use ELPH_POT_GENERATE=True, VASP needs to be compiled with HDF5 support.")
#endif
   end subroutine

   subroutine elphon_deriv_io(io, settings, prim_cell, t_info, wdes, &
      Dij, Qij, V_grid, dV_du, dDij_du, force_const, polar_data)
      use elphon_util, only: polar_data_t

      type(in_struct), intent(in) :: io
      type(elphon_deriv_settings), intent(in) :: settings
      type(prim_cell_t), target, intent(in) :: prim_cell
      type(type_info), intent(in) :: t_info
      type(wavedes), intent(in) :: wdes
      OVERLAP, intent(in) :: Dij(:, :, :, :)
      OVERLAP, intent(in) :: Qij(:, :, :, :)
      type(grid_3d), intent(in) :: V_grid
      complex(q), intent(in) :: dV_du(:, :, :, :)
      OVERLAP, intent(in) :: dDij_du(:, :, :, :, :, :)
      real(q), intent(in) :: force_const(:, :)
      type(polar_data_t), intent(in) :: polar_data

      type(prim_cell_t), pointer :: prim_cell_ptr

      if (settings%use_prim_latt) then
         allocate(prim_cell_ptr)
         call prim_cell_ptr%init_from_latt(prim_cell, settings%prim_latt)
      else
         prim_cell_ptr => prim_cell
      endif

#ifdef VASP_HDF5
      call write_phelel_params( &
         io, prim_cell_ptr, t_info, wdes, settings%fft_mesh, Dij, Qij, &
         V_grid, dV_du, dDij_du, force_const, polar_data)
#endif

      if (.not. associated(prim_cell_ptr, prim_cell)) deallocate(prim_cell_ptr)
   end subroutine

#ifdef VASP_HDF5
   subroutine write_phelel_params(io, prim_cell, t_info, wdes, fft_mesh, &
      Dij, Qij, V_grid, dV_du, dDij_du, force_const, polar_data, filename)
      use elphon_util, only: polar_data_t

      type(in_struct), intent(in) :: io
      type(prim_cell_t), intent(in) :: prim_cell
      type(type_info), intent(in) :: t_info
      type(wavedes), intent(in) :: wdes
      integer, intent(in) :: fft_mesh(3)
      OVERLAP, intent(in) :: Dij(:, :, :, :)
      OVERLAP, intent(in) :: Qij(:, :, :, :)
      type(grid_3d), intent(in) :: V_grid
      complex(q), intent(in) :: dV_du(:, :, :, :)
      OVERLAP, intent(in) :: dDij_du(:, :, :, :, :, :)
      real(q), intent(in) :: force_const(:, :)
      type(polar_data_t), intent(in) :: polar_data
      character(*), optional, intent(in) :: filename

      integer(hid_t) :: file_id
      character(:), allocatable :: act_filename
      integer :: supercell_matrix(3, 3)
      integer :: super_target_to_source(3, 3)
      integer :: prim_source_to_target(3, 3)
      type(grid_interpolator_t) :: grid_interpol
      integer, allocatable :: lattice_points(:, :)
      real(q), allocatable :: grid_points(:, :)
      integer :: interpol_fft_mesh(3)

      act_filename = "phelel_params.hdf5"
      if (present(filename)) act_filename = filename

      supercell_matrix = nint(prim_cell%super_to_prim)

      call find_target_grid( &
         supercell_matrix, fft_mesh, interpol_fft_mesh, &
         super_target_to_source, prim_source_to_target)

      call grid_interpol%init(V_grid, prim_cell%super_latt%a, super_target_to_source, interpol_fft_mesh)

      call generate_lattice_points(supercell_matrix, lattice_points)
      call generate_grid_points(lattice_points, fft_mesh, grid_points)

      ! _test_mani warn about FFT grid change

      VH5_CHECK(vh5_file_create_or_overwrite(act_filename, file_id))

      call write_primitive_cell_info(io, file_id, prim_cell, t_info)
      call write_paw_overlaps(io, file_id, prim_cell, wdes, Dij, Qij, dDij_du)
      call write_local_pot_derivative(io, file_id, &
         grid_interpol, grid_points, supercell_matrix, super_target_to_source, dV_du)
      call write_force_constants(io, file_id, force_const)
      call write_polar_data(io, file_id, polar_data)
      call write_grid_dimensions(io, file_id, &
         lattice_points, grid_points, fft_mesh, prim_cell%prim_to_super)

      VH5_CHECK(vh5_file_close_writing(file_id))
   end subroutine

   subroutine write_primitive_cell_info(io, group_id, prim_cell, t_info)
      type(in_struct), intent(in) :: io
      integer(hid_t), intent(in) :: group_id
      type(prim_cell_t), intent(in) :: prim_cell
      type(type_info), intent(in) :: t_info

      integer, allocatable :: prim_indices(:)
      real(q), allocatable :: prim_positions(:, :)

      if (io%iu0 < 0) return

      allocate(prim_indices(prim_cell%num_prim_atoms))
      prim_indices = prim_cell%super_idx(1, :)

      VH5_CHECK(vh5_write(group_id, "primitive_lattice", transpose(prim_cell%prim_latt%a)))
      VH5_CHECK(vh5_write(group_id, "supercell_lattice", transpose(prim_cell%super_latt%a)))
      VH5_CHECK(vh5_write(group_id, "unitcell_lattice", transpose(prim_cell%prim_latt%a)))
      VH5_CHECK(vh5_write(group_id, "supercell_matrix", transpose(nint(prim_cell%super_to_prim))))
      VH5_CHECK(vh5_write(group_id, "supercell_positions", prim_cell%posion_super))
      VH5_CHECK(vh5_write(group_id, "p2s_map", prim_indices-1))
      VH5_CHECK(vh5_write(group_id, "atom_indices_in_derivatives", prim_indices-1))
      VH5_CHECK(vh5_write(group_id, "primitive_masses", t_info%pomass(t_info%ityp(prim_indices))))

      allocate(prim_positions(3, prim_cell%num_prim_atoms))
      prim_positions = matmul(prim_cell%super_to_prim, prim_cell%posion_super(:, prim_indices))
      VH5_CHECK(vh5_write(group_id, "primitive_positions", prim_positions))
   end subroutine

   subroutine write_paw_overlaps(io, group_id, prim_cell, wdes, Dij, Qij, dDij_du)
      type(in_struct), intent(in) :: io
      integer(hid_t), intent(in) :: group_id
      type(prim_cell_t), intent(in) :: prim_cell
      type(wavedes), intent(in) :: wdes
      OVERLAP, intent(in) :: Dij(:, :, :, :)
      OVERLAP, intent(in) :: Qij(:, :, :, :)
      OVERLAP, intent(in) :: dDij_du(:, :, :, :, :, :)

      integer :: num_chan, ncdij
      integer :: num_prim_atoms, num_displacements
      integer :: total_size
      complex(q), allocatable :: overlap_out(:, :, :, :)
      complex(q), allocatable :: overlap_deriv_out(:, :, :, :, :, :)

      num_chan = size(Dij, 1)
      ncdij = size(Dij, 4)
      num_displacements = size(dDij_du, 5)
      num_prim_atoms = size(dDij_du, 6)

      allocate(overlap_out(num_chan, num_chan, prim_cell%num_prim_atoms, wdes%ncdij))
      total_size = size(overlap_out)

      call gather_overlap(Dij)
      if (io%iu0 >= 0) VH5_CHECK(vh5_write(group_id, "Dij", overlap_out))
      call gather_overlap(Qij)
      if (io%iu0 >= 0) VH5_CHECK(vh5_write(group_id, "qij", overlap_out))

      deallocate(overlap_out)
      allocate(overlap_deriv_out( &
         num_chan, num_chan, prim_cell%num_super_atoms, &
         num_displacements, num_prim_atoms, wdes%ncdij))

      call gather_overlap_derivative(dDij_du)
      if (io%iu0 >= 0) VH5_CHECK(vh5_write(group_id, "dDijdu", overlap_deriv_out))
      overlap_deriv_out = 0
      if (io%iu0 >= 0) VH5_CHECK(vh5_write(group_id, "dqijdu", overlap_deriv_out))

   contains

      subroutine gather_overlap(local_overlap)
         use wave, only: ni_global
         OVERLAP, target, intent(in) :: local_overlap(:, :, :, :)

         integer :: local_idx, global_idx
         integer :: prim_idx

         overlap_out = 0
         do prim_idx = 1, prim_cell%num_prim_atoms
            do local_idx = 1, wdes%nions
               global_idx = ni_global(local_idx, wdes%comm_inb)
               if (global_idx == prim_cell%super_idx(1, prim_idx)) then
                  overlap_out(:, :, prim_idx, :) = local_overlap(:, :, local_idx, :)
               endif
            enddo
         enddo
         CALLMPI(M_sum_z(wdes%comm_inb, overlap_out, size(overlap_out)))
      end subroutine

      subroutine gather_overlap_derivative(local_overlap)
         use wave, only: ni_global
         OVERLAP, target, intent(in) :: local_overlap(:, :, :, :, :, :)

         integer :: local_idx, global_idx
         integer :: icdij

         overlap_deriv_out = 0
         do icdij = 1, wdes%ncdij
            do local_idx = 1, wdes%nions
               global_idx = ni_global(local_idx, wdes%comm_inb)
               overlap_deriv_out(:, :, global_idx, :, :, icdij) = &
                  local_overlap(:, :, local_idx, icdij, :, :)
            enddo
         enddo
         CALLMPI(M_sum_z(wdes%comm_inb, overlap_deriv_out, size(overlap_deriv_out)))
      end subroutine
   end subroutine

   subroutine write_local_pot_derivative(io, group_id, &
      grid_interpol, grid_points, supercell_matrix, super_target_to_source, dV_du)

      type(in_struct), intent(in) :: io
      integer(hid_t), intent(in) :: group_id
      type(grid_interpolator_t), intent(inout) :: grid_interpol
      real(q), intent(in) :: grid_points(:, :)
      integer, intent(in) :: supercell_matrix(3, 3)
      integer, intent(in) :: super_target_to_source(3, 3)
      complex(q), intent(in) :: dV_du(:, :, :, :)

      integer :: num_prim_atoms, ncdij
      integer :: prim_idx, icdij, cart_idx
      complex(q), allocatable :: interpol_buffer(:)
      complex(q), allocatable :: dV_du_out(:, :, :, :)

      num_prim_atoms = size(dV_du, 4)
      ncdij = size(dV_du, 2)

      allocate(interpol_buffer(grid_interpol%target_grid%mplwv))
      allocate(dV_du_out(size(grid_points, 2), 3, num_prim_atoms, ncdij))

      ! Iterate over all degrees of freedom
      do prim_idx = 1, num_prim_atoms
         do cart_idx = 1, 3
            do icdij = 1, ncdij
               ! Interpolate the data to the interpolation grid
               call grid_interpol%interpolate( &
                  dV_du(:, icdij, cart_idx, prim_idx), interpol_buffer)
               ! Transfer potential to output grid at relevant grid points
               call transfer_to_grid_points( &
                  grid_interpol%target_grid, grid_points, supercell_matrix, super_target_to_source, &
                  interpol_buffer, dV_du_out(:, cart_idx, prim_idx, icdij))
            enddo
         enddo
      enddo

      if (io%iu0 >= 0) VH5_CHECK(vh5_write(group_id, "dVdu", dV_du_out))
   end subroutine

   subroutine write_force_constants(io, group_id, force_const)
      use string, only: str

      type(in_struct), intent(in) :: io
      integer(hid_t), intent(in) :: group_id
      real(q), intent(in) :: force_const(:, :)

      real(q), allocatable :: force_const_blocks(:, :, :, :)
      integer :: array_size
      integer :: num_atoms
      integer :: atom_idx_1, atom_idx_2

      array_size = size(force_const, 1)
      num_atoms = array_size / 3
      if (num_atoms * 3 /= array_size .or. array_size /= size(force_const, 2)) &
         call vtutor%bug("IFC matrix has wrong size: " // &
         str(shape(force_const)) // " " // str(num_atoms), __FILE__, __LINE__)

      allocate(force_const_blocks(3, 3, num_atoms, num_atoms))
      do atom_idx_2 = 1, num_atoms
         do atom_idx_1 = 1, num_atoms
            force_const_blocks(:, :, atom_idx_1, atom_idx_2) = &
               force_const(atom_idx_1*3-2 : atom_idx_1*3, atom_idx_2*3-2 : atom_idx_2*3)
         enddo
      enddo

      if (io%iu0 >= 0) VH5_CHECK(vh5_write(group_id, "force_constants", force_const_blocks))
   end subroutine

   subroutine write_polar_data(io, group_id, polar_data)
      use elphon_util, only: polar_data_t

      type(in_struct), intent(in) :: io
      integer(hid_t), intent(in) :: group_id
      type(polar_data_t), intent(in) :: polar_data

      real(q), allocatable :: buffer(:, :, :)
      integer :: num_prim, prim_idx

      if (allocated(polar_data%dielectric_tensor)) then
         VH5_CHECK(vh5_write(group_id, "dielectric_constant", transpose(polar_data%dielectric_tensor)))
      endif
      if (allocated(polar_data%born_eff_charges)) then
         num_prim = size(polar_data%born_eff_charges, 3)
         allocate(buffer(3, 3, num_prim))
         do prim_idx = 1, num_prim
            buffer(:, :, prim_idx) = transpose(polar_data%born_eff_charges(:, :, prim_idx))
         enddo
         VH5_CHECK(vh5_write(group_id, "born_effective_charges", buffer))
      endif
   end subroutine

   subroutine write_grid_dimensions(io, group_id, &
      lattice_points, grid_points, prim_fft_mesh, prim_to_super)

      type(in_struct), intent(in) :: io
      integer(hid_t), intent(in) :: group_id
      integer, intent(in) :: lattice_points(:, :)
      real(q), intent(in) :: grid_points(:, :)
      integer, intent(in) :: prim_fft_mesh(3)
      real(q), intent(in) :: prim_to_super(3, 3)

      real(q), allocatable :: grid_points_out(:, :)

      allocate(grid_points_out(3, size(grid_points, 2)))

      grid_points_out = matmul(prim_to_super, grid_points)

      VH5_CHECK(vh5_write(group_id, "lattice_point", lattice_points))
      VH5_CHECK(vh5_write(group_id, "grid_point", grid_points_out))
      VH5_CHECK(vh5_write(group_id, "FFT_mesh", prim_fft_mesh))
   end subroutine
#endif

   subroutine find_target_grid(supercell_matrix, prim_fft_mesh, &
      target_fft_mesh, super_target_to_source, prim_source_to_target)
      use mkpoints, only: SmithNormalForm

      integer, intent(in) :: supercell_matrix(3, 3)
      integer, intent(in) :: prim_fft_mesh(3)
      integer, intent(out) :: target_fft_mesh(3)
      integer, intent(out) :: super_target_to_source(3, 3)
      integer, intent(out) :: prim_source_to_target(3, 3)

      integer :: i
      integer :: matrix(3, 3)
      integer :: normal_form(3, 3)

      if (is_diagonal_and_positive(supercell_matrix)) then
         super_target_to_source = 0
         prim_source_to_target = 0
         do i = 1, 3
            target_fft_mesh(i) = supercell_matrix(i, i) * prim_fft_mesh(i)
            super_target_to_source(i, i) = 1
            prim_source_to_target(i, i) = 1
         enddo
      else
         do i = 1, 3
            matrix(:, i) = prim_fft_mesh(:) * supercell_matrix(:, i)
         enddo
         call SmithNormalForm(matrix, prim_source_to_target, normal_form, super_target_to_source)
         do i = 1, 3
            target_fft_mesh(i) = normal_form(i, i)
         enddo
      endif

   contains
      pure logical function is_diagonal_and_positive(m)
         integer, intent(in) :: m(3, 3)
         integer :: a, b
         is_diagonal_and_positive = .false.
         do a = 1, 3
            do b = 1, 3
               if (a == b) then
                  if (m(a, b) < 0) return
               else
                  if (m(a, b) /= 0) return
               endif
            enddo
         enddo
         is_diagonal_and_positive = .true.
      end function
   end subroutine

   subroutine grid_interpolator_init(this, initial_grid, initial_lattice, target_to_initial, target_fft_mesh)
      use mgrid, only: inilgrd, gen_rc_grid, gen_rc_sub_grid

      class(grid_interpolator_t), intent(out) :: this
      type(grid_3d), intent(in) :: initial_grid
      real(q), intent(in) :: initial_lattice(3, 3)
      integer, intent(in) :: target_to_initial(3, 3)
      integer, intent(in) :: target_fft_mesh(3)

      integer, parameter :: identity(3, 3) = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])

      this%source_grid = initial_grid
      this%recip_target_lattice = transpose(matrix_inverse_3x3_r( &
         matmul(initial_lattice, target_to_initial)))
      this%recip_source_lattice = transpose(matrix_inverse_3x3_r(initial_lattice))

      associate(tgrid => this%target_grid)
         ! Allocate storage for the FFT.
         allocate(this%buffer(initial_grid%mplwv))
         this%target_to_source = target_to_initial

         ! Set target FFT mesh.
         tgrid%ngptar = target_fft_mesh
         ! call FFTCHK(tgrid%ngptar)

         ! Initialize target grid.
         call inilgrd(tgrid%ngptar(1), tgrid%ngptar(2), tgrid%ngptar(3), tgrid)

         this%is_aligned = all(target_to_initial == identity)

         if (this%is_aligned) then
            ! Point to the same communicators as the initial grid.
            tgrid%comm => initial_grid%comm
            tgrid%comm_kin => initial_grid%comm_kin
            tgrid%comm_kinter => initial_grid%comm_kinter
            ! Create target grid as sub grid of the initial grid.
            ! This also sets the transition table between the two.
            call gen_rc_sub_grid(tgrid, initial_grid, this%aligned_transit, .true., .true.)
         else
            ! Create new communicator for target grid.
            call init_grid_comm(tgrid)
            ! Create target grid as a separate grid.
            call gen_rc_grid(tgrid)
            ! Create map between the two grids.
            call this%set_general_transit()
         endif

         ! Setup communication part of target grid.
         CALLMPI(mapset(tgrid))
      end associate

   contains

      subroutine init_grid_comm(grid)
#if defined(MPI) || defined(MPI_CHAIN)
         use mpi, only: MPI_COMM_SELF
         use mpimy, only: M_initc
#endif
         type(grid_3d), intent(inout) :: grid

         ! Initialize communicator of target grid to be fully serial.
         ! _test_mani this can be improved in the future
         allocate(grid%comm)
#if defined(MPI) || defined(MPI_CHAIN)
         grid%comm%mpi_comm = MPI_COMM_SELF
         CALL M_initc(grid%comm)
#endif
         ! Just point everything to the same serial communicator.
         ! _test_mani this can be improved in the future
         grid%comm_kin => grid%comm
         grid%comm_kinter => grid%comm
      end subroutine
   end subroutine

   subroutine grid_interpolator_interpolate(this, initial_values, interpol_values)
      class(grid_interpolator_t), intent(inout) :: this
      complex(q), intent(in) :: initial_values(:)
      complex(q), intent(out) :: interpol_values(:)

      integer :: pw_idx, mapped_pw_idx
      integer :: start_idx
      complex(q) :: coeff

      ! Copy initial values into FFT buffer
      this%buffer(:) = initial_values(:)
      ! FFT values into reciprocal space
      call FFT3D(this%buffer, this%source_grid, -1)
      ! Scale with number of points
      this%buffer = this%buffer / this%source_grid%nplwv

      ! Initialize plane-wave components to zero
      interpol_values = 0
      ! Copy plane-wave components into buffer of target FFT grid
      if (this%is_aligned) then
         ! Simple case, just copy plane-wave components
         call cp_grid(this%source_grid, this%target_grid, &
            this%aligned_transit, this%buffer, interpol_values)
      else
         start_idx = 1
         if (this%source_grid%real2cplx) start_idx = -this%source_grid%rc%np
         do pw_idx = start_idx, this%source_grid%rc%np
            mapped_pw_idx = this%general_transit(pw_idx)
            if (mapped_pw_idx /= 0) then
               coeff = this%buffer(abs(pw_idx))
               if (pw_idx < 0 .neqv. mapped_pw_idx < 0) coeff = conjg(coeff)
               interpol_values(abs(mapped_pw_idx)) = interpol_values(abs(mapped_pw_idx)) + coeff
            endif
         enddo

         do pw_idx = 1, this%target_grid%rc%np
            if (this%plane_wave_weight(pw_idx) > 1) then
               interpol_values(pw_idx) = interpol_values(pw_idx) / this%plane_wave_weight(pw_idx)
            endif
         enddo

         ! Accumulate all plane-wave components onto one rank
         ! _test_mani right now this is done for all ranks
         CALLMPI(M_sum_z(this%source_grid%comm, interpol_values, this%target_grid%mplwv))
      endif

      ! FFT values into real space
      call FFT3D(interpol_values, this%target_grid, 1)
   end subroutine

   subroutine grid_interpolator_set_general_transit(this)
      class(grid_interpolator_t), intent(inout) :: this

      integer, allocatable :: g_map(:, :, :)
      integer :: lower_bounds_1(3)
      integer :: lower_bounds_2(3)
      integer :: upper_bounds_2(3)
      integer :: column_idx
      integer :: x_idx, y_idx, z_idx
      integer :: g_vec(3)
      integer :: shift(3)
      integer :: pw_idx
      integer :: gz
      real(q) :: min_dist
      real(q), allocatable :: distances(:)
      type(min_image_mapper_t) :: min_image_mapper

      associate (gr1 => this%source_grid, gr2 => this%target_grid)
         ! Get min and max G-vector elements of target grid in each direction.
         lower_bounds_1 = [minval(gr1%lpctx), minval(gr1%lpcty), minval(gr1%lpctz)]
         lower_bounds_2 = [minval(gr2%lpctx), minval(gr2%lpcty), minval(gr2%lpctz)]
         upper_bounds_2 = [maxval(gr2%lpctx), maxval(gr2%lpcty), maxval(gr2%lpctz)]

         ! This array maps all G-vectors of the target grid to the corresponding
         ! plane-wave indices in the source grid. The indices of the map array
         ! are the G-vector components in direct coordinates.
         ! This mapping is done locally, i.e., no MPI communication takes place.
         ! _test_mani likely this has to change in the more general case.
         allocate(g_map( &
            lower_bounds_2(1) : upper_bounds_2(1), &
            lower_bounds_2(2) : upper_bounds_2(2), &
            lower_bounds_2(3) : upper_bounds_2(3)))
         ! G-vectors that do not exist in either grid are mapped to the value 0.
         g_map = 0

         allocate(distances(gr2%rc%np))
         distances = huge(1._q)

         ! We set this up purely to get minimum Cartesian distances for the G-vectors
         call min_image_mapper%init(this%recip_target_lattice * spread(gr2%ngptar, 1, 3))

         ! Loop over all G-vectors of the target grid
         pw_idx = 0
         do column_idx = 1, gr2%rc%ncol
            y_idx = gr2%rc%i2(column_idx)
            z_idx = gr2%rc%i3(column_idx)
            do x_idx = 1, gr2%rc%nrow
               ! Get plane-wave index of target grid
               pw_idx = pw_idx + 1
               ! Get G-vector of target grid in direct coordinates
               g_vec = [gr2%lpctx(x_idx), gr2%lpcty(y_idx), gr2%lpctz(z_idx)]
               ! Map G-vector to plane-wave index
               g_map(g_vec(1), g_vec(2), g_vec(3)) = pw_idx
               ! Set minimum distance for this G-vector
               shift = min_image_mapper%get_shift(real(g_vec, q) / gr2%ngptar, distances(pw_idx))

               ! Map -G explicitly if half-grid mode is used
               if (gr2%real2cplx .and. g_vec(3) > 0) then
                  g_vec = -g_vec
                  g_vec = modulo(g_vec - lower_bounds_2, gr2%ngptar) + lower_bounds_2
                  ! Sometimes we map back to positive z so we need to check
                  if (g_vec(3) < 0) g_map(g_vec(1), g_vec(2), g_vec(3)) = -pw_idx
               endif
            enddo
         enddo

         allocate(this%plane_wave_weight(gr2%rc%np))
         this%plane_wave_weight = 0

         if (gr1%real2cplx) then
            allocate(this%general_transit(-gr1%rc%np : gr1%rc%np))
         else
            allocate(this%general_transit(gr1%rc%np))
         endif
         this%general_transit = 0

         call min_image_mapper%init(this%recip_source_lattice * spread(gr1%ngptar, 1, 3))

         ! Loop over all G-vectors of the source grid
         pw_idx = 0
         do column_idx = 1, gr1%rc%ncol
            y_idx = gr1%rc%i2(column_idx)
            z_idx = gr1%rc%i3(column_idx)
            do x_idx = 1, gr1%rc%nrow
               ! Get plane-wave index of source grid
               pw_idx = pw_idx + 1
               ! Get G-vector of source grid in direct coordinates
               g_vec = [gr1%lpctx(x_idx), gr1%lpcty(y_idx), gr1%lpctz(z_idx)]
               ! Shift G-vector to have smallest distance from origin in Cartesian coordinates
               shift = min_image_mapper%get_shift(real(g_vec, q) / gr1%ngptar, min_dist)
               g_vec = g_vec + shift * gr1%ngptar
               ! Get corresponding G-vector in target grid
               g_vec = matmul(g_vec, this%target_to_source)
               ! Map plane-wave index of source grid to plane-wave index of target grid
               this%general_transit(pw_idx) = map_g_vec(g_vec)

               ! Map -G explicitly if half-grid mode is used
               if (gr1%real2cplx .and. gr1%lpctz(z_idx) > 0) then
                  gz = -gr1%lpctz(z_idx)
                  gz = modulo(gz - lower_bounds_1(3), gr1%ngz) + lower_bounds_1(3)
                  ! Sometimes we map back to positive z so we need to check
                  if (gz < 0) this%general_transit(-pw_idx) = map_g_vec(-g_vec)
               endif
            enddo
         enddo
      end associate

      CALLMPI(M_sum_i(this%source_grid%comm, this%plane_wave_weight, this%target_grid%rc%np))

   contains
      integer function map_g_vec(g_in) result(pw_idx)
         integer, intent(in) :: g_in(3)
         integer :: g(3)
         integer :: abs_idx

         g = modulo(g_in - lower_bounds_2, this%target_grid%ngptar) + lower_bounds_2
         pw_idx = g_map(g(1), g(2), g(3))

         abs_idx = abs(pw_idx)
         if (min_dist < distances(abs_idx) + 1e-7_q) then
            this%plane_wave_weight(abs_idx) = this%plane_wave_weight(abs_idx) + 1
         else
            pw_idx = 0
         endif
      end function
   end subroutine

   pure integer function det_3x3_i(a) result(det)
      integer, intent(in) :: a(3, 3)

      det   = a(1, 1) * a(2, 2) * a(3, 3) &
            + a(2, 1) * a(3, 2) * a(1, 3) &
            + a(1, 2) * a(2, 3) * a(3, 1) &
            - a(1, 3) * a(2, 2) * a(3, 1) &
            - a(2, 1) * a(1, 2) * a(3, 3) &
            - a(3, 2) * a(2, 3) * a(1, 1)
   end function

   pure function matrix_inverse_3x3_i(a) result(b)
      integer, intent(in) :: a(3, 3)
      integer :: b(3, 3)
      integer :: det_inv

      ! calculate inverse determinant.
      det_inv = 1 / det_3x3_i(a)

      ! calculate inverse of the matrix.
      b(1, 1) = +det_inv * (a(2, 2) * a(3, 3) - a(2, 3) * a(3, 2))
      b(2, 1) = -det_inv * (a(2, 1) * a(3, 3) - a(2, 3) * a(3, 1))
      b(3, 1) = +det_inv * (a(2, 1) * a(3, 2) - a(2, 2) * a(3, 1))
      b(1, 2) = -det_inv * (a(1, 2) * a(3, 3) - a(1, 3) * a(3, 2))
      b(2, 2) = +det_inv * (a(1, 1) * a(3, 3) - a(1, 3) * a(3, 1))
      b(3, 2) = -det_inv * (a(1, 1) * a(3, 2) - a(1, 2) * a(3, 1))
      b(1, 3) = +det_inv * (a(1, 2) * a(2, 3) - a(1, 3) * a(2, 2))
      b(2, 3) = -det_inv * (a(1, 1) * a(2, 3) - a(1, 3) * a(2, 1))
      b(3, 3) = +det_inv * (a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1))
   end function

   pure real(q) function det_3x3_r(a) result(det)
      real(q), intent(in) :: a(3, 3)

      det   = a(1, 1) * a(2, 2) * a(3, 3) &
            + a(2, 1) * a(3, 2) * a(1, 3) &
            + a(1, 2) * a(2, 3) * a(3, 1) &
            - a(1, 3) * a(2, 2) * a(3, 1) &
            - a(2, 1) * a(1, 2) * a(3, 3) &
            - a(3, 2) * a(2, 3) * a(1, 1)
   end function

   pure function matrix_inverse_3x3_r(a) result(b)
      real(q), intent(in) :: a(3, 3)
      real(q) :: b(3, 3)
      real(q) :: det_inv

      ! calculate inverse determinant.
      det_inv = 1._q / det_3x3_r(a)

      ! calculate inverse of the matrix.
      b(1, 1) = +det_inv * (a(2, 2) * a(3, 3) - a(2, 3) * a(3, 2))
      b(2, 1) = -det_inv * (a(2, 1) * a(3, 3) - a(2, 3) * a(3, 1))
      b(3, 1) = +det_inv * (a(2, 1) * a(3, 2) - a(2, 2) * a(3, 1))
      b(1, 2) = -det_inv * (a(1, 2) * a(3, 3) - a(1, 3) * a(3, 2))
      b(2, 2) = +det_inv * (a(1, 1) * a(3, 3) - a(1, 3) * a(3, 1))
      b(3, 2) = -det_inv * (a(1, 1) * a(3, 2) - a(1, 2) * a(3, 1))
      b(1, 3) = +det_inv * (a(1, 2) * a(2, 3) - a(1, 3) * a(2, 2))
      b(2, 3) = -det_inv * (a(1, 1) * a(2, 3) - a(1, 3) * a(2, 1))
      b(3, 3) = +det_inv * (a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1))
   end function

   subroutine generate_lattice_points(supercell_matrix, lattice_points)
      use mkpoints, only: HermiteNormalForm
      integer, intent(in) :: supercell_matrix(3, 3)
      integer, allocatable, intent(out) :: lattice_points(:, :)

      integer :: normal_form(3, 3), tmp(3, 3)
      integer :: supercell_size(3)
      integer :: point_idx
      integer :: x, y, z
      integer :: i

      call HermiteNormalForm(supercell_matrix, normal_form, tmp)
      do i = 1, 3
         supercell_size(i) = normal_form(i, i)
      enddo

      allocate(lattice_points(3, product(supercell_size)))

      point_idx = 0
      do x = 0, supercell_size(1) - 1
         do y = 0, supercell_size(2) - 1
            do z = 0, supercell_size(3) - 1
               point_idx = point_idx + 1
               lattice_points(:, point_idx) = [x, y, z]
            enddo
         enddo
      enddo
   end subroutine

   subroutine generate_grid_points(lattice_points, fft_mesh, grid_points)
      integer, intent(in) :: lattice_points(:, :)
      integer, intent(in) :: fft_mesh(3)
      real(q), allocatable, intent(out) :: grid_points(:, :)

      integer :: i1, i2, i3
      integer :: point_idx
      integer :: lattice_point_count
      integer :: lattice_point_idx

      lattice_point_count = size(lattice_points, 2)

      allocate(grid_points(3, lattice_point_count * product(fft_mesh)))

      point_idx = 0
      do lattice_point_idx = 1, lattice_point_count
         do i3 = 0, fft_mesh(3) - 1
            do i2 = 0, fft_mesh(2) - 1
               do i1 = 0, fft_mesh(1) - 1
                  point_idx = point_idx + 1
                  grid_points(:, point_idx) = &
                     lattice_points(:, lattice_point_idx) + &
                     [real(q):: i1, i2, i3] / fft_mesh
               enddo
            enddo
         enddo
      enddo
   end subroutine

   subroutine transfer_to_grid_points(grid, grid_points, &
      supercell_matrix, super_target_to_source, pot_in, pot_out)

      type(grid_3d), intent(in) :: grid
      real(q), intent(in) :: grid_points(:, :)
      integer, intent(in) :: supercell_matrix(3, 3)
      integer, intent(in) :: super_target_to_source(3, 3)
      complex(q), intent(in) :: pot_in(:)
      complex(q), intent(out) :: pot_out(:)

      integer :: num_out_points, point_idx
      integer :: fast_idx, slow_idx_1, slow_idx_2
      integer :: column_idx, pot_idx
      real(q) :: vec(3)
      integer :: idx_vec(3)
      real(q) :: trafo(3, 3)

      num_out_points = size(pot_out)

      trafo = matrix_inverse_3x3_r(real(matmul(supercell_matrix, super_target_to_source), q))

      fast_idx = grid%rl%nfast
      slow_idx_1 = modulo(fast_idx, 3) + 1
      slow_idx_2 = modulo(slow_idx_1, 3) + 1

      pot_out = 0

      do point_idx = 1, num_out_points
         vec = matmul(trafo, grid_points(:, point_idx))
         idx_vec = modulo(nint(vec * grid%ngptar), grid%ngptar)

         column_idx = grid%rl%index(idx_vec(slow_idx_1), idx_vec(slow_idx_2))
         if (column_idx == 0) cycle
         pot_idx = idx_vec(fast_idx) + (column_idx - 1) * grid%ngptar(fast_idx) + 1
         pot_out(point_idx) = get_val(pot_in, pot_idx)
      enddo

      CALLMPI(M_sum_z(grid%comm, pot_out, num_out_points))

   contains

      pure RGRID function get_val(arr, complex_idx)
         complex(q), intent(in) :: arr(:)
         integer, intent(in) :: complex_idx
#ifdef realmode
         integer :: array_idx

         array_idx = (complex_idx + 1) / 2
         if (mod(complex_idx, 2) == 1) then
            get_val = real(arr(array_idx), q)
         else
            get_val = aimag(arr(array_idx))
         endif
#else
         get_val = arr(complex_idx)
#endif
      end function
   end subroutine

   ! _test_mani this should be removed when the new finite-difference code is here
   subroutine calc_elphon_derivatives(findif_progress, V, Dij, dV_du, dDij_du)
      use finite_differences, only: findif_progress_t

      type(findif_progress_t), intent(in) :: findif_progress
      complex(q), intent(in) :: V(:, :)
      OVERLAP, intent(in) :: Dij(:, :, :, :)
      complex(q), intent(out) :: dV_du(:, :, :, :)
      OVERLAP, intent(out) :: dDij_du(:, :, :, :, :, :)

      integer :: prim_idx, cart_idx

      prim_idx = (findif_progress%cur_dof + 2) / 3
      cart_idx = findif_progress%cur_cart_idx

      if (findif_progress%cur_displ == 1) then
         dV_du(:, :, cart_idx, prim_idx) = V(:, :)
         dDij_du(:, :, :, :, cart_idx, prim_idx) = Dij(:, :, :, :)
      else
         dV_du(:, :, cart_idx, prim_idx) = &
            (dV_du(:, :, cart_idx, prim_idx) - V(:, :)) / &
            (2 * findif_progress%cur_step)

         dDij_du(:, :, :, :, cart_idx, prim_idx) = &
            (dDij_du(:, :, :, :, cart_idx, prim_idx) - Dij(:, :, :, :)) / &
            (2 * findif_progress%cur_step)
      endif
   end subroutine
end module
