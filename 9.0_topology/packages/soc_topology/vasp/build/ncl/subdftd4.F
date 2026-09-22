!> Wrapper around the DFTD4 library
!>
!> Translates the VASP structural data into a format compatible with the DFTD4
!> library. DFTD4 is an external library and needs to be installed separately.
!> Please check the [documentation](https://github.com/dftd4/dftd4) for more
!> information.
module vdwd4

   use prec, only: q
   use base, only: in_struct
#ifdef DFTD4
   use dftd4, only: rational_damping_param
#endif

   implicit none

   integer, parameter :: unknown_functional = 1, unknown_atom = -1, vdw_dftd4 = 13

#ifdef DFTD4
   logical, parameter :: dftd4_available = .true.
#else
   logical, parameter :: dftd4_available = .false.
   type rational_damping_param
   end type rational_damping_param
#endif

   type vdw_settings
      integer ivdw
      type(in_struct) io
      type(rational_damping_param) parameters
   end type vdw_settings

   type vdw_structure
      real(q), pointer :: cell(:,:) => null()
      integer, pointer :: types(:) => null()
      character(len=2), pointer :: elements(:) => null()
      real(q), pointer :: positions(:,:) => null()
   end type vdw_structure

   type vdw_interaction
      real(q) energy, stress(3,3)
      real(q), allocatable :: forces(:,:)
   end type vdw_interaction

   interface vdw_structure
      module procedure make_structure
   end interface vdw_structure

   private
   public vdw_settings, default_settings, vdw_forces_d4, unknown_functional, &
      vdw_interaction, vdw_structure, vdw_reader, vdw_dftd4, dftd4_available, &
      vdw_writer

contains

   !> Merge information from various types used by VASP into a single type used
   !> by the vdW routine.
   type(vdw_structure) function make_structure(latt_cur, dyn, t_info, elem) result (structure)
      use lattice, only: latt
      use poscar, only: dynamics, type_info
      type(latt), target, intent(in) :: latt_cur
      type(type_info), target, intent(in) :: t_info
      type(dynamics), target, intent(in) :: dyn
      character(len=2), target, intent(in) :: elem(:)
      structure%cell => latt_cur%a
      structure%types => t_info%iTyp
      structure%elements => elem
      structure%positions => dyn%posIon
   end function make_structure

   !> Compute the vdw interaction with DFT-D4
   type(vdw_interaction) function vdw_forces_d4(settings, structure) result (results)
#ifdef DFTD4
      use dftd4, only: structure_type, d4_model, realspace_cutoff, new_d4_model, get_dispersion
#else
      use tutor, only: vtutor
#endif
      type(vdw_settings), intent(in) :: settings
      type(vdw_structure), intent(in) :: structure
#ifdef DFTD4
      type(structure_type) struct
      type(d4_model) model
      !
      allocate(results%forces(3, size(structure%types)))
      struct = convert_structure_to_dftd4(structure)
      call new_d4_model(model, struct)
      call get_dispersion(struct, model, settings%parameters, realspace_cutoff(), &
         results%energy, results%forces, results%stress)
      call convert_units_of_results(results)
#else
      call vtutor%error("The code was not compiled with DFTD4 support.")
#endif
   end function vdw_forces_d4

#ifdef DFTD4

   !> Convert VASP data to the format used by DFTD4
   type (structure_type) function convert_structure_to_dftd4(structure) result (struct)
      use constant, only: AuToA
      use dftd4, only: structure_type, new
      use tutor, only: vtutor
      type(vdw_structure), intent(in) :: structure
      integer, allocatable :: atomic_numbers(:)
      real(q), allocatable :: positions(:,:)
      real(q) lattice(3,3)
      integer num_atom
      !
      num_atom = size(structure%types)
      allocate(atomic_numbers(num_atom), positions(3, num_atom))
      atomic_numbers = atomic_number(structure%elements(structure%types))
      if (any(atomic_numbers == unknown_atom)) &
         call vtutor%bug("Could not convert the elements to atomic numbers", __FILE__, __LINE__)
      lattice = structure%cell / AuToA
      positions = matmul(lattice, structure%positions)
      call new(struct, atomic_numbers, positions, lattice=lattice)
   end function convert_structure_to_dftd4

   !> Map an ion type onto the atomic number
   elemental integer function atomic_number(element)
      use string, only: lowercase
      character(len=2), intent(in) :: element
      character(len=2), parameter :: elements(94) = ['h ', 'he', &
         'li', 'be', 'b ', 'c ', 'n ', 'o ', 'f ', 'ne', &
         'na', 'mg', 'al', 'si', 'p ', 's ', 'cl', 'ar', &
         'k ', 'ca', 'sc', 'ti', 'v ', 'cr', 'mn', 'fe', 'co', 'ni', 'cu', &
            'zn', 'ga', 'ge', 'as', 'se', 'br', 'kr', &
         'rb', 'sr', 'y ', 'zr', 'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag', &
            'cd', 'in', 'sn', 'sb', 'te', 'i ', 'xe', &
         'cs', 'ba', 'la', 'ce', 'pr', 'nd', 'pm', 'sm', 'eu', 'gd', 'tb', 'dy', &
            'ho', 'er', 'tm', 'yb', 'lu', 'hf', 'ta', 'w ', 're', 'os', 'ir', 'pt', &
            'au', 'hg', 'tl', 'pb', 'bi', 'po', 'at', 'rn', &
         'fr', 'ra', 'ac', 'th', 'pa', 'u ', 'np', 'pu']
      character(len=2) element_
      integer i
      element_ = lowercase(element)
      atomic_number = unknown_atom
      do i = 1, size(elements)
         if (element_ == elements(i)) then
            atomic_number = i
            exit
         end if
      end do
   end function atomic_number

   !> Convert between atomic units used by DFTD4 and Angstroem-eV units used by VASP.
   subroutine convert_units_of_results(results)
      use constant, only: RyToEv, AuToA
      type(vdw_interaction), intent(inout) :: results
      real(q), parameter :: HtrToEv = 2 * RyToEv
      results%energy = results%energy * HtrToEv
      results%forces = -results%forces * HtrToEv / AuToA
      results%stress = -results%stress * HtrToEv
   end subroutine convert_units_of_results
#endif

   !> Read all settings in the INCAR file related to the DFTD4 library
   subroutine vdw_reader(io, xc_functional, settings)
      use base, only: in_struct
      use reader_tags, only: writexmlincar, process_incar, open_incar_if_found, close_incar_if_found
      use tutor, only: vtutor
      type(in_struct), intent(in) :: io
      character(len=*), intent(in) :: xc_functional
      type(vdw_settings), intent(out) :: settings
      integer ierr, ivdw
      logical lopen
      !
      ivdw = 0
      call open_incar_if_found(io%iu5, lopen)
      call process_incar(lopen, io%iu0, io%iu5, 'IVDW', ivdw, ierr, writexmlincar)
      if (ivdw == 13) call read_dftd4_tags(settings)
      call close_incar_if_found(io%iu5)
      settings%ivdw = ivdw
      settings%io = io
      !
   contains
      !
      subroutine read_dftd4_tags(settings)
         type(vdw_settings), intent(out) :: settings
#ifdef DFTD4
         character(len=:), allocatable :: user_defined
         real(q) rdum
         logical unknown
         !
         settings = default_settings(xc_functional, ierr)
         unknown = ierr == unknown_functional
         user_defined = ''
         !
         associate (param => settings%parameters)
            call process_incar(lopen, io%iu0, io%iu5, 'VDW_S6', rdum, ierr, writexmlincar)
            if (ierr == 0) call overwrite_if_sane('VDW_S6', [0.0_q, 2.0_q], rdum, param%s6, user_defined)
            call process_incar(lopen, io%iu6, io%iu5, 'VDW_S8', rdum, ierr, writexmlincar)
            if (ierr == 0) call overwrite_if_sane('VDW_S8', [0.0_q, 10.0_q], rdum, param%s8, user_defined)
            call process_incar(lopen, io%iu6, io%iu5, 'VDW_A1', rdum, ierr, writexmlincar)
            if (ierr == 0) call overwrite_if_sane('VDW_A1', [0.0_q, 10.0_q], rdum, param%a1, user_defined)
            call process_incar(lopen, io%iu6, io%iu5, 'VDW_A2', rdum, ierr, writexmlincar)
            if (ierr == 0) call overwrite_if_sane('VDW_A2', [0.0_q, 10.0_q], rdum, param%a2, user_defined)
         end associate
         !
         if (unknown) then
            if (len_trim(user_defined) > 0) then
               call vtutor%warning("Using user supplied parameters" // user_defined // &
                  "\nbecause the xc functional is not supported. Required parameters not &
                  &defined in the INCAR file will be replaced with PBE defaults ones.")
            else
               call vtutor%error('Unsupported xc functional, please specify VDW_A1, &
                  &VDW_A2, VDW_S6, and VDW_S8 in the INCAR file.')
            end if
         else if (len_trim(user_defined) > 0) then
            call vtutor%alert("Specifying the vdW parameters" // user_defined // &
               "\nin the INCAR file will overwrite the defaults for " // trim(xc_functional) // &
               ". Please make sure that is what you intended.")
         end if
#else
         call vtutor%error("The code was not compiled with DFTD4 support. Please &
            recompile the code adding -DDFTD4 to the CPP_OPTIONS and add the &
            necessary include and link parameters.")
#endif
      end subroutine read_dftd4_tags
      !
      subroutine overwrite_if_sane(tag, sane_interval, new_value, use_value, user_defined)
         use string, only: str
         character(len=*), intent(in) :: tag
         real(q), intent(in) :: sane_interval(2), new_value
         real(q), intent(inout) :: use_value
         character(len=:), allocatable, intent(inout) :: user_defined
         if ((minval(sane_interval) <= new_value).and.(new_value <= maxval(sane_interval))) then
            user_defined = user_defined // '\n' // tag // ' = ' // str(new_value)
            use_value = new_value
         else
            call vtutor%alert(tag // ' is not reasonable. Taking default.')
         end if
      end subroutine overwrite_if_sane
      !
   end subroutine vdw_reader

   subroutine vdw_writer(settings)
       type(vdw_settings), intent(in) :: settings
       character(len=*), parameter :: fmt_int = '(3x,a,t13,i6,4x,a)'
       integer unit_
       unit_ = settings%io%iu6
       if (unit_ < 0) &
           return
       write(unit_,*) 'Van der Waals corrections'
       write(unit_,fmt_int) 'IVDW    =', settings%ivdw, 'specifies the selected vdW correction'
       if (settings%ivdw == 13) then
#ifdef DFTD4
          call write_dftd4_info(unit_, settings%parameters)
#else
          write(unit_,'(3x,a)') 'DFT-D4  = not available'
#endif
       end if
       write(unit_,'(a)')
   end subroutine vdw_writer

#ifdef DFTD4
   subroutine write_dftd4_info(unit_, parameters)
       use dftd4, only: get_dftd4_version
       integer, intent(in) :: unit_
       type(rational_damping_param), intent(in) :: parameters
       character(len=:), allocatable :: version
       character(len=*), parameter :: fmt_real = '(3x,a,t13,f8.4,2x,a)'
       call get_dftd4_version(string=version)
       write(unit_,'(3x,a)') 'DFT-D4  = v' // version
       write(unit_,fmt_real) 'VDW_S6  =', parameters%s6, 'S6 parameter of DFT-D4'
       write(unit_,fmt_real) 'VDW_S8  =', parameters%s8, 'S8 parameter of DFT-D4'
       write(unit_,fmt_real) 'VDW_S9  =', parameters%s9, 'S9 parameter of DFT-D4'
       write(unit_,fmt_real) 'VDW_A1  =', parameters%a1, 'A1 parameter of DFT-D4'
       write(unit_,fmt_real) 'VDW_A2  =', parameters%a2, 'A2 parameter of DFT-D4'
       write(unit_,fmt_real) 'VDW_ALP =', parameters%alp, 'ALP parameter of DFT-D4'
   end subroutine write_dftd4_info
#endif


   type(vdw_settings) function default_settings(xc_functional, ierr) result (settings)
#ifdef DFTD4
      use dftd4, only: get_rational_damping, damping_param
      use tutor, only: vtutor
      class(damping_param), allocatable :: param
#endif
      character(len=*), intent(in) :: xc_functional
      integer, intent(out) :: ierr
      ierr = 0
#ifdef DFTD4
      call get_rational_damping(xc_functional, param)
      if (.not.allocated(param)) then
          ierr = unknown_functional
          call get_rational_damping('pbe', param)
      end if
      select type (param)
      type is (rational_damping_param)
          settings%parameters = param
      class default
          call vtutor%bug("These damping parameters are not implemented", __FILE__, __LINE__)
      end select
#endif
   end function default_settings

end module vdwd4
