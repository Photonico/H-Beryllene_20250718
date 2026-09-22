# 1 "subdftd4.F"
!> Wrapper around the DFTD4 library
!>
!> Translates the VASP structural data into a format compatible with the DFTD4
!> library. DFTD4 is an external library and needs to be installed separately.
!> Please check the [documentation](https://github.com/dftd4/dftd4) for more
!> information.
module vdwd4

   use prec, only: q
   use base, only: in_struct
# 13


   implicit none

   integer, parameter :: unknown_functional = 1, unknown_atom = -1, vdw_dftd4 = 13

# 21

   logical, parameter :: dftd4_available = .false.
   type rational_damping_param
   end type rational_damping_param


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
# 75

      use tutor, only: vtutor

      type(vdw_settings), intent(in) :: settings
      type(vdw_structure), intent(in) :: structure
# 90

      call vtutor%error("The code was not compiled with DFTD4 support.")

   end function vdw_forces_d4

# 154


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
# 213

         call vtutor%error("The code was not compiled with DFTD4 support. Please &
            recompile the code adding -DDFTD4 to the CPP_OPTIONS and add the &
            necessary include and link parameters.")

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
# 248

          write(unit_,'(3x,a)') 'DFT-D4  = not available'

       end if
       write(unit_,'(a)')
   end subroutine vdw_writer

# 271



   type(vdw_settings) function default_settings(xc_functional, ierr) result (settings)
# 279

      character(len=*), intent(in) :: xc_functional
      integer, intent(out) :: ierr
      ierr = 0
# 295

   end function default_settings

end module vdwd4
