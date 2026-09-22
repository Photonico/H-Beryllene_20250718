# 1 "libmbd.F"
!> Wrapper around the libMBD library
!>
!> Translates the VASP structural data into a format compatible with the libMBD
!> library. libMBD is an external library and needs to be installed separately.
!> Please check the [documentation](https://github.com/libmbd/libmbd) for more
!> information.
!>
!> The RPA-related options are disabled (see commented lines below) since
!> it is not clear how to use them.
module libmbd

   use prec, only: q
# 15

   implicit none

   type libmbd_structure
      integer, pointer :: nions => null()
      real(q), pointer :: cell(:,:) => null()
      real(q), pointer :: cell_recip(:,:) => null()
      integer, pointer :: types(:) => null()
      character(len=2), pointer :: elements(:) => null()
      real(q), pointer :: positions(:,:) => null()
   end type libmbd_structure

   type libmbd_interaction
      real(q) energy, stress(3,3)
      real(q), allocatable :: forces(:,:)
!      real(q), allocatable :: spectrum(:), modes(:,:)
   end type libmbd_interaction

   contains

   type(libmbd_structure) function structure_settings_libmbd(latt_cur,dyn,t_info,elem) result (structure)
      use lattice, only: latt
      use poscar, only: dynamics, type_info
      type(latt), target, intent(in) :: latt_cur
      type(dynamics), target, intent(in) :: dyn
      type(type_info), target, intent(in) :: t_info
      character(len=2), target, intent(in) :: elem(:)
      structure%nions => t_info%nions
      structure%cell => latt_cur%a
      structure%cell_recip => latt_cur%b
      structure%types => t_info%iTyp
      structure%elements => elem
      structure%positions => dyn%posIon
   end function structure_settings_libmbd

# 276


# 352


end module libmbd
