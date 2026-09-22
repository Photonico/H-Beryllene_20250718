#include "symbol.inc"
#ifdef libvaspml
module mvaspml

   use iso_c_binding, only: c_ptr, c_null_ptr, c_associated
   use ml_ff_prec
   use tutor, only: vtutor
   use mpi_help, only : ml_mpi_par
   use verlet_cell_list_mod, only : verlet_cell_list_type
   implicit none
   private
   public :: vaspml_wrapper

   !================================================================================================
   !> Wrapper around vaspml Frame class.
   !>
   !> This derived type stores a C pointer to an instance of a vasmpl Frame class. The instance is
   !> created when the init() member function is called and destroyed with finalize(). Member
   !> functions are used to pass data between VASP and vaspml.
   !================================================================================================
   type vaspml_wrapper

      private
      !> C-pointer to instance of vaspml Frame class.
      type(c_ptr) :: frame = c_null_ptr
      type( verlet_cell_list_type ) :: nearest_neighbor_radial
      type( verlet_cell_list_type ) :: nearest_neighbor_angular
      real( kind = q ),allocatable  :: posion( : , : )
      real( kind = q )              :: lattice( 1:3 , 1:3 )
      real( kind = q )              :: inverse_lattice( 1:3 , 1:3 )
      integer                       :: nions
      integer                       :: ntype
      integer,allocatable           :: itype( : )
      type( ml_mpi_par )            :: communicator
      real( kind = q )              :: w1
      real( kind = q )              :: w2
      real( kind = q )              :: rcut1
      real( kind = q )              :: rcut2
      integer                       :: number_atoms_per_node
      integer,allocatable           :: number_atoms_per_type_per_node( : )
      character( len=: ),allocatable      :: types
      real( kind = q )              :: omega
      

      contains
      procedure, public, pass( this ) :: init => wrapper_init
      procedure, public, pass( this ) :: compute => wrapper_main
      procedure, public, pass( this ) :: finalize => wrapper_finalize
      procedure, public, pass( this ) :: fill_neighbor_list_arrays => fill_neighbor_list_arrays
      procedure, public, pass( this ) :: rescale_atoms => rescale_atoms
      procedure, public, pass( this ) :: define_elements_per_thread => define_elements_per_thread
      procedure, public, pass( this ) :: prepare_types => prepare_types 
      procedure, public, pass( this ) :: fill_force_array => fill_force_array
      procedure, public, pass( this ) :: get_energy => get_energy
      procedure, public, pass( this ) :: get_stress_tensor => get_stress_tensor
   end type vaspml_wrapper

   interface
      !=============================================================================================
      !> Interface to C function which creates Frame instance.
      !=============================================================================================
      function create_frame() result(ptr) bind(c, name="createFrame")
         use iso_c_binding
         type(c_ptr) :: ptr
      end function create_frame
      !=============================================================================================
      !> Interface to C function which reads in the force field data.
      !=============================================================================================
      subroutine setup_force_field(ptr,mpi_comm) bind(c, name="setupForceField")
         use iso_c_binding
         type(c_ptr), intent(inout)     :: ptr
         integer( c_int ),optional, intent( in ) :: mpi_comm
      end subroutine setup_force_field
      !=============================================================================================
      !> Interface to C function which returns ML_W1.
      !=============================================================================================
      function get_w1(ptr) bind( c, name="get_W1") result( w1 )
         use iso_c_binding
         use ml_ff_prec
         type( c_ptr ), intent( inout ) :: ptr
         real( kind = c_double )        :: w1
      end function get_w1
      !=============================================================================================
      !> Interface to C function which returns ML_W2.
      !=============================================================================================
      function get_w2(ptr) bind( c, name="get_W2" ) result( w2 )
         use iso_c_binding
         use ml_ff_prec
         type( c_ptr ), intent( inout ) :: ptr
         real( kind = c_double )        :: w2
      end function get_w2
      !=============================================================================================
      !> Interface to C function which returns ML_RCUT1.
      !=============================================================================================
      function get_rcut1( ptr ) bind( c, name="get_RCUT1" ) result( rcut1 )
         use iso_c_binding
         use ml_ff_prec
         type( c_ptr ), intent( inout ) :: ptr
         real( kind = c_double )        :: rcut1
      end function get_rcut1
      !=============================================================================================
      !> Interface to C function which returns ML_RCUT2.
      !=============================================================================================
      function get_rcut2( ptr ) bind( c, name="get_RCUT2") result( rcut2 )
         use iso_c_binding
         use ml_ff_prec
         type( c_ptr ), intent( inout ) :: ptr
         real( kind = c_double )        :: rcut2
      end function get_rcut2
      !=============================================================================================
      !> Interface to C function which assigns number of atoms per core to c++ neighbor list.
      !=============================================================================================
      subroutine resize_neighbor_arrays( ptr, nions, key ) bind( c, name="resizeNeighborArrays")
         use iso_c_binding
         type(c_ptr), intent( inout ) :: ptr
         integer( kind = c_int )      :: nions
         character(kind=c_char), dimension(*), intent(in) :: key
      end subroutine resize_neighbor_arrays
      !=============================================================================================
      !> Interface to C function which destroys Frame instance.
      !=============================================================================================
      subroutine destroy_frame(ptr) bind(c, name="destroyFrame")
         use iso_c_binding
         type(c_ptr), intent(inout) :: ptr
      end subroutine destroy_frame
      !=============================================================================================
      !> Interface to C function which sets the number of atoms per type per node for the 
      !> c++ neighbor list
      !> 
      !> ntypes number of types determines length of number_atoms_per_type_per_node
      !> number_atoms_per_type_per_node stores the number of elements for certain type for 
      !> a single thread
      !=============================================================================================
      subroutine set_types_per_atom( ptr, ntypes, number_atoms_per_type_per_node, key ) bind(c, name="set_nAtomsType")
         use iso_c_binding
         type(c_ptr), intent(inout) :: ptr
         integer( kind = c_int ), intent( in ) :: ntypes
         integer( kind = c_int ), intent( in ) :: number_atoms_per_type_per_node( 1:ntypes )
         character(kind=c_char), dimension(*), intent(in) :: key
      end subroutine set_types_per_atom
      subroutine fill_neighbor_list_single_atom( ptr,&
                                                 number_neighbors,&
                                                 number_element, &
                                                 central_type,&
                                                 neighbors,&
                                                 neighbors_types,&
                                                 neighbor_dist,&
                                                 neighbor_connect,&
                                                 key ) bind( c, name = "fillNeighhborArrays" )
         use iso_c_binding
         type(c_ptr), intent(inout) :: ptr
         integer( kind = c_int ),intent( in ) :: number_neighbors
         integer( kind = c_int ),intent( in ) :: number_element
         integer( kind = c_int ),intent( in ) :: central_type
         integer( kind = c_int ),intent( in ) :: neighbors( 1:number_neighbors )
         integer( kind = c_int ),intent( in ) :: neighbors_types( 1:number_neighbors )
         real( kind = c_double ),intent( in ) :: neighbor_dist( 1:number_neighbors )
         real( kind = c_double ),intent( in ) :: neighbor_connect( 1:3, number_neighbors )
         character( kind = c_char ), dimension(*), intent(in) :: key
      end subroutine fill_neighbor_list_single_atom
      !=============================================================================================
      !> update the frame instance such that forces are computed. These
      !> forces can then be sent to vasp for ionic update
      !=============================================================================================
      subroutine update_vaspml( ptr, omega ) bind( c, name = "update" )
         use iso_c_binding
         type(c_ptr), intent(inout) :: ptr
         real( kind = c_double ), intent( in ) :: omega
      end subroutine update_vaspml

      subroutine set_type_map( ptr, types ) bind( c, name = "set_typeMap" )
         use iso_c_binding
         type(c_ptr), intent(inout) :: ptr
         character( kind = c_char ), dimension(*), intent( in ) :: types
      end subroutine set_type_map
            
      subroutine fill_force_single_atom( ptr, nions, inions, number_element, tifor ) bind( c, name = "fillForceSingleAtom" )
         use iso_c_binding
         type(c_ptr), intent(inout) :: ptr
         integer( kind = c_int ), intent( in )   :: nions
         integer( kind = c_int ), intent( in )   :: inions
         integer( kind = c_int ), intent( in )   :: number_element
         real( kind = c_double ),intent( inout ) :: tifor( 1:3, 1:nions )
      end subroutine fill_force_single_atom

      subroutine get_potential_energy( ptr, toten ) bind( c, name = "getPotentialEnergy" )
         use iso_c_binding
         type(c_ptr), intent(inout) :: ptr
         real( kind = c_double ),intent( inout ) :: toten
      end subroutine get_potential_energy
      
      subroutine get_stress_tensor_vaspml( ptr, tsif ) bind( c, name = "getStressTensor" )
         use iso_c_binding
         type(c_ptr), intent(inout) :: ptr
         real( kind = c_double ),intent( inout ) :: tsif( 1:3, 1:3 )
      end subroutine get_stress_tensor_vaspml
   end interface

   contains
   !================================================================================================
   !> Initialize vaspml, create Frame instance and set up force field.
   !>
   !> communicator mpi world communicator from vasp
   !================================================================================================
   subroutine wrapper_init(this, communicator, posion, lattice, inverse_lattice, nions, ntype, itype )

      use ml_ff_constant, only : autoa, pi
      use ml_main_subroutines, only : initialize_cell_list_algorithm
      class(vaspml_wrapper), intent(inout) :: this
      type( ml_mpi_par ), intent( in )     :: communicator
      real( kind = q ), intent( in )       :: posion( :, : )
      real( kind = q ), intent( in )       :: lattice( :, : )
      real( kind = q ), intent( in )       :: inverse_lattice( :, : )
      integer, intent( in )                :: nions
      integer, intent( in )                :: ntype
      integer, intent( in )                :: itype( : )

      if (c_associated(this%frame)) then
         call vtutor%bug("vaspml_wrapper: Trying to create another Frame &
                         &instance although pointer is already associated.", &
                         __FILE__, __LINE__)
      end if
      this%frame = create_frame()
      if (.not. c_associated(this%frame)) then
         call vtutor%bug("vaspml_wrapper: Could not create Frame instance.", &
                         __FILE__, __LINE__)
      end if
      this%communicator         = communicator
      this%communicator%node_me = this%communicator%node_me - 1
      this%nions        = nions
      this%ntype        = ntype
      this%lattice      = lattice / autoa
      this%inverse_lattice = inverse_lattice * autoa * 2.0_q * pi
      call setup_force_field(this%frame,this%communicator%mpi_comm)

      if ( allocated( this%posion ) ) deallocate( this%posion )
      allocate( this%posion( 1:3, nions ) )
      if ( allocated( this%itype ) ) deallocate( this%itype )
      allocate( this%itype( 1:nions ) )
      this%itype = itype
      ! rescale positions to cartesian coordinates
      call this%rescale_atoms( posion, this%lattice )
      this%w1    =  get_w1( this%frame )
      this%w2    =  get_w2( this%frame )
      this%rcut1 =  get_rcut1( this%frame )
      this%rcut2 =  get_rcut2( this%frame )
      call initialize_cell_list_algorithm( this%nearest_neighbor_radial,&
                                           this%nearest_neighbor_angular,&
                                           this%rcut1,&
                                           this%w1,&
                                           this%rcut2,&
                                           this%w2,&
                                           this%nions,&
                                           this%lattice,&
                                           this%inverse_lattice,&
                                           this%itype,&
                                           this%posion,&
                                           this%ntype,&
                                           this%communicator%ncpu,&
                                           this%communicator%node_me )
   end subroutine wrapper_init
   !================================================================================================
   !> rescale 3xN atoms array by a 3x3 matrix
   !> can be used to convert atomic positions from direct to cartesian coordinates and vice versa
   !>
   !> posion atomic positions stored as (xyz, number ions )
   !> lattice defines the matrix by which the positions are multiplied
   !================================================================================================
   subroutine rescale_atoms( this, posion, lattice )

      class(vaspml_wrapper), intent(inout) :: this
      real( kind = q ), intent( in ) :: posion( :, : )
      real( kind = q ), intent( in ) :: lattice( :, : )

      integer :: inions
      
      do inions = 1, this%nions
         this%posion( 1:3, inions )= &
         posion( 1, inions ) * lattice( 1:3 , 1 ) + &
         posion( 2, inions ) * lattice( 1:3 , 2 ) + &
         posion( 3, inions ) * lattice( 1:3 , 3 )
      enddo

   end subroutine rescale_atoms

   subroutine wrapper_main( this, posion, lattice, inverse_lattice, nions, ntype, itype, types, omega, toten, tifor, tsif )

      use iso_c_binding
      use ml_ff_constant, only : autoa, pi
      use force_field, only : compute_nearest_neighbor_arrays
      class(vaspml_wrapper), intent(inout) :: this
      real( kind = q ), intent( in )       :: posion( :, : )
      real( kind = q ), intent( in )       :: lattice( :, : )
      real( kind = q ), intent( in )       :: inverse_lattice( :, : )
      integer, intent( in )                :: nions
      integer, intent( in )                :: ntype
      integer, intent( in )                :: itype( : )
      character(len=2),intent(in)          :: types( : )
      real( kind = q ), intent( in )       :: omega
      real( kind = q ), intent( out )      :: toten
      real( kind = q ), intent( out )      :: tifor( :, : )
      real( kind = q ), intent( out )      :: tsif( :, : )

      integer :: inions

      integer :: mneib_est1_dummy
      integer :: mneib_est2_dummy

      this%nions        = nions
      this%ntype        = ntype
      this%lattice      = lattice / autoa
      this%inverse_lattice = inverse_lattice * autoa * 2.0_q * pi
      this%omega        = omega / autoa**3
      if ( allocated( this%posion ) ) deallocate( this%posion )
      allocate( this%posion( 1:3, nions ) )
      if ( allocated( this%itype ) ) deallocate( this%itype )
      allocate( this%itype( 1:nions ) )
      this%nions        = nions
      this%itype = itype
      call this%prepare_types( types )
      call set_type_map( this%frame, this%types )
      ! rescale positions to cartesian coordinates
      call this%rescale_atoms( posion, this%lattice )

      call compute_nearest_neighbor_arrays( this%nearest_neighbor_radial, &
                                            this%nearest_neighbor_angular, &
                                            this%lattice, &
                                            this%inverse_lattice, &
                                            this%itype, &
                                            this%ntype, &
                                            this%nions, &
                                            this%nions, &
                                            this%posion, &
                                            this%communicator,&
                                            ! radial variables
                                            mneib_est1_dummy, &
                                            this%w1, &
                                            this%rcut1, &
                                            ! angular variables
                                            mneib_est2_dummy, &
                                            this%w2, &
                                            this%rcut2 )
      !! fill neighbor arrays
      call this%fill_neighbor_list_arrays
      !!! call interface update
      call update_vaspml( this%frame, this%omega )
      !!! fill toten
      call this%get_energy( toten )
      !!! fill tifor
      call this%fill_force_array( tifor )
      !!! fill tisif
      call this%get_stress_tensor( tsif )

   end subroutine wrapper_main


   subroutine prepare_types( this, types )
      
      use iso_c_binding
      class(vaspml_wrapper), intent(inout) :: this
      character(len=2),intent(in)          :: types( : )

      integer :: i

      this%types = ""
      do i = 1, this%ntype
         this%types =  this%types // types( i ) //";"
      end do
      this%types = this%types // C_NULL_CHAR
   end subroutine prepare_types

   !================================================================================================
   !> Finalize vaspml, destroy Frame instance.
   !================================================================================================
   subroutine fill_neighbor_list_arrays( this )
      use iso_c_binding
      class(vaspml_wrapper), intent(inout) :: this
      integer :: inions
      integer :: inions_tmp
      integer :: number_element
      character( len =  7 ) :: key

      call this%define_elements_per_thread
      if ( this%w1 > 0 ) then
         call resize_neighbor_arrays( this%frame, this%number_atoms_per_node, "2-body"//C_NULL_CHAR )
         call set_types_per_atom( this%frame, this%ntype, this%number_atoms_per_type_per_node, "2-body"//C_NULL_CHAR )
      end if
      if ( this%w2 > 0 ) then
         call resize_neighbor_arrays( this%frame, this%number_atoms_per_node, "3-body"//C_NULL_CHAR )
         call set_types_per_atom( this%frame, this%ntype, this%number_atoms_per_type_per_node, "3-body"//C_NULL_CHAR )
      end if

      number_element = 0
      do inions_tmp = 1, this%nions, this%communicator%ncpu
         inions = inions_tmp + this%communicator%node_me
         if ( inions .le. this%nions ) then
            ! fill 2-body neighbor list
            if ( this%w1 > 0  ) then
               call fill_neighbor_list_single_atom( this%frame,&
                                                    size( this%nearest_neighbor_radial%nn_list( inions )%neighbors ),&
                                                    number_element,& 
                                                    ! counting starts at zero in c++ -> minus 1
                                                    this%itype( inions )-1,&
                                                    ! counting starts at zero in c++ -> minus 1
                                                    this%nearest_neighbor_radial%nn_list( inions )%neighbors-1,&
                                                    ! counting starts at zero in c++ -> minus 1
                                                    this%nearest_neighbor_radial%nn_list( inions )%neighbors_types-1,&
                                                    this%nearest_neighbor_radial%nn_list( inions )%neighbor_dist,&
                                                    -this%nearest_neighbor_radial%nn_list( inions )%neighbor_connect,&
                                                    "2-body"//C_NULL_CHAR )
            end if
            ! fill 3-body neighbor list
            if ( this%w2 > 0  ) then
               call fill_neighbor_list_single_atom( this%frame,&
                                                    size( this%nearest_neighbor_angular%nn_list( inions )%neighbors ),&
                                                    number_element,& 
                                                    ! counting starts at zero in c++ -> minus 1
                                                    this%itype( inions )-1,&
                                                    ! counting starts at zero in c++ -> minus 1
                                                    this%nearest_neighbor_angular%nn_list( inions )%neighbors-1,&
                                                    ! counting starts at zero in c++ -> minus 1
                                                    this%nearest_neighbor_angular%nn_list( inions )%neighbors_types-1,&
                                                    this%nearest_neighbor_angular%nn_list( inions )%neighbor_dist,&
                                                    -this%nearest_neighbor_angular%nn_list( inions )%neighbor_connect,&
                                                    "3-body"//C_NULL_CHAR )
            end if
            number_element = number_element + 1
         end if
      end do
   end subroutine fill_neighbor_list_arrays
   
   
   subroutine fill_force_array( this, tifor )
   
      use ml_ff_constant, onlY : funit
      use mpi_help
      class(vaspml_wrapper), intent(inout) :: this
      real( kind = q ),intent( inout )     :: tifor( :, : )

      integer :: number_element
      integer :: inions_tmp
      integer :: inions
      integer :: ierror
      
      tifor = 0.0_q
      number_element = 0
      do inions_tmp = 1, this%nions, this%communicator%ncpu
         inions = inions_tmp + this%communicator%node_me
         if ( inions .le. this%nions ) then
            call fill_force_single_atom( this%frame, this%nions, inions-1, number_element, tifor( :, : ) )
            number_element = number_element + 1
         end if
      end do
      tifor( :, : ) = tifor( :, : ) * funit
      call MPI_allreduce( mpi_in_place, tifor, 3*this%nions, mpi_double, mpi_sum, this%communicator%mpi_comm, ierror )
      
   end subroutine fill_force_array


   !================================================================================================
   !> Finalize vaspml, destroy Frame instance.
   !================================================================================================
   subroutine wrapper_finalize(this)

      class(vaspml_wrapper), intent(inout) :: this

      if (.not. c_associated(this%frame)) then
         call vtutor%bug("vaspml_wrapper: Trying to destroy Frame instance &
                         &although pointer is not yet associated.", &
                         __FILE__, __LINE__)
      end if
      call destroy_frame(this%frame)
      this%frame = c_null_ptr

   end subroutine wrapper_finalize
   !================================================================================================
   !> determine how many atoms each thread on the c++ side will get. Also determine how many
   !> atoms per type per thread will be passed to the c++ code.
   !>
   !> computes 
   !> number_atoms_per_node storing the total number of atoms after round robin distribution
   !> is done
   !> number_atoms_per_type_per_node the total number of atoms per type after round robin 
   !> distribution of atoms.
   !================================================================================================
   subroutine define_elements_per_thread( this )
      class(vaspml_wrapper), intent(inout) :: this

      integer :: inions_tmp
      integer :: inions

      if ( allocated( this%number_atoms_per_type_per_node ) ) deallocate( this%number_atoms_per_type_per_node )
      allocate( this%number_atoms_per_type_per_node( 1:this%ntype ) )
      this%number_atoms_per_type_per_node = 0
      this%number_atoms_per_node = 0
      do inions_tmp = 1, this%nions, this%communicator%ncpu
         inions = inions_tmp + this%communicator%node_me
         if ( inions .le. this%nions ) then
            this%number_atoms_per_type_per_node( this%itype( inions ) ) = &
                          this%number_atoms_per_type_per_node( this%itype( inions ) ) + 1
            this%number_atoms_per_node = this%number_atoms_per_node + 1
         end if
      end do
   end subroutine define_elements_per_thread
   

   subroutine get_energy( this, toten )
      
      use ml_ff_constant, onlY : eunit
      use mpi_help
      class(vaspml_wrapper), intent(inout) :: this
      real( kind = q ),intent( inout )     :: toten
      
      integer :: ierror

      toten = 0.0_q
      call get_potential_energy( this%frame, toten )
      toten = toten * eunit 
      call MPI_allreduce( mpi_in_place, toten, 1, mpi_double, mpi_sum, this%communicator%mpi_comm, ierror )
   end subroutine get_energy
   
   subroutine get_stress_tensor( this, tsif )
      
      use ml_ff_constant, onlY : sunit, evtoj
      use mpi_help
      class(vaspml_wrapper), intent(inout) :: this
      real( kind = q ),intent( inout )     :: tsif( 1:3, 1:3 )
      
      integer :: ierror
      real( kind = q ) :: fac
      tsif = 0.0_q
      
      call get_stress_tensor_vaspml( this%frame, tsif )
      fac = evtoj*1d22/this%omega/autoa**3
      tsif = tsif * sunit / fac
      call MPI_allreduce( mpi_in_place, tsif, 9, mpi_double, mpi_sum, this%communicator%mpi_comm, ierror )

   end subroutine get_stress_tensor

end module mvaspml
#endif
