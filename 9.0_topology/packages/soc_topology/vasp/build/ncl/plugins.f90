# 1 "plugins.F"
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


# 2 "plugins.F" 2 
!> Interface VASP calculation with external plugins.
!>
!> This module contains the tools to convert data from VASP types to C structs
!> and calls a external C code that can then interface with other libraries or
!> codes.
!>
!> There are a few design choices here that came up during the development
!> process:
!>
!> * All preprocessor statements are limited to this file so that the rest of
!>   the code can remain without modification. There is a constant `use_plugins`
!>   defined in this module that can be used in place of the preprocessor
!>   statements in other module if it turns out to be necessary. To facilitate
!>   this design all public procedures will turn to empty placeholders if the
!>   code is compiled without the preprocessor flag.
!>
!> * The interface routines map closely to current functionality of VASP. For
!>   example, the local-potential update is called from within `POTLOK` and the
!>   force update from within `FORCE_AND_STRESS`. With this choice the plugins
!>   are less likely to break existing VASP functionality.
!>
!> * Return values of plugins are additions to the current values. This choice
!>   allows multiple plugins running on the same interface at the same time
!>   without conflict (unless the plugins themselves break this functionality).
!>
!> Adding a new interface requires the following steps:
!>
!> * Create a new function in this module to convert the VASP data to C structs.
!> * Define the C structs in both the C header and the Fortran file (folder plugins)
!> * Define the dataclasses in the Python package
!> * Please keep the interfaces in alphabetical order.
!>
!> All interface routines have a common logic to them:
!>
!> * They perform a sanity check.
!> * They allocate the necessary memory in a buffer.
!> * They encode the immutable VASP data in a object called constants.
!> * They encode the mutable VASP data in a object called additions.
!> * They call the interface to the plugin.
!> * They handle the error if it is nonzero.
!> * They decode the returned additions object and apply the changes to the VASP data.
!>
module plugins

    use base, only: q, plugins_settings, info_struct
    use iso_c_binding, only: c_int, c_ptr
    use mgrid_struct_def, only: grid_3d
    use mpimy, only: communic
    use poscar_struct_def, only: latt, type_info
    use pseudo_struct_def, only: potcar

# 55


    implicit none

!> Use this in code outside of this module instead of the preprocessor flags.
# 62

    logical, parameter, public :: use_plugins = .false.


! internal variables to handle errors when the module is initialized/finalized twice
    logical :: initialized = .false., finalized = .false.

!
! private helper types to contain the memory during the communication
!

    type, abstract :: buffer_base
        type(communic) comm
        logical broadcast, active_rank
    end type buffer_base

    type, extends(buffer_base) :: buffer_force_and_stress
        integer(c_int), allocatable :: atomic_numbers(:)
        real(q), allocatable :: ZVAL(:)
        real(q), allocatable :: charge_density(:,:,:)
        real(q), allocatable :: forces(:,:)
        real(q), allocatable :: stress(:,:)
    end type buffer_force_and_stress

    type, extends(buffer_base) :: buffer_local_potential
        integer(c_int), allocatable :: atomic_numbers(:)
        real(q), allocatable :: ZVAL(:)
        real(q), allocatable :: charge_density(:,:,:)
        real(q), allocatable :: hartree_potential(:,:,:)
        real(q), allocatable :: ion_potential(:,:,:)
        real(q), allocatable :: potential(:,:,:)
        real(q), allocatable :: dipole_moment(:)
    end type buffer_local_potential

    type, extends(buffer_base) :: buffer_structure
        real(q), allocatable :: lattice_vectors(:,:)
        real(q), allocatable :: positions(:,:)
        integer(c_int), allocatable :: atomic_numbers(:)
        real(q), allocatable :: charge_density(:,:,:)
    end type buffer_structure

    type, extends(buffer_base) :: buffer_occupancies
    end type buffer_occupancies

    private
    public start_plugins, stop_plugins, plugins_reader, write_plugins, &
        plugins_force_and_stress, plugins_local_potential, plugins_machine_learning, &
        plugins_occupancies, plugins_structure

contains

!> This routine must be called before any of the interface routines. You cannot call
!> this routine a second time to reinitialize after a prior `stop_plugins` call.
    subroutine start_plugins()
        use tutor, only:vtutor
        integer(c_int) :: error
        if (initialized) then
            call vtutor%bug("Plugins are already initialized, do not initialize them another time.", &
                "plugins.F", 120)
        end if
# 126

        initialized = .true.
    end subroutine start_plugins


!> This routine should be called at the end of the code. You must not call further
!> interface routines after calling this function.
    subroutine stop_plugins()
        use tutor, only: vtutor
        integer(c_int) :: error
        if (.not.initialized) then
            call vtutor%bug("Finalize called before initialize.", "plugins.F", 137)
        end if
        if (finalized) then
            call vtutor%bug("Plugins are already finalized, do not finalize them another time.", &
                "plugins.F", 141)
        end if
# 147

    end subroutine stop_plugins


!> Read the INCAR file and detect which interfaces should be activated.
!> You may use existing INCAR tags to activate plugins where this makes
!> sense like it is 1._q for IBRION.
    subroutine plugins_reader(incar, ibrion, settings)
        use incar_reader, only: incar_file, process_incar
        use string, only: lowercase
        use tutor, only: vtutor
        type(incar_file), intent(inout) :: incar
        integer, intent(inout) :: ibrion
        type(plugins_settings), intent(out) :: settings
        logical update_structure
!
        settings%mode = "--"
        call process_incar(incar, "PLUGINS/MODE", settings%mode)
        settings%mode = lowercase(settings%mode)
        call read_flag(incar, "PLUGINS/FORCE_AND_STRESS", settings%force_and_stress)
        call read_flag(incar, "PLUGINS/LOCAL_POTENTIAL", settings%local_potential)

        call read_flag(incar, "PLUGINS/MACHINE_LEARNING", settings%machine_learning)
        call process_incar(incar, "PLUGINS/ML_OUTBLOCK", settings%outblock)
        if (settings%outblock<1) then
           call vtutor%error("plugins_reader: &
                PLUGINS/OUTBLOCK must be larger than 0.")
        endif
! Set default of PLUGINS/ML_OUTPUT_MODE here
        if (settings%machine_learning) then
           settings%output_mode = 0
        else
           settings%output_mode = 1
        endif
        call process_incar(incar, "PLUGINS/ML_OUTPUT_MODE", settings%output_mode)

        call read_flag(incar, "PLUGINS/OCCUPANCIES", settings%occupancies)
        update_structure = .false.
        call read_flag(incar, "PLUGINS/STRUCTURE", update_structure)
        settings%structure = ibrion == 12 .or. update_structure
!
! sanity checks
        select case(settings%mode)
        case ("--")
            settings%mode = "serial"
        case ("serial", "parallel", "vasp")
            if (.not.use_plugins) &
                call error_no_plugin_support("PLUGINS/MODE")
        case default
            call vtutor%error("PLUGINS/MODE = " // settings%mode // " is not " // &
                "implemented. Please check the spelling and correct the selection.")
        end select
        if (.not.use_plugins.and.settings%structure) then
            call vtutor%error("IBRION = 12 requires PLUGINS. Please recompile " // &
                "the code with -DPLUGINS and rerun the calculation!")
        end if
        if (update_structure .and. (ibrion /= 0 .and. ibrion /= 12)) &
            call vtutor%error("IBRION must not be set if 'PLUGINS/STRUCTURE' is set")
        if (settings%structure) ibrion = 12
        if (.not.initialized) &
            call start_plugins()
    end subroutine plugins_reader

    subroutine read_flag(incar, key, flag)
        use incar_reader, only: incar_file, process_incar
        type(incar_file), intent(inout) :: incar
        character(len=*), intent(in) :: key
        logical, intent(inout) :: flag
!
        call process_incar(incar, key, flag)
        if (.not.use_plugins.and.flag) &
            call error_no_plugin_support(key)
    end subroutine read_flag

    subroutine error_no_plugin_support(key)
        use tutor, only: vtutor
        character(len=*), intent(in) :: key
        call vtutor%error(key // " was set in the INCAR file but the code " // &
            "was not compiled with PLUGINS support. Please recompile with " // &
            "-DPLUGINS and rerun the calculation!")
    end subroutine error_no_plugin_support

!> Report the settings to the OUTCAR file
    subroutine write_plugins(unit_, settings)
        integer, intent(in) :: unit_
        type(plugins_settings), intent(in) :: settings
        character(len=*), parameter :: format_ = "(3x,a22,l6,4x,a)"
        character(len=:), allocatable :: explanation
!
        write(unit_, *) "Python interface"
        if (use_plugins) then
            write(unit_, format_) "PLUGINS active      = ", .true., "the code was linked to PLUGINS"
            select case (settings%mode)
            case ("serial")
                explanation = "run plugins only on root process"
            case ("parallel")
                explanation = "run plugins in parallel but do not distribute data"
!            case ("vasp")
!                explanation = "run plugins in parallel with VASP's data distribution"
            end select
            write(unit_, '(3x,a22,a,t36,a)') "  -mode             = ", settings%mode, explanation
            write(unit_, format_) "  -force_and_stress = ", settings%force_and_stress, "update the forces and the stress externally"
            write(unit_, format_) "  -local_potential  = ", settings%local_potential, "update the potential externally"

            write(unit_, format_) "  -machine_learning = ", settings%machine_learning, "use machine_learning to change VASP flow"

            write(unit_, format_) "  -occupancies      = ", settings%occupancies, "update quantities related to the occupancies"
            write(unit_, format_) "  -structure        = ", settings%structure, "overwrite the structure externally"
        else
            write(unit_, format_) "PLUGINS active      = ", .false., "recompile with -DPLUGINS to activate"
        end if
    end subroutine write_plugins

!
! The interfaces to the plugins
!

!> This routine allows to modify the force and the stress before the ionic update
!>
!> Note: if the charge density is not present, ML interface is called instead
    subroutine plugins_force_and_stress(info, t_info, latt_cur, p, gridc, &
            energy, forces, stress, charge_density)
        type(info_struct), intent(in) :: info
        type(type_info), intent(in) :: t_info
        type(latt), intent(in) :: latt_cur
        type(potcar), intent(in) :: p(:)
        type(grid_3d), intent(in) :: gridc
        real(q), intent(inout) :: energy
        real(q), intent(inout) :: forces(:,:)
        real(q), intent(inout) :: stress(:,:)
        complex(q), intent(in), optional :: charge_density(:)
        integer(c_int) error
!
# 297

    end subroutine plugins_force_and_stress


!> This interface allows to add an external local potential in real space
    subroutine plugins_local_potential(info, t_info, latt_cur, p, gridc, charge_density, &
            potential, hartree_potential, ion_potential, dipole_moment, energy)
        type(info_struct), intent(in) :: info
        type(type_info), intent(in) :: t_info
        type(latt), intent(in) :: latt_cur
        type(potcar), intent(in) :: p(:)
        type(grid_3d), intent(in) :: gridc
        complex(q), intent(in) :: charge_density(:)
        complex(q), intent(inout) :: potential(:)
        complex(q), intent(in) :: hartree_potential(:)
        complex(q), intent(in) :: ion_potential(:)
        real(q), intent(in)    :: dipole_moment(:)
        real(q), intent(inout) :: energy
# 333

    end subroutine plugins_local_potential


!> This interface allows to provide an external machine learned force field
    subroutine plugins_machine_learning(info, t_info, latt_cur, p, gridc, dyn, &
        symm, isymop, nrotk, ldo_ab_initio, energy, forces, stress, press)
        use base, only: symmetry
        use constant, only: EVTOJ
        use poscar_struct_def, only: dynamics
!
        type(info_struct), intent(in) :: info
        type(type_info), intent(in) :: t_info
        type(latt), intent(in) :: latt_cur
        type(potcar), intent(in) :: p(:)
        type(grid_3d), intent(in) :: gridc
        type(dynamics), intent(in) :: dyn
        type(symmetry), intent(in) :: symm
        integer, intent(in) :: isymop(:,:,:), nrotk
        logical, intent(in) :: ldo_ab_initio
        real(q), intent(out) :: energy, forces(:,:), stress(:,:), press
!
# 374

!
    end subroutine plugins_machine_learning


!> This interface allow to change quantities related to the occupancies at the end of each ionic step
    subroutine plugins_occupancies(info, comm, nelect, efermi, efermi_set, sigma, ismear, emin, emax, nupdown)
        type(info_struct), intent(in) :: info
        type(communic), intent(in) :: comm
        real(q), intent(inout) :: nelect
        real(q), intent(in) :: efermi
        real(q), intent(inout) :: efermi_set
        real(q), intent(inout) :: sigma
        integer, intent(inout) :: ismear
        real(q), intent(inout) :: emin
        real(q), intent(inout) :: emax
        real(q), intent(inout) :: nupdown
        integer(c_int)         :: error
# 407

    end subroutine plugins_occupancies


!> This interface allows to update the structure at the end of an ionic step
!> replacing the _normal_ update procedure.
    subroutine plugins_structure(gridc, info, t_info, latt_cur, p, comm, energy, forces, stress, charge_density)
        type(grid_3d), intent(in) :: gridc
        type(info_struct), intent(in) :: info
        type(type_info), intent(inout) :: t_info
        type(latt), intent(inout) :: latt_cur
        type(potcar), intent(in) :: p(:)
        type(communic), intent(in) :: comm
        real(q), intent(in) :: energy
        real(q), intent(inout) :: forces(:,:)
        real(q), intent(inout) :: stress(:,:)
        complex(q), intent(in) :: charge_density(:)
        complex(q), allocatable :: cwork(:)
        integer(c_int) error
# 453

    end subroutine plugins_structure

!
! common helper routines used by all interfaces
!

    subroutine sanity_check_on_data_types_and_status
        use iso_c_binding, only: c_double
        use tutor, only: vtutor
        real(q) double_fortran
        real(c_double) double_cpp
        integer int_fortran
        integer(c_int) int_cpp
!
        if (.not.initialized) then
            call vtutor%bug("Please call `start_plugins` before calling any interface to a plugin", &
                "plugins.F", 470)
        end if
        if (finalized) then
            call vtutor%bug("Cannot call any interface after finalizing the plugins with `stop_plugins`", &
                "plugins.F", 474)
        end if
        if (kind(double_fortran) /= kind(double_cpp)) then
            call vtutor%bug("The bytesize of floating-point numbers is different. " &
                // "Please check the setup of the VASP datatypes!", "plugins.F", 478)
        end if
        if (kind(int_fortran) /= kind(int_cpp)) then
            call vtutor%bug("The bytesize of integer numbers is different. " &
                // "Please check the setup of the VASP datatypes!", "plugins.F", 482)
        end if
    end subroutine sanity_check_on_data_types_and_status

    subroutine error_handling(buffer, routine, error_this_rank)
        use string, only: str
        use tutor, only: vtutor
        class(buffer_base), intent(in) :: buffer
        character(len=*), intent(in) :: routine
        integer(c_int), intent(in) :: error_this_rank
        integer error_any_rank
!
# 523

    end subroutine error_handling


    subroutine merge_grid_quantity(mode, gridc, distributed, gathered)
        use tutor, only: vtutor
        character(len=*), intent(in) :: mode
        type(grid_3d), intent(in) :: gridc
        complex(q), intent(in) :: distributed(:)
        real(q), intent(out), allocatable :: gathered(:,:,:)
        COMPLEX(q), allocatable :: work(:)
        real(q), allocatable :: rwork(:,:)
        integer :: nz
!
        select case (mode)
        case ("serial", "parallel")
            allocate(work(gridc%mplwv))
            allocate(gathered(gridc%ngz, gridc%ngy, gridc%ngx))
            allocate(rwork(gridc%ngx,gridc%ngy))
            work = 0._q
            call rl_add(distributed, 1.0_q, distributed, 0.0_q, work, gridc)
            do nz = 1, gridc%ngz
               call mrg_grid_rl_plane(gridc, rwork, work, nz)
               gathered(nz,:,:) = transpose(rwork)
            end do
        case default
            call vtutor%bug("Merging grid quantities is not implemented for mode '" // mode // "'.", "plugins.F", 549)
        end select
    end subroutine merge_grid_quantity

    subroutine setup_buffer_base(buffer, mode, comm)
        class(buffer_base), intent(inout) :: buffer
        character(len=*), intent(in) :: mode
        type(communic), intent(in) :: comm
        buffer%comm = comm
        buffer%broadcast = mode == "serial"
        buffer%active_rank = mode /= "serial" .or. comm%node_me == 1
    end subroutine setup_buffer_base

# 885


! this routine should perhaps be moved to poscar.F

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
        integer, parameter :: unknown_atom = -1
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

end module
