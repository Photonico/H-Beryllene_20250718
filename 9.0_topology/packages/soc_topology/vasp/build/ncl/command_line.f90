# 1 "command_line.F"
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


# 2 "command_line.F" 2 
!> Read arguments from the command line
!>
!> The point of command line arguments is that they facilitate getting some quicker
!> answers then reading the INCAR file would, because it is processed before the
!> INCAR file.
!>
!> Currently the following command line arguments are implemented
!>
!> # --version | -v
!>
!> Return some information about the version of VASP that is used.
!>
!> # --dry-run | -n
!>
!> halt program and wait to attach debugger
!>
!> # --debug | -g
!>
!> Execute only a part of the code to test whether the setup would fit into
!> memory and the input can be processed.
module command_line

    use base, only: in_struct
    use tutor, only: vtutor, isAlert, isAdvice, isError
    use build_info, only : cpp_options, link_line

    implicit none

    integer, parameter :: default_status = -1
    character(len=*), parameter :: dry_run = "VASP will execute in dry-run mode. &
        &This means that only a small subset of the code will be run through to &
        &get some insight whether the setup is reasonable. The input files will &
        &be read and some sanity checks are conducted on them. After addressing &
        &any potential issues rerun VASP removing the dry-run command line option."
# 42

    type parsed_argument
        character(len=:), allocatable :: description
        logical :: should_stop = .false.
        integer :: stop_status = 0
        integer :: status_ = default_status
    end type parsed_argument

contains

!> Parse any potential command line arguments and react to them.
    subroutine parse_command_line(io, should_write)
!> Some command line arguments modify the settings.
        type(in_struct), intent(inout) :: io
!> Switch writing to the standard output on or off. Useful to limit output to a single rank.
        logical, intent(in) :: should_write
        character(len=256) text
        type(parsed_argument) argument
        integer index_
        do index_ = 1, command_argument_count()
            call get_command_argument(index_, text)
            call parse_argument(trim(text), io, argument)
            if (should_write) call log_description(argument%description, argument%status_)
            if (argument%should_stop) call vtutor%stopCode(argument%stop_status)
        end do
    end subroutine parse_command_line

    pure subroutine parse_argument(text, io, argument)
        use version
        character(len=*), intent(in) :: text
        type(in_struct), intent(inout) :: io
        type(parsed_argument), intent(out) :: argument
        select case(text)
        case ('--version', '-v')
            argument%description = vasp()
            argument%should_stop = .true.
        case ('--cpp-options', '-c')
            argument%description = cpp_options
            argument%should_stop = .true.
        case ('--link-line', '-l')
            argument%description = link_line
            argument%should_stop = .true.
        case ('--dry-run', '-n')
            argument%description = dry_run
            argument%status_ = isAdvice
            io%dry_run = .true.
# 93

        case default
            argument%description = "Command line argument '" // trim(text) // "' was not understood."
            argument%should_stop = .true.
            argument%stop_status = isError
        end select
    end subroutine parse_argument

    subroutine log_description(description, status_)
        use iso_fortran_env, only: output_unit
        character(len=*), intent(in) :: description
        integer, intent(in) :: status_
        select case (status_)
        case (default_status)
            write(output_unit, '(a)') description
        case (isAlert)
            call vtutor%alert(description)
        case (isAdvice)
            call vtutor%advice(description)
        case (isError)
            call vtutor%error(description)
        end select
    end subroutine log_description

end module command_line
