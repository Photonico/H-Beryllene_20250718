# 1 "reader_base.F"
!> Backwards compatible wrapper around new incar reader
!>
!> Generally for new developments you should prefer to use the new routines
!> in the #incar_reader module. However, where this is not easily possible
!> or when you write code adjacent to old process_incar statement, you could
!> use the old style interface.
module reader_tags
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


# 9 "reader_base.F" 2 
  use prec

  use vhdf5_base, only: incar_found

  use incar_reader, only: incar_file, process_incar, count_elements, not_found_error

  implicit none

!> Keep the INCAR file in memory
  type(incar_file) incar_f

  logical, save :: writexmlincar

!> Interface between the old implementation and the new (1._q,0._q).
!>
!> This implementation uses several deprecated flags to read data from the INCAR file.
!> We still accept this flags for backwards compatiblity but they have no impact on
!> the result. The behavior is dictated by the underlying routine in #incar_reader.
  interface process_incar
     module procedure read_incar_logical_scalar
     module procedure read_incar_logical_array
     module procedure read_incar_integer_scalar
     module procedure read_incar_integer_array
     module procedure read_incar_double_scalar
     module procedure read_incar_double_array
     module procedure read_incar_double_array_2d
     module procedure read_incar_double_array_3d
     module procedure read_incar_string
  end interface process_incar

contains

!> deprecated routine for historic compatibility
  subroutine open_incar_if_found(iu5, lopen)
    integer :: iu5
    logical, optional :: lopen
  end subroutine

!> deprecated routine for historic compatibility
  subroutine close_incar_if_found(iu5)
    integer :: iu5
  end subroutine

  subroutine check_error_and_number(flag_name, iu0, ierr, n, numb, lcont)
    use string, only: str
    use tutor, only: vtutor
    character(len=*), intent(in)  :: flag_name
    integer, intent(in) :: iu0, n, numb
    integer, intent(inout) :: ierr
    logical, intent(inout) :: lcont

! add ierr=7 to denote n<numb
    if ((ierr==0).and.(n<numb)) ierr=7

! if an error was signaled and lcont=.TRUE. we return in error
    if ((ierr/=0).and.lcont) then
       lcont=.false. ; return
    endif

! stop if lcont=.FALSE. and the error is not "tag-not-present" (ierr/=3)
    if ((ierr/=0) .and.(ierr/=3)) then
       CALL vtutor%error("Error reading item " // flag_name // " from file INCAR. \n Error code &
          &was IERR= " // str(ierr) // " ... . Found N= " // str(n) // " data.")
    endif

! if the error is "tag-not-present" (ierr=3) we return in error
    if (ierr/=0) then
       lcont=.false. ; return
    endif

! if no error was signaled
    lcont=.true.

    return
  end subroutine check_error_and_number

!
! logical scalar
!
  subroutine read_incar_logical_scalar(lopen, iu0, iu5, flag_name, flag_value, ierr, lwritexml, lcontinue, foundnumber)
    logical, intent(inout) :: flag_value
# 1 "./reader_base.inc" 1 
    logical, intent(in) :: lopen ! ignored
    integer, intent(in) :: iu0, iu5 ! ignored
    character(len=*), intent(in)  :: flag_name
    logical, optional :: lwritexml ! ignored
    logical, optional :: lcontinue
    integer, optional :: foundnumber
    logical ignore_errors
    integer ierr, num_elements
# 15

!
    num_elements = count_elements(incar_f, flag_name)
    if (present(foundnumber)) foundnumber = num_elements
    if (num_elements == 0) then
        ierr = not_found_error
        return
    end if
!
    ierr = 0
    ignore_errors = .false.
    if (present(lcontinue)) ignore_errors = lcontinue
    if (ignore_errors) then
        call process_incar(incar_f, flag_name, flag_value, ierr)
    else
        call process_incar(incar_f, flag_name, flag_value)
    end if
# 91 "reader_base.F" 2 
  end subroutine read_incar_logical_scalar
!
! logical array 1d
!
  subroutine read_incar_logical_array(lopen, iu0, iu5, flag_name, flag_value_, numb, ierr, lwritexml, lcontinue, foundnumber)
    logical, dimension(:), target, contiguous, intent(inout) :: flag_value_


# 1 "./reader_base.inc" 1 
    logical, intent(in) :: lopen ! ignored
    integer, intent(in) :: iu0, iu5 ! ignored
    character(len=*), intent(in)  :: flag_name
    logical, optional :: lwritexml ! ignored
    logical, optional :: lcontinue
    integer, optional :: foundnumber
    logical ignore_errors
    integer ierr, num_elements

    integer, intent(in) :: numb
    logical, pointer :: flag_value(:)
    flag_value(1:numb) => flag_value_



!
    num_elements = count_elements(incar_f, flag_name)
    if (present(foundnumber)) foundnumber = num_elements
    if (num_elements == 0) then
        ierr = not_found_error
        return
    end if
!
    ierr = 0
    ignore_errors = .false.
    if (present(lcontinue)) ignore_errors = lcontinue
    if (ignore_errors) then
        call process_incar(incar_f, flag_name, flag_value, ierr)
    else
        call process_incar(incar_f, flag_name, flag_value)
    end if
# 100 "reader_base.F" 2 
  end subroutine read_incar_logical_array
!
! integer scalar
!
  subroutine read_incar_integer_scalar(lopen, iu0, iu5, flag_name, flag_value, ierr, lwritexml, lcontinue, foundnumber)
    integer, intent(inout) :: flag_value
# 1 "./reader_base.inc" 1 
    logical, intent(in) :: lopen ! ignored
    integer, intent(in) :: iu0, iu5 ! ignored
    character(len=*), intent(in)  :: flag_name
    logical, optional :: lwritexml ! ignored
    logical, optional :: lcontinue
    integer, optional :: foundnumber
    logical ignore_errors
    integer ierr, num_elements
# 15

!
    num_elements = count_elements(incar_f, flag_name)
    if (present(foundnumber)) foundnumber = num_elements
    if (num_elements == 0) then
        ierr = not_found_error
        return
    end if
!
    ierr = 0
    ignore_errors = .false.
    if (present(lcontinue)) ignore_errors = lcontinue
    if (ignore_errors) then
        call process_incar(incar_f, flag_name, flag_value, ierr)
    else
        call process_incar(incar_f, flag_name, flag_value)
    end if
# 107 "reader_base.F" 2 
  end subroutine read_incar_integer_scalar
!
! integer array 1d
!
  subroutine read_incar_integer_array(lopen, iu0, iu5, flag_name, flag_value_, numb, ierr, lwritexml, lcontinue, foundnumber)
    integer, dimension(:), target, contiguous, intent(inout) :: flag_value_


# 1 "./reader_base.inc" 1 
    logical, intent(in) :: lopen ! ignored
    integer, intent(in) :: iu0, iu5 ! ignored
    character(len=*), intent(in)  :: flag_name
    logical, optional :: lwritexml ! ignored
    logical, optional :: lcontinue
    integer, optional :: foundnumber
    logical ignore_errors
    integer ierr, num_elements

    integer, intent(in) :: numb
    integer, pointer :: flag_value(:)
    flag_value(1:numb) => flag_value_



!
    num_elements = count_elements(incar_f, flag_name)
    if (present(foundnumber)) foundnumber = num_elements
    if (num_elements == 0) then
        ierr = not_found_error
        return
    end if
!
    ierr = 0
    ignore_errors = .false.
    if (present(lcontinue)) ignore_errors = lcontinue
    if (ignore_errors) then
        call process_incar(incar_f, flag_name, flag_value, ierr)
    else
        call process_incar(incar_f, flag_name, flag_value)
    end if
# 116 "reader_base.F" 2 
  end subroutine read_incar_integer_array
!
! real scalar
!
  subroutine read_incar_double_scalar(lopen, iu0, iu5, flag_name, flag_value, ierr, lwritexml, lcontinue, foundnumber)
    real(q), intent(inout) :: flag_value
# 1 "./reader_base.inc" 1 
    logical, intent(in) :: lopen ! ignored
    integer, intent(in) :: iu0, iu5 ! ignored
    character(len=*), intent(in)  :: flag_name
    logical, optional :: lwritexml ! ignored
    logical, optional :: lcontinue
    integer, optional :: foundnumber
    logical ignore_errors
    integer ierr, num_elements
# 15

!
    num_elements = count_elements(incar_f, flag_name)
    if (present(foundnumber)) foundnumber = num_elements
    if (num_elements == 0) then
        ierr = not_found_error
        return
    end if
!
    ierr = 0
    ignore_errors = .false.
    if (present(lcontinue)) ignore_errors = lcontinue
    if (ignore_errors) then
        call process_incar(incar_f, flag_name, flag_value, ierr)
    else
        call process_incar(incar_f, flag_name, flag_value)
    end if
# 123 "reader_base.F" 2 
  end subroutine read_incar_double_scalar
!
! real array 1d
!
  subroutine read_incar_double_array(lopen, iu0, iu5, flag_name, flag_value_, numb, ierr, lwritexml, lcontinue, foundnumber)
    real(q), dimension(:), target, contiguous, intent(inout) :: flag_value_


# 1 "./reader_base.inc" 1 
    logical, intent(in) :: lopen ! ignored
    integer, intent(in) :: iu0, iu5 ! ignored
    character(len=*), intent(in)  :: flag_name
    logical, optional :: lwritexml ! ignored
    logical, optional :: lcontinue
    integer, optional :: foundnumber
    logical ignore_errors
    integer ierr, num_elements

    integer, intent(in) :: numb
    real(q), pointer :: flag_value(:)
    flag_value(1:numb) => flag_value_



!
    num_elements = count_elements(incar_f, flag_name)
    if (present(foundnumber)) foundnumber = num_elements
    if (num_elements == 0) then
        ierr = not_found_error
        return
    end if
!
    ierr = 0
    ignore_errors = .false.
    if (present(lcontinue)) ignore_errors = lcontinue
    if (ignore_errors) then
        call process_incar(incar_f, flag_name, flag_value, ierr)
    else
        call process_incar(incar_f, flag_name, flag_value)
    end if
# 132 "reader_base.F" 2 
  end subroutine read_incar_double_array

  subroutine read_incar_double_array_2d(lopen, iu0, iu5, flag_name, flag_value_, numb, ierr, lwritexml, lcontinue, foundnumber)
    real(q), dimension(:,:), target, contiguous, intent(inout) :: flag_value_


# 1 "./reader_base.inc" 1 
    logical, intent(in) :: lopen ! ignored
    integer, intent(in) :: iu0, iu5 ! ignored
    character(len=*), intent(in)  :: flag_name
    logical, optional :: lwritexml ! ignored
    logical, optional :: lcontinue
    integer, optional :: foundnumber
    logical ignore_errors
    integer ierr, num_elements

    integer, intent(in) :: numb
    real(q), pointer :: flag_value(:)
    flag_value(1:numb) => flag_value_



!
    num_elements = count_elements(incar_f, flag_name)
    if (present(foundnumber)) foundnumber = num_elements
    if (num_elements == 0) then
        ierr = not_found_error
        return
    end if
!
    ierr = 0
    ignore_errors = .false.
    if (present(lcontinue)) ignore_errors = lcontinue
    if (ignore_errors) then
        call process_incar(incar_f, flag_name, flag_value, ierr)
    else
        call process_incar(incar_f, flag_name, flag_value)
    end if
# 139 "reader_base.F" 2 
  end subroutine read_incar_double_array_2d

  subroutine read_incar_double_array_3d(lopen, iu0, iu5, flag_name, flag_value_, numb, ierr, lwritexml, lcontinue, foundnumber)
    real(q), dimension(:,:,:), target, contiguous, intent(inout) :: flag_value_


# 1 "./reader_base.inc" 1 
    logical, intent(in) :: lopen ! ignored
    integer, intent(in) :: iu0, iu5 ! ignored
    character(len=*), intent(in)  :: flag_name
    logical, optional :: lwritexml ! ignored
    logical, optional :: lcontinue
    integer, optional :: foundnumber
    logical ignore_errors
    integer ierr, num_elements

    integer, intent(in) :: numb
    real(q), pointer :: flag_value(:)
    flag_value(1:numb) => flag_value_



!
    num_elements = count_elements(incar_f, flag_name)
    if (present(foundnumber)) foundnumber = num_elements
    if (num_elements == 0) then
        ierr = not_found_error
        return
    end if
!
    ierr = 0
    ignore_errors = .false.
    if (present(lcontinue)) ignore_errors = lcontinue
    if (ignore_errors) then
        call process_incar(incar_f, flag_name, flag_value, ierr)
    else
        call process_incar(incar_f, flag_name, flag_value)
    end if
# 146 "reader_base.F" 2 
  end subroutine read_incar_double_array_3d

  subroutine read_incar_string(lopen, iu0, iu5, flag_name, flag_value_, numb, ierr, lwritexml, lcontinue, foundnumber)
    character(len=*), intent(inout) :: flag_value_
    character(len=:), allocatable :: flag_value
    integer, intent(in) :: numb
# 1 "./reader_base.inc" 1 
    logical, intent(in) :: lopen ! ignored
    integer, intent(in) :: iu0, iu5 ! ignored
    character(len=*), intent(in)  :: flag_name
    logical, optional :: lwritexml ! ignored
    logical, optional :: lcontinue
    integer, optional :: foundnumber
    logical ignore_errors
    integer ierr, num_elements
# 15

!
    num_elements = count_elements(incar_f, flag_name)
    if (present(foundnumber)) foundnumber = num_elements
    if (num_elements == 0) then
        ierr = not_found_error
        return
    end if
!
    ierr = 0
    ignore_errors = .false.
    if (present(lcontinue)) ignore_errors = lcontinue
    if (ignore_errors) then
        call process_incar(incar_f, flag_name, flag_value, ierr)
    else
        call process_incar(incar_f, flag_name, flag_value)
    end if
# 153 "reader_base.F" 2 
    if (allocated(flag_value)) flag_value_ = flag_value
  end subroutine read_incar_string

end module reader_tags
