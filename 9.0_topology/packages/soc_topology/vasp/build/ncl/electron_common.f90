# 1 "electron_common.F"
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


# 2 "electron_common.F" 2 

module electron_common

    use base, only: in_struct, info_struct, q
    use mpimy, only: communic

    implicit none

    private
    public checkAbort, testBreakCondition, testNumberOfStep

contains

    subroutine testBreakCondition(energyDiff, evalDiff, numStep, noMixing, info, abortWithoutConv)
        real(q), intent(in) :: energyDiff, evalDiff
        integer, intent(in) :: numStep
        logical, intent(in) :: noMixing
        type(info_struct), intent(inout) :: info
        logical, intent(out) :: abortWithoutConv
! eigenvalues and energy must be converged
        info%lAbort = max(abs(energyDiff), abs(evalDiff)) < info%eDiff
! charge-density not constant and in last cycle no change of charge
        if (.not.info%lMix .and. .not.info%lChCon .and. .not.noMixing) info%lAbort = .false.
        call testNumberOfStep(numStep, info%nElm, info, abortWithoutConv)
    end subroutine testBreakCondition

    subroutine testNumberOfStep(numStep, maxStep, info, abortWithoutConv)
        integer, intent(in) :: numStep, maxStep
        type(info_struct), intent(inout) :: info
        logical, intent(out) :: abortWithoutConv
        abortWithoutConv = .false.
! do not stop during the non-selfconsistent startup phase
        if (numStep <= abs(info%nElmDl)) info%lAbort = .false.
! do not stop before minimum number of iterations is reached
        if (numStep < min(abs(info%nElMin), maxStep)) info%lAbort = .false.
! but stop after INFO%NELM steps no matter where we are now
        if (numStep >= maxStep) then
! if abort was not set the code did not converged except if we
! performed a single shot calculation
             abortWithoutConv = .not.info%lAbort .and. (maxStep > 1)
             info%lAbort = .true.
        end if
    end subroutine testNumberOfStep

    subroutine checkAbort(abortWithoutConv, io, comm, info)
        logical, intent(in) :: abortWithoutConv
        type(in_struct), intent(in) :: io
        type(communic), intent(in) :: comm
        type(info_struct), intent(inout) :: info
        if (info%lAbort) then
            call writeAbortReason(abortWithoutConv, io, comm)
            return
        end if
        info%lSoft = readStopcar(io, comm)
        if (info%lSoft) call writeHardStop(io, comm)
    end subroutine checkAbort

    subroutine writeAbortReason(abortWithoutConv, io, comm)
        logical, intent(in) :: abortWithoutConv
        type(in_struct), intent(in) :: io
        type(communic), intent(in) :: comm
        integer node_me, ionode

        node_me = comm%node_me
        ionode = comm%ionode

        IF (NODE_ME==IONODE) THEN
        write(io%iu6,'(a)') new_line('n')
        if (.not.abortWithoutConv) then
            write(io%iu6, '(a)') '------------------------ aborting loop because&
                & EDIFF is reached ----------------------------------------'
        else
            write(io%iu6, '(a)') '------------------------ aborting loop EDIFF &
                &was not reached (unconverged)  ----------------------------'
        end if
        write(io%iu6,'(a)') new_line('n')
        ENDIF
    end subroutine writeAbortReason

    logical function readStopcar(io, comm) result (lSoft)
        type(in_struct), intent(in) :: io
        type(communic), intent(in) :: comm
        integer idum, itmp, ierr, ncount
        real(q) rdum
        complex(q) cdum
        character charac
        lSoft = .false.


        call RDATAB(io%lOpen, 'STOPCAR', 99, 'LABORT', '=', '#', ';', 'L', &
            idum, rdum, cdum, lSoft, charac, ncount, 1, ierr)
        itmp = 0
        if (lSoft) itmp = 1
        CALL M_sum_i(comm, itmp, 1)
        lSoft = itmp > 0


    end function readStopcar

    subroutine writeHardStop(io, comm)
        type(in_struct), intent(in) :: io
        type(communic), intent(in) :: comm
        integer node_me, ionode

        node_me = comm%node_me
        ionode = comm%ionode

        IF (NODE_ME==IONODE) THEN
        if (io%iu0 >= 0) write(io%iu0,*) 'hard stop encountered!  aborting job ...'
        write(io%iu6,'(a)') new_line('n')
        write(io%iu6,'(a)') '------------------------ aborting loop because hard &
            &stop was set ---------------------------------------'
        write(io%iu6,'(a)') new_line('n')
        ENDIF
    end subroutine writeHardStop

end module electron_common
