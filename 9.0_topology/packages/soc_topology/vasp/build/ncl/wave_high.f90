# 1 "wave_high.F"
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


# 2 "wave_high.F" 2 


!***********************************************************************
!
!> this module implements high level routines to operate on wavefunctions
!>
!> care is taken to avoid indexing pointer arrays since
!> this incures performance penalties on all Intel compilers
!> (but particularly on the efc compiler)
!> (constructs such as W1%CPTWFP(I, ..) are avoided)
!> instead low level F77 routines are used, which are implemented
!> at the end of the routine
!
!***********************************************************************

MODULE wave_high
  USE prec
  USE wave
  USE wave_mpi
!***********************************************************************
!
!> interfaces for functions that accept only F77 style arrays
!> these F77 routines would cause a problem if the first element
!> of a pointer array is passed to the F77 routines
!
!***********************************************************************

  INTERFACE
     SUBROUTINE ARRAY_TO_W1( W1, C, CPROJ)
       USE wave
       TYPE (wavefun1)    W1
       COMPLEX(q):: C
       COMPLEX(q), OPTIONAL :: CPROJ
     END SUBROUTINE ARRAY_TO_W1
  END INTERFACE

  INTERFACE
     SUBROUTINE W1_TO_ARRAY( W1, C, CPROJ)
       USE wave
       TYPE (wavefun1)    W1
       COMPLEX(q):: C
       COMPLEX(q), OPTIONAL :: CPROJ
     END SUBROUTINE W1_TO_ARRAY
  END INTERFACE

  INTERFACE
     SUBROUTINE ECCP_NL(LMDIM,LMMAXC,CDIJ,CPROJ1,CPROJ2,CNL)
!$ACC ROUTINE VECTOR
       USE prec
       COMPLEX(q)      CNL
       INTEGER LMDIM, LMMAXC
       COMPLEX(q) CDIJ
       COMPLEX(q) CPROJ1,CPROJ2
     END SUBROUTINE ECCP_NL
  END INTERFACE

  INTERFACE
     SUBROUTINE OVERL(WDES1, LOVERL, LMDIM, CQIJ, CPROF, CRESUL)
       USE wave
       TYPE (wavedes1) WDES1
       LOGICAL LOVERL
       INTEGER LMDIM
       COMPLEX(q) CQIJ,CDIJ
       COMPLEX(q) CRESUL,CPROF
     END SUBROUTINE OVERL
  END INTERFACE

  INTERFACE
     SUBROUTINE OVERL1(WDES1, LMDIM, CDIJ, CQIJ, EVALUE, CPROF,CRESUL)
       USE wave
       TYPE (wavedes1) WDES1
       INTEGER LMDIM
       COMPLEX(q) CQIJ,CDIJ
       REAL(q) :: EVALUE
       COMPLEX(q) CRESUL,CPROF
     END SUBROUTINE OVERL1
  END INTERFACE

  INTERFACE
     SUBROUTINE OVERL1_C(WDES1, LMDIM, CDIJ, CQIJ, EVALUE, CPROF,CRESUL)
       USE wave
       TYPE (wavedes1) WDES1
       INTEGER LMDIM
       COMPLEX(q) CQIJ,CDIJ
       COMPLEX(q) EVALUE
       COMPLEX(q) CRESUL,CPROF
     END SUBROUTINE OVERL1_C
  END INTERFACE

  INTERFACE
     SUBROUTINE OVERL1_CCDIJ(WDES1, LMDIM, CDIJ, CQIJ, EVALUE, CPROF,CRESUL)
       USE wave
       TYPE (wavedes1) WDES1
       INTEGER LMDIM
       COMPLEX(q) CQIJ,CDIJ
       REAL(q) :: EVALUE
       COMPLEX(q) CRESUL,CPROF
     END SUBROUTINE OVERL1_CCDIJ
  END INTERFACE

  INTERFACE ELEMENT
     MODULE PROCEDURE W1_FROM_W
     MODULE PROCEDURE W1_FROM_WA
     MODULE PROCEDURE W1_FROM_WA2
  END INTERFACE

  INTERFACE ELEMENTS
     MODULE PROCEDURE WA_FROM_W
     MODULE PROCEDURE WA_FROM_WA
     MODULE PROCEDURE WA_FROM_WA2
  END INTERFACE

  INTERFACE REDISTRIBUTE_PROJ
     MODULE PROCEDURE W_REDIS_PROJ
     MODULE PROCEDURE WA_REDIS_PROJ
     MODULE PROCEDURE W1_REDIS_PROJ
  END INTERFACE

  INTERFACE REDISTRIBUTE_PW
     MODULE PROCEDURE W_REDIS_PW
     MODULE PROCEDURE WA_REDIS_PW
     MODULE PROCEDURE W1_REDIS_PW
  END INTERFACE

  INTERFACE ASSIGNMENT
     MODULE PROCEDURE W1_COPY_REVERSE_ARG
  END INTERFACE

  INTERFACE CONNECT
     MODULE PROCEDURE WA_CONNECT2_W1
     MODULE PROCEDURE WA_CONNECT2_AR
     MODULE PROCEDURE W1_CONNECT2_AR
  END INTERFACE

  INTERFACE DISCONNECT
     MODULE PROCEDURE WA_NULLIFY
     MODULE PROCEDURE W1_NULLIFY
  END INTERFACE

  CONTAINS

!***********************************************************************
!
!> FFT of a wavefunction to real space
!
!***********************************************************************

    SUBROUTINE FFTWAV_W1( W1)
      IMPLICIT NONE
      INTEGER ISPINOR
      TYPE (wavefun1)    W1


      
      IF (W1%WDES1%NK==0) THEN
         CALL vtutor%bug("internal error in FFTWAV_W1: NK is set to zero: WDES not set up properly", "wave_high.F", 157)
      ENDIF
      IF (W1%WDES1%LUSEINV) THEN
         IF (.NOT. ASSOCIATED(W1%WDES1%FFTSCA) .OR. .NOT. ASSOCIATED(W1%WDES1%NINDPW_INV)) THEN
            CALL vtutor%bug("internal error in FFTWAV_W1: FFTSCA is not associated and LUSEINV is &
               &set " // str(ASSOCIATED(W1%WDES1%FFTSCA)) // " " // str(ASSOCIATED(W1%WDES1%NINDPW_INV)), &
               "wave_high.F", 163)
         ENDIF
         DO ISPINOR=0,W1%WDES1%NRSPINORS-1
            CALL FFTWAV_USEINV(W1%WDES1%NGVECTOR, W1%WDES1%NINDPW(1), W1%WDES1%NINDPW_INV(1), W1%WDES1%FFTSCA(1,2), &
                 W1%CR(1+ISPINOR*W1%WDES1%GRID%MPLWV), &
                 W1%CPTWFP(1+ISPINOR*W1%WDES1%NGVECTOR),W1%WDES1%GRID)
         ENDDO
      ELSE
         DO ISPINOR=0,W1%WDES1%NRSPINORS-1
            CALL FFTWAV(W1%WDES1%NGVECTOR, W1%WDES1%NINDPW(1), &
                 W1%CR(1+ISPINOR*W1%WDES1%GRID%MPLWV), &
                 W1%CPTWFP(1+ISPINOR*W1%WDES1%NGVECTOR),W1%WDES1%GRID)
         ENDDO
      ENDIF

      
    END SUBROUTINE FFTWAV_W1

!
!> FFT of a wavefunction to real space
!> similar version as above, however, receives WDES1 and explicit "WORK" arrays CR and CPTWFP
!> this routine is new in vasp.6 and should consistently replace
!> FFTWAV
!

    SUBROUTINE FFTWAV_WDES1( WDES1, CR, CPTWFP)
      IMPLICIT NONE
!> descriptor for (1._q,0._q) k-point
      TYPE (wavedes1)    WDES1
!> input:  orbitals in real space
      COMPLEX(q):: CR(WDES1%GRID%RL%NP)
!> output: orbitals in reciprocal space
      COMPLEX(q):: CPTWFP(WDES1%NRPLWV)

      IF (WDES1%NK==0) THEN
         CALL vtutor%bug("internal error in FFTWAV_WDES1: NK is set to zero: WDES not set up properly", "wave_high.F", 198)
      ENDIF
      IF (WDES1%LUSEINV) THEN
         IF (.NOT. ASSOCIATED(WDES1%FFTSCA) .OR. .NOT. ASSOCIATED(WDES1%NINDPW_INV)) THEN
            CALL vtutor%bug("internal error in FFTWAV_WDES1: FFTSCA is not associated and LUSEINV &
               &is set " // str(ASSOCIATED(WDES1%FFTSCA)) // " " // str(ASSOCIATED(WDES1%NINDPW_INV)), &
               "wave_high.F", 204)
         ENDIF
         CALL FFTWAV_USEINV(WDES1%NGVECTOR, WDES1%NINDPW(1), WDES1%NINDPW_INV(1), WDES1%FFTSCA(1,2), &
              CR(1),CPTWFP(1),WDES1%GRID)
      ELSE
         CALL FFTWAV(WDES1%NGVECTOR, WDES1%NINDPW(1), &
              CR(1),CPTWFP(1),WDES1%GRID)
      ENDIF
    END SUBROUTINE FFTWAV_WDES1


!***********************************************************************
!
!> FFT of a wavefunction to real space
!> a WDES1 must be supplied
!> as well as the input and output array
!> this routine is new in vasp.6 and should consistently replace
!> FFTEXT
!
!***********************************************************************

    SUBROUTINE FFTEXT_WDES1(WDES1, CR, CPTWFP, LADD)
      IMPLICIT NONE
      INTEGER ISPINOR
!> descriptor for (1._q,0._q) k-point
      TYPE (wavedes1)    WDES1
!> input:  orbitals in real space
      COMPLEX(q):: CR(WDES1%GRID%RL%NP)
!> output: orbitals in reciprocal space
      COMPLEX(q):: CPTWFP(WDES1%NRPLWV)
!> add results to CPTWFP and does not clear CPTWFP beforehand
      LOGICAL    :: LADD

      IF (WDES1%NK==0) THEN
         CALL vtutor%bug("internal error in FFTWAV_WDES1: NK is set to zero: WDES not set up properly", "wave_high.F", 238)
      ENDIF
      IF (WDES1%LUSEINV) THEN
         IF (.NOT. ASSOCIATED(WDES1%FFTSCA) .OR. .NOT. ASSOCIATED(WDES1%NINDPW_INV)) THEN
            CALL vtutor%bug("internal error in FFTEXT_WDES1: FFTSCA is not associated and LUSEINV &
               &is set " // str(ASSOCIATED(WDES1%FFTSCA)) // " " // str(ASSOCIATED(WDES1%NINDPW_INV)), &
               "wave_high.F", 244)
         ENDIF
         CALL FFTEXT_USEINV(WDES1%NGVECTOR, WDES1%NINDPW(1), WDES1%FFTSCA(1,1), &
              CR(1), CPTWFP(1), WDES1%GRID, LADD)
      ELSE
         CALL FFTEXT(WDES1%NGVECTOR, WDES1%NINDPW(1), &
              CR(1), CPTWFP(1), WDES1%GRID, LADD)
      ENDIF
    END SUBROUTINE FFTEXT_WDES1


!***********************************************************************
!
!> This routine accomplishes the same thing as the assignment:
!> ~~~
!>   W1_LHS = W1_RHS
!> ~~~
!> would, but allows for the use of "fast"-memory (in the sense that
!> the fm_*arrays are excluded from the assignment).
!
!***********************************************************************

  SUBROUTINE W1_ASSIGN(W1_LHS,W1_RHS)
    IMPLICIT NONE
    TYPE (wavefun1), INTENT(IN) :: W1_RHS
    TYPE (wavefun1) :: W1_LHS

    W1_LHS%WDES1  =>W1_RHS%WDES1
    W1_LHS%FERWE  = W1_RHS%FERWE
    W1_LHS%AUX    = W1_RHS%AUX
    W1_LHS%CELEN  = W1_RHS%CELEN
    W1_LHS%NB     = W1_RHS%NB
    W1_LHS%ISP    = W1_RHS%ISP
    W1_LHS%LDO    = W1_RHS%LDO

    W1_LHS%CPTWFP =>W1_RHS%CPTWFP
    W1_LHS%CPROJ  =>W1_RHS%CPROJ
    W1_LHS%CR     =>W1_RHS%CR

  END SUBROUTINE W1_ASSIGN


!***********************************************************************
!
!> copy a W1 structure
!>  ~~~
!>  W2 = W1  (W1 -> W2)
!>  ~~~
!> the argument arrangement is similar to  DCOPY, ZCOPY in BLAS level 1
!> for syntactic sugar the reverse operation is also present
!
!> @details @ref openmp :
!> BLAS1 calls were replaced by explicit OMP PARALLEL loops
!> (for some reason BLAS1 call thread really badly).
!
!***********************************************************************

  SUBROUTINE W1_COPY( W1, W2)
    IMPLICIT NONE
    TYPE (wavefun1), INTENT(IN) :: W1
    TYPE (wavefun1) :: W2
    INTEGER I

    


# 312

    CALL ZCOPY( W1%WDES1%NRPLWV, W1%CPTWFP(1), 1, W2%CPTWFP(1), 1)


    IF (W1%WDES1%LGAMMA) THEN
       CALL DCOPY( W1%WDES1%NPROD, W1%CPROJ(1), 1, W2%CPROJ(1), 1)
    ELSE
       CALL ZCOPY( W1%WDES1%NPROD, W1%CPROJ(1), 1, W2%CPROJ(1), 1)
    ENDIF

    IF (ASSOCIATED(W1%CR) .AND. ASSOCIATED(W2%CR)) THEN
       IF (SIZE(W1%CR) /=W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS) THEN
          CALL vtutor%bug("W1_COPY: real space allocation is not correct","wave_high.F",324)
       ENDIF
       IF (SIZE(W1%CR) /=SIZE(W2%CR)) THEN
          CALL vtutor%bug("W1_COPY: real space allocation is different " // &
             str(SIZE(W1%CR)) // " " // str(SIZE(W2%CR)),"wave_high.F",328)
       ENDIF
       CALL ZCOPY( W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS, W1%CR(1), 1, W2%CR(1), 1)
    ENDIF
# 362

    

  END SUBROUTINE W1_COPY

  SUBROUTINE W1_COPY_NOCR( W1, W2)
    IMPLICIT NONE
    TYPE (wavefun1), INTENT(IN) :: W1
    TYPE (wavefun1) :: W2

# 374

    CALL ZCOPY( W1%WDES1%NRPLWV, W1%CPTWFP(1), 1, W2%CPTWFP(1), 1)


    IF (W1%WDES1%LGAMMA) THEN
       CALL DCOPY( W1%WDES1%NPROD, W1%CPROJ(1), 1,  W2%CPROJ(1), 1)
    ELSE
       CALL ZCOPY( W1%WDES1%NPROD, W1%CPROJ(1),  1, W2%CPROJ(1), 1)
    ENDIF
  END SUBROUTINE W1_COPY_NOCR

  SUBROUTINE W1_COPY_CPROJ( W1, W2)
    IMPLICIT NONE
    TYPE (wavefun1), INTENT(IN) :: W1
    TYPE (wavefun1) :: W2
    IF (W1%WDES1%LGAMMA) THEN
       CALL DCOPY( W1%WDES1%NPROD, W1%CPROJ(1), 1,  W2%CPROJ(1), 1)
    ELSE
       CALL ZCOPY( W1%WDES1%NPROD, W1%CPROJ(1),  1, W2%CPROJ(1), 1)
    ENDIF
  END SUBROUTINE W1_COPY_CPROJ

  SUBROUTINE W1_COPY_REVERSE_ARG( W2, W1 )
    IMPLICIT NONE
    TYPE (wavefun1), INTENT(IN) :: W1
    TYPE (wavefun1) :: W2

# 403

    CALL ZCOPY( W1%WDES1%NRPLWV, W1%CPTWFP(1), 1, W2%CPTWFP(1), 1)

    IF (W1%WDES1%LGAMMA) THEN
       CALL DCOPY( W1%WDES1%NPROD, W1%CPROJ(1), 1,  W2%CPROJ(1), 1)
    ELSE
       CALL ZCOPY( W1%WDES1%NPROD, W1%CPROJ(1),  1, W2%CPROJ(1), 1)
    ENDIF

    IF (ASSOCIATED(W1%CR) .AND. ASSOCIATED(W2%CR)) THEN
       IF (SIZE(W1%CR) /=W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS) THEN
          CALL vtutor%bug("internal error in W1_COPY: real space allocation is not correct", "wave_high.F", 414)
       ENDIF
       IF (SIZE(W1%CR) /=SIZE(W2%CR)) THEN
          CALL vtutor%bug("internal error in W1_COPY: real space allocation is different " // &
             str(SIZE(W1%CR)) // " " // str(SIZE(W2%CR)), "wave_high.F", 418)
       ENDIF

       CALL ZCOPY( W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS, W1%CR(1), 1, W2%CR(1), 1)
    ENDIF

  END SUBROUTINE W1_COPY_REVERSE_ARG


!***********************************************************************
!
!> copy a WA structure
!> ~~~
!>  W2 = W1  (W1 -> W2)
!> ~~~
!> the argument arrangement is similar to  DCOPY, ZCOPY in BLAS level 1
!> the redistributed wavefunctions pointers in the destination
!> are also properly set
!
!***********************************************************************

  SUBROUTINE WA_COPY( W1, W2 )
    IMPLICIT NONE
    TYPE (wavefuna), INTENT(IN) :: W1
    TYPE (wavefuna) :: W2

    IF (SIZE(W1%CPTWFP) /= SIZE(W2%CPTWFP)) THEN
       CALL vtutor%bug("WA_COPY: size mismatch in CW " // str(SIZE(W1%CPTWFP)) // &
          " " // str(SIZE(W2%CPTWFP)), "wave_high.F", 446)
    ENDIF

    IF (SIZE(W1%CPROJ) /= SIZE(W2%CPROJ)) THEN
       CALL vtutor%bug("WA_COPY: size mismatch in CPROJ " // str(SIZE(W1%CPROJ)) &
          // " " // str(SIZE(W2%CPROJ)), "wave_high.F", 451)
    ENDIF

    CALL ZCOPY( SIZE(W1%CPTWFP), W1%CPTWFP(1,1), 1, W2%CPTWFP(1,1), 1)

    IF (W1%WDES1%LGAMMA) THEN
       CALL DCOPY( SIZE(W1%CPROJ), W1%CPROJ(1,1), 1,  W2%CPROJ(1,1), 1)
    ELSE
       CALL ZCOPY( SIZE(W1%CPROJ), W1%CPROJ(1,1), 1,  W2%CPROJ(1,1), 1)
    ENDIF
! remember WDES1
    W2%WDES1 =>W1%WDES1
! (1._q,0._q) dimensional indexing assumed
    W2%FIRST_DIM=0
! remember spin index
    W2%ISP =W1%ISP
! set redistributed wavefunction indices
    IF (W2%WDES1%DO_REDIS) THEN
       CALL SET_WPOINTER(W2%CW_RED,    W2%WDES1%NRPLWV_RED, W2%WDES1%NB_TOT, W2%CPTWFP(1,1))
       CALL SET_GPOINTER(W2%CPROJ_RED, W2%WDES1%NPROD_RED,  W2%WDES1%NB_TOT, W2%CPROJ(1,1))
    ELSE
       W2%CW_RED=>W2%CPTWFP
       W2%CPROJ_RED=>W2%CPROJ
    ENDIF
  END SUBROUTINE WA_COPY

  SUBROUTINE WA_COPY_CW( W1, W2 )
    IMPLICIT NONE
    TYPE (wavefuna), INTENT(IN) :: W1
    TYPE (wavefuna) :: W2

    IF (SIZE(W1%CPTWFP) /= SIZE(W2%CPTWFP)) THEN
       CALL vtutor%bug("internal error in WA_COPY: size mismatch in CW " // str(SIZE(W1%CPTWFP)) // &
          " " // str(SIZE(W2%CPTWFP)), "wave_high.F", 484)
    ENDIF

    CALL ZCOPY( SIZE(W1%CPTWFP), W1%CPTWFP(1,1), 1, W2%CPTWFP(1,1), 1)

! remember WDES1
    W2%WDES1 =>W1%WDES1
! (1._q,0._q) dimensional indexing assumed
    W2%FIRST_DIM=0
! remember spin index
    W2%ISP =W1%ISP
! set redistributed wavefunction indices
    IF (W2%WDES1%DO_REDIS) THEN
       CALL SET_WPOINTER(W2%CW_RED,    W2%WDES1%NRPLWV_RED, W2%WDES1%NB_TOT, W2%CPTWFP(1,1))
    ELSE
       W2%CW_RED=>W2%CPTWFP
    ENDIF
  END SUBROUTINE WA_COPY_CW

  SUBROUTINE WA_COPY_CPROJ( W1, W2 )
    IMPLICIT NONE
    TYPE (wavefuna), INTENT(IN) :: W1
    TYPE (wavefuna) :: W2

    IF (SIZE(W1%CPROJ) /= SIZE(W2%CPROJ)) THEN
       CALL vtutor%bug("internal error in WA_COPY: size mismatch in CPROJ " // str(SIZE(W1%CPROJ)) &
          // " " // str(SIZE(W2%CPROJ)), "wave_high.F", 510)
    ENDIF

    IF (W1%WDES1%LGAMMA) THEN
       CALL DCOPY( SIZE(W1%CPROJ), W1%CPROJ(1,1), 1,  W2%CPROJ(1,1), 1)
    ELSE
       CALL ZCOPY( SIZE(W1%CPROJ), W1%CPROJ(1,1), 1,  W2%CPROJ(1,1), 1)
      ENDIF
! remember WDES1
    W2%WDES1 =>W1%WDES1
! (1._q,0._q) dimensional indexing assumed
    W2%FIRST_DIM=0
! remember spin index
    W2%ISP =W1%ISP
! set redistributed wavefunction indices
    IF (W2%WDES1%DO_REDIS) THEN
       CALL SET_GPOINTER(W2%CPROJ_RED, W2%WDES1%NPROD_RED,  W2%WDES1%NB_TOT, W2%CPROJ(1,1))
    ELSE
       W2%CPROJ_RED=>W2%CPROJ
    ENDIF
  END SUBROUTINE WA_COPY_CPROJ


!***********************************************************************
!
!> redistribute the wavefunction character for a W1 or WA array
!
!***********************************************************************

  SUBROUTINE W1_REDIS_PROJ( W1)
    IMPLICIT NONE
    TYPE (wavefun1), INTENT(IN) :: W1

    IF (W1%WDES1%DO_REDIS) CALL REDIS_PROJ(W1%WDES1, 1, W1%CPROJ(1))
  END SUBROUTINE W1_REDIS_PROJ

  SUBROUTINE WA_REDIS_PROJ( WA)
    IMPLICIT NONE
    TYPE (wavefuna), INTENT(IN) :: WA

    IF (WA%WDES1%DO_REDIS) CALL REDIS_PROJ(WA%WDES1, SIZE(WA%CPROJ,2), WA%CPROJ(1,1))
  END SUBROUTINE WA_REDIS_PROJ

  SUBROUTINE W_REDIS_PROJ( W)
    IMPLICIT NONE
    TYPE (wavespin), INTENT(IN) :: W
    TYPE (wavedes1)    WDES1          ! descriptor for (1._q,0._q) k-point
    INTEGER :: K1, ISP

    IF (W%WDES%DO_REDIS) THEN
       DO K1=1,W%WDES%NKPTS

          IF (MOD(K1-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1) CYCLE

          CALL SETWDES(W%WDES,WDES1,K1)
          DO ISP=1,W%WDES%ISPIN
             CALL WA_REDIS_PROJ( ELEMENTS( W, WDES1, ISP))
          ENDDO
       ENDDO
    ENDIF
  END SUBROUTINE W_REDIS_PROJ


!***********************************************************************
!
!> redistribute the plane wave coefficients for a W1 or WA array
!
!***********************************************************************

  SUBROUTINE W1_REDIS_PW( W1)
    IMPLICIT NONE
    TYPE (wavefun1), INTENT(IN) :: W1

    IF (W1%WDES1%DO_REDIS) CALL REDIS_PW(W1%WDES1, 1, W1%CPTWFP(1))
  END SUBROUTINE W1_REDIS_PW

  SUBROUTINE WA_REDIS_PW( WA)
    IMPLICIT NONE
    TYPE (wavefuna), INTENT(IN) :: WA

    IF (WA%WDES1%DO_REDIS) CALL REDIS_PW(WA%WDES1, SIZE(WA%CPTWFP,2), WA%CPTWFP(1,1))
  END SUBROUTINE WA_REDIS_PW


  SUBROUTINE W_REDIS_PW( W)
    IMPLICIT NONE
    TYPE (wavespin), INTENT(IN) :: W
    TYPE (wavedes1)    WDES1          ! descriptor for (1._q,0._q) k-point
    INTEGER :: K1, ISP

    IF (W%WDES%DO_REDIS) THEN
       DO K1=1,W%WDES%NKPTS

          IF (MOD(K1-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1) CYCLE

          CALL SETWDES(W%WDES,WDES1,K1)
          DO ISP=1,W%WDES%ISPIN
             CALL WA_REDIS_PW( ELEMENTS( W, WDES1, ISP))
          ENDDO
       ENDDO
    ENDIF
  END SUBROUTINE W_REDIS_PW


!***********************************************************************
!
!> update of vector
!>  ~~~
!>  W2 = W1*a + W2
!>  ~~~
!>
!> @details @ref openmp :
!> BLAS1 calls were replaced by explicit OMP PARALLEL loops
!> (for some reason BLAS1 call thread really badly).
!
!***********************************************************************

  SUBROUTINE W1_DAXPY( W1, SCALE, W2)
    IMPLICIT NONE
    REAL(q) SCALE
    TYPE (wavefun1)    W1, W2
!$  INTEGER I

    


    CALL DAXPY( W1%WDES1%NPL*2, SCALE, W1%CPTWFP(1), 1, W2%CPTWFP(1), 1)

    IF (W1%WDES1%LGAMMA) THEN
       CALL DAXPY( W1%WDES1%NPRO, SCALE, W1%CPROJ(1), 1,  W2%CPROJ(1), 1)
    ELSE
       CALL DAXPY( W1%WDES1%NPRO*2, SCALE, W1%CPROJ(1),  1, W2%CPROJ(1), 1)
    ENDIF

    IF (ASSOCIATED(W1%CR) .AND. ASSOCIATED(W2%CR)) THEN
       IF (SIZE(W1%CR) /=W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS) THEN
          CALL vtutor%bug("internal error in W1_COPY: real space allocation is not correct", "wave_high.F", 646)
       ENDIF

! real space wavefunction complex, SCALE real
       CALL DAXPY( W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS*2, SCALE, W1%CR(1), 1, W2%CR(1), 1)
    ENDIF
# 678


    

  END SUBROUTINE W1_DAXPY

!
!> @details @ref openmp :
!> BLAS1 calls were replaced by explicit OMP PARALLEL loops
!> (for some reason BLAS1 call thread really badly).
!
  SUBROUTINE W1_GAXPY( W1, SCALE, W2)
    IMPLICIT NONE
    COMPLEX(q) SCALE
    TYPE (wavefun1)    W1, W2
!$  INTEGER I

    


    IF (W1%WDES1%LGAMMA) THEN
! wavefunction complex, SCALE real
       CALL DAXPY( W1%WDES1%NPL*2, SCALE, W1%CPTWFP(1), 1, W2%CPTWFP(1), 1)
       CALL DAXPY( W1%WDES1%NPRO, SCALE, W1%CPROJ(1), 1,  W2%CPROJ(1), 1)
    ELSE
! wavefunction complex, SCALE complex
       CALL ZAXPY( W1%WDES1%NPL, SCALE, W1%CPTWFP(1), 1, W2%CPTWFP(1), 1)
       CALL ZAXPY( W1%WDES1%NPRO, SCALE, W1%CPROJ(1),  1, W2%CPROJ(1), 1)
    ENDIF

    IF (ASSOCIATED(W1%CR) .AND. ASSOCIATED(W2%CR)) THEN
       IF (SIZE(W1%CR) /=W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS) THEN
          CALL vtutor%bug("internal error in W1_COPY: real space allocation is not correct", "wave_high.F", 710)
       ENDIF

       IF (W1%WDES1%LGAMMA) THEN
! real space wavefunction complex, SCALE real
          CALL DAXPY( W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS*2, SCALE, W1%CR(1), 1, W2%CR(1), 1)
       ELSE
! real space wavefunction complex, SCALE complex
          CALL ZAXPY( W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS, SCALE, W1%CR(1), 1, W2%CR(1), 1)
       ENDIF
    ENDIF
# 747


    

  END SUBROUTINE W1_GAXPY


!***********************************************************************
!
!> update of vector
!>  ~~~
!>  W1 = W1*a
!>  ~~~
!>
!> @details @ref openmp :
!> BLAS1 calls were replaced by explicit OMP PARALLEL loops
!> (for some reason BLAS1 call thread really badly).
!
!***********************************************************************

  SUBROUTINE W1_DSCAL( W1, SCALE)
    IMPLICIT NONE
    REAL(q) SCALE
    TYPE (wavefun1)    W1
!$  INTEGER I

    

! since this function can be used to (0._q,0._q) out an array we operate
! on all elements (dimension) and not only on those that are
! actually used

    CALL DSCAL( W1%WDES1%NRPLWV*2, SCALE, W1%CPTWFP(1), 1)

    IF (W1%WDES1%LGAMMA) THEN
       CALL DSCAL( W1%WDES1%NPROD, SCALE, W1%CPROJ(1), 1)
    ELSE
       CALL DSCAL( W1%WDES1%NPROD*2, SCALE, W1%CPROJ(1),  1)
    ENDIF

    IF (ASSOCIATED(W1%CR)) THEN
       IF (SIZE(W1%CR) /=W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS) THEN
          CALL vtutor%bug("internal error in W1_COPY: real space allocation is not correct", "wave_high.F", 789)
       ENDIF
! real space wavefunction complex, SCALE real
       CALL DSCAL( W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS*2, SCALE, W1%CR(1), 1)
    ENDIF
# 820


    

  END SUBROUTINE W1_DSCAL


!***********************************************************************
!
!> calculate the dot product between two wavefunctions
!> ~~~
!>  C=   W1^* x W2
!> ~~~
!> this is a substitue for the routine CINDPROD but mind
!> the W1 and W2 are interchanged
!>
!> @details @ref openmp :
!> a ZDOTC call is replaced by an explicit OMP PARALLEL loop
!> (for some reason BLAS1 call thread really badly).
!> The nested loop over \"types\" + \"ions-of-type\" is replaced by
!> a loop over \"all ions\", that is distributed over all available
!> threads.
!
!***********************************************************************

  FUNCTION W1_DOT( W1, W2, CQIJ) RESULT (C)

!! USE moffload

    IMPLICIT NONE
    TYPE (wavefun1)   :: W1, W2
    COMPLEX(q)              :: C
    COMPLEX(q), OPTIONAL :: CQIJ(:,:,:,:) ! optional overlap operator
! local
    COMPLEX(q)    :: CNL, CTMP
    INTEGER :: LMDIM, NPRO, NPRO_, ISPINOR, ISPINOR_, LMMAXC, NT, NI

    REAL(q),    EXTERNAL :: DDOT
    COMPLEX(q), EXTERNAL :: ZDOTC

!$  INTEGER :: I
!$ACC ROUTINE(ECCP_NL) VECTOR

    


    IF (W1%WDES1%LGAMMA) THEN
       C=DDOT(2*W1%WDES1%NPL,W1%CPTWFP(1),1,W2%CPTWFP(1),1)
    ELSE
       C=ZDOTC( W1%WDES1%NPL,W1%CPTWFP(1),1,W2%CPTWFP(1),1)
    ENDIF
# 878

    IF (PRESENT(CQIJ) .AND. W1%WDES1%LOVERL .AND.  W1%WDES1%NPROD>0 ) THEN
       CNL=0
       LMDIM=SIZE(CQIJ,1)
# 886

!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) DEFAULT(NONE) &
!$OMP SHARED(LMDIM,W1,W2,CQIJ) &
!$OMP PRIVATE(ISPINOR,ISPINOR_,NI,NT,LMMAXC,NPRO,NPRO_,CTMP) &
!$OMP REDUCTION(+:CNL)

       spinor: DO ISPINOR=0,W1%WDES1%NRSPINORS-1
       DO ISPINOR_=0,W1%WDES1%NRSPINORS-1
          DO NI=1,W1%WDES1%NIONS
             NT=W1%WDES1%ITYP(NI)
             LMMAXC=W1%WDES1%LMMAX(NT)
             IF (LMMAXC==0) CYCLE
             NPRO =ISPINOR *(W1%WDES1%NPRO/2)+W1%WDES1%LMBASE(NI)
             NPRO_=ISPINOR_*(W1%WDES1%NPRO/2)+W1%WDES1%LMBASE(NI)

             CTMP=0; CALL ECCP_NL(LMDIM,LMMAXC,CQIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),W2%CPROJ(NPRO_+1),W1%CPROJ(NPRO+1),CTMP)
             CNL=CNL+CTMP
          ENDDO
       ENDDO
       ENDDO spinor
# 909

!$OMP END PARALLEL DO

       C=C+CNL
    ENDIF
!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
    CALL M_sum_z(W1%WDES1%COMM_INB, C, 1)

    

  END FUNCTION W1_DOT


!***********************************************************************
!
!> calculate the inproduct between (1._q,0._q) wavefunction and
!> a set of wavefunctions with possible scales and add the result
!> to a third vector
!>  ~~~
!>  C=   WA^* x W1*SCALEA+C * SCALEC
!>  ~~~
!
!***********************************************************************

  SUBROUTINE W1_GEMV( SCALEA, WA, W1, SCALEC, C, NC , CQIJ)

!! USE moffload

    IMPLICIT NONE
!> scaleing constant for WA x W1
    COMPLEX(q):: SCALEA
    TYPE (wavefuna)    WA
    TYPE (wavefun1)    W1
!> scaling constant for C
    COMPLEX(q):: SCALEC
!> result
    COMPLEX(q):: C(*)
!> stride of C
    INTEGER NC
!> optional overlap operator
    COMPLEX(q), OPTIONAL:: CQIJ(:,:,:,:)
! local
    INTEGER :: LMDIM, I, N
    COMPLEX(q) CRESUL(WA%WDES1%NPRO)
    COMPLEX(q) CTMP(SIZE(WA%CPTWFP,2))

    

    N=SIZE(WA%CPTWFP,2)

!$ACC ENTER DATA CREATE(CTMP) 
    IF (W1%WDES1%LGAMMA) THEN
       CALL DGEMV( 'T', 2* WA%WDES1%NPL, N, 1._q , WA%CPTWFP(1,1) , &
                        2* WA%WDES1%NRPLWV, W1%CPTWFP(1) , 1 , 0._q,  CTMP(1), 1)
    ELSE
       CALL ZGEMV( 'C',  WA%WDES1%NPL, N, (1._q,0._q) , WA%CPTWFP(1,1) , &
                         WA%WDES1%NRPLWV, W1%CPTWFP(1) , 1 , (0._q,0._q),  CTMP(1), 1)
    ENDIF

    IF (PRESENT(CQIJ) .AND. WA%WDES1%LOVERL .AND.  WA%WDES1%NPROD>0 ) THEN
       LMDIM = SIZE(CQIJ,1)
!$ACC ENTER DATA CREATE(CRESUL) 
!       CALL OVERL1(WA%WDES1, LMDIM, CQIJ(1,1,1,1), CQIJ(1,1,1,1), 0.0_q, WA%CPROJ(1,1), CRESUL(1))
       CALL OVERL1(WA%WDES1, LMDIM, CQIJ(1,1,1,1), CQIJ(1,1,1,1), 0.0_q, W1%CPROJ(1), CRESUL(1))
       IF (W1%WDES1%LGAMMA) THEN
          CALL DGEMV( 'T', WA%WDES1%NPRO,  N, 1._q ,WA%CPROJ(1,1), &
                           WA%WDES1%NPROD, CRESUL , 1 , 1._q,  CTMP(1), 1)
       ELSE
          CALL ZGEMV( 'C', WA%WDES1%NPRO,  N, (1._q,0._q) ,WA%CPROJ(1,1), &
                           WA%WDES1%NPROD, CRESUL , 1 , (1._q,0._q),  CTMP(1), 1)
       ENDIF
!$ACC EXIT DATA DELETE(CRESUL) 
    ENDIF
    CALL M_sum_z(WA%WDES1%COMM_INB, CTMP, N)

    IF (SCALEC==0) THEN
!$ACC ENTER DATA CREATE(C(1:(N-1)*NC+1)) 
!$ACC PARALLEL LOOP PRESENT(C,CTMP) 
       DO I=0,N-1
          C(I*NC+1)=CTMP(I+1)*SCALEA
       ENDDO
    ELSE
!$ACC ENTER DATA COPYIN(C(1:(N-1)*NC+1)) 
!$ACC PARALLEL LOOP PRESENT(C,CTMP) 
       DO I=0,N-1
          C(I*NC+1)=C(I*NC+1)*SCALEC+CTMP(I+1)*SCALEA
       ENDDO
    ENDIF
!$ACC EXIT DATA COPYOUT(C(1:(N-1)*NC+1)) DELETE(CTMP) 

    

  END SUBROUTINE W1_GEMV


!***********************************************************************
!
!> W1 descriptor from W array
!>
!> alternative to SETWAV (returns a W1 descriptor)
!> it is somewhat slimmed down to optimize performance
!> and returns only the wavefunction and wavefunction character pointers
!> the full version remains SETWAV
!
!***********************************************************************

  FUNCTION W1_FROM_W( W, WDES1, NB, ISP) RESULT (W1)
    IMPLICIT NONE
    INTEGER NB, ISP
    TYPE (wavespin) W
    TYPE (wavefun1) W1
    TYPE (wavedes1), TARGET :: WDES1
    INTEGER NK

    NK=WDES1%NK

    IF (NB<=0 .OR. NB> SIZE(W%CPTWFP,2)) THEN
       CALL vtutor%bug("internal error in W1_FROM_W: bounds exceed " // str(NB) // " " // &
          str(SIZE(W%CPTWFP,2)), "wave_high.F", 1027)
    ENDIF

    W1%CPTWFP=>W%CPTWFP(:,NB,NK,ISP)
    W1%CPROJ =>W%CPROJ(:,NB,NK,ISP)
    W1%WDES1 => WDES1
    NULLIFY(W1%CR)
    W1%LDO=.TRUE.
  END FUNCTION W1_FROM_W


!***********************************************************************
!
!> W1 descriptor from WA descriptor
!
!***********************************************************************

  FUNCTION W1_FROM_WA( WA, N1) RESULT (W1)
    IMPLICIT NONE
    INTEGER N1
    TYPE (wavefun1) W1
    TYPE (wavefuna) WA

    IF (N1<=0 .OR. N1> SIZE(WA%CPTWFP,2)) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1) // " " // &
          str(SIZE(WA%CPTWFP,2)), "wave_high.F", 1052)
    ENDIF

    W1%CPTWFP=>WA%CPTWFP(:,N1)
    W1%CPROJ =>WA%CPROJ(:,N1)
    W1%WDES1 =>WA%WDES1

    IF (ASSOCIATED(WA%CR)) THEN
       W1%CR =>WA%CR(:,N1)
    ELSE
       NULLIFY(W1%CR)
    ENDIF

    W1%LDO=.TRUE.
  END FUNCTION W1_FROM_WA


!***********************************************************************
!
!> subindex a WA array
!> return WA(N1:N2)
!
!***********************************************************************

  FUNCTION WA_FROM_WA( WA, N1, N2) RESULT (W1)
    IMPLICIT NONE
    INTEGER N1, N2
    TYPE (wavefuna) W1
    TYPE (wavefuna) WA

    IF (N1<=0 .OR. N1> SIZE(WA%CPTWFP,2)) THEN
       CALL vtutor%bug("internal error in WA_FROM_WA: bounds exceed " // str(N1) // " " // &
          str(SIZE(WA%CPTWFP,2)), "wave_high.F", 1084)
    ENDIF
    IF (N2<N1 .OR. N2> SIZE(WA%CPTWFP,2)) THEN
       CALL vtutor%bug("internal error in WA_FROM_WA: bounds exceed " // str(N1) // " " // str(N2) &
          // " " // str(SIZE(WA%CPTWFP,2)), "wave_high.F", 1088)
    ENDIF

    W1%CPTWFP=>WA%CPTWFP(:,N1:N2)
    W1%CPROJ =>WA%CPROJ(:,N1:N2)
    W1%WDES1 =>WA%WDES1
  END FUNCTION WA_FROM_WA


!***********************************************************************
!
!> (1._q,0._q) element from a WA array using two indices
!> return WA(N1, N2)
!
!***********************************************************************

  FUNCTION W1_FROM_WA2( WA, N1, N2) RESULT (W1)
    IMPLICIT NONE
    INTEGER N1, N2
    TYPE (wavefun1) W1
    TYPE (wavefuna) WA

    IF (N1<=0 .OR. N1> WA%FIRST_DIM) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1) // " " // &
          str(WA%FIRST_DIM), "wave_high.F", 1112)
    ENDIF
    IF (N1+ (N2-1)*WA%FIRST_DIM<=0 .OR. N1+ (N2-1)*WA%FIRST_DIM> SIZE(WA%CPTWFP,2)) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1+(N2-1)) // " " // &
          str(SIZE(WA%CPTWFP,2)), "wave_high.F", 1116)
    ENDIF

    W1%CPTWFP=>WA%CPTWFP(:,N1+ (N2-1)*WA%FIRST_DIM)
    W1%CPROJ =>WA%CPROJ (:,N1+ (N2-1)*WA%FIRST_DIM)
    W1%WDES1 =>WA%WDES1
    NULLIFY(W1%CR)
    W1%LDO=.TRUE.
  END FUNCTION W1_FROM_WA2


!***********************************************************************
!
!> subindex a WA array using two indices
!> return WA(N1:N12, N2)
!
!***********************************************************************

  FUNCTION WA_FROM_WA2( WA, N1, N12, N2) RESULT (W1)
    IMPLICIT NONE
    INTEGER N1, N12, N2
    TYPE (wavefuna) W1
    TYPE (wavefuna) WA

    IF (N1<=0 .OR. N1> WA%FIRST_DIM) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1) // " " // &
          str(WA%FIRST_DIM), "wave_high.F", 1142)
    ENDIF
    IF (N1+ (N2-1)*WA%FIRST_DIM<=0 .OR. N1+ (N2-1)*WA%FIRST_DIM> SIZE(WA%CPTWFP,2)) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1+(N2-1)) // " " // &
          str(SIZE(WA%CPTWFP,2)), "wave_high.F", 1146)
    ENDIF
    IF (N12<N1 .OR. N12> WA%FIRST_DIM) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1) // " " // str(N12) &
          // " " // str(WA%FIRST_DIM), "wave_high.F", 1150)
    ENDIF
    IF (N12+ (N2-1)*WA%FIRST_DIM<N1+ (N2-1)*WA%FIRST_DIM .OR. N12+ (N2-1)*WA%FIRST_DIM> SIZE(WA%CPTWFP,2)) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N12+(N2-1)) // " " // &
          str(SIZE(WA%CPTWFP,2)), "wave_high.F", 1154)
    ENDIF


    W1%CPTWFP=>WA%CPTWFP(:,N1+ (N2-1)*WA%FIRST_DIM:N12+ (N2-1)*WA%FIRST_DIM)
    W1%CPROJ =>WA%CPROJ (:,N1+ (N2-1)*WA%FIRST_DIM:N12+ (N2-1)*WA%FIRST_DIM)
    W1%WDES1 =>WA%WDES1
  END FUNCTION WA_FROM_WA2


!***********************************************************************
!
!> get a WA structure from a wavefunction array
!> this is equivalent to the low level routine SETWAVA
!
!***********************************************************************

  FUNCTION WA_FROM_W( W, WDES1, ISP ) RESULT (WA)
    IMPLICIT NONE
    TYPE (wavespin), INTENT(IN) :: W
    TYPE (wavedes1), TARGET ::  WDES1
    INTEGER ISP
    TYPE (wavefuna) ::  WA

    WA%CPTWFP=>W%CPTWFP(:,:,WDES1%NK,ISP)
    WA%CPROJ =>W%CPROJ(:,:,WDES1%NK,ISP)
    WA%FERWE =>W%FERWE(:,WDES1%NK,ISP)
    WA%AUX   =>W%AUX  (:,WDES1%NK,ISP)
    WA%CELEN =>W%CELEN(:,WDES1%NK,ISP)
    WA%WDES1 =>WDES1
! remember spin index
    WA%ISP =ISP
! set redistributed wavefunction indices
    IF (WDES1%DO_REDIS) THEN
       CALL SET_WPOINTER(WA%CW_RED,    WDES1%NRPLWV_RED, W%WDES%NB_TOT, WA%CPTWFP(1,1))
       CALL SET_GPOINTER(WA%CPROJ_RED, WDES1%NPROD_RED,  W%WDES%NB_TOT, WA%CPROJ(1,1))
    ELSE
       WA%CW_RED=>WA%CPTWFP
       WA%CPROJ_RED=>WA%CPROJ
    ENDIF
  END FUNCTION WA_FROM_W


!************************* SUBROUTINE NEWWAVA **************************
!
!>  create storage for a wavefunction array WA
!
!***********************************************************************

  SUBROUTINE NEWWAVA(WA, WDES1, NDIM, NDIM2, ALLOC_REAL)

!! USE moffload

    IMPLICIT NONE
    TYPE (wavefuna), TARGET :: WA
    TYPE (wavedes1), TARGET :: WDES1
    INTEGER :: NDIM
    INTEGER, OPTIONAL :: NDIM2
    LOGICAL, OPTIONAL :: ALLOC_REAL
! local variables
    INTEGER :: NDIMTOT

    NDIMTOT=NDIM
    WA%FIRST_DIM=0

    IF (PRESENT(NDIM2)) THEN
       NDIMTOT=NDIM*NDIM2
       WA%FIRST_DIM=NDIM
    ENDIF

    ALLOCATE(WA%CPTWFP(WDES1%NRPLWV,NDIMTOT), WA%CPROJ(WDES1%NPROD,NDIMTOT), &
             WA%CELEN(NDIMTOT), WA%FERWE(NDIMTOT), WA%AUX(NDIMTOT))
!$ACC ENTER DATA CREATE(WA%CPTWFP,WA%CPROJ) IF(OFFLOAD_ON.AND.ACC_IS_PRESENT(WA,1)) ASYNC(ACC_ASYNC_Q)

! WA%CPTWFP=0; WA%CPROJ=0 ! Breaks ???.

! set redistributed wavefunction indices
    IF (WDES1%DO_REDIS) THEN
       CALL SET_WPOINTER(WA%CW_RED,    WDES1%NRPLWV_RED, NDIM*WDES1%NB_PAR, WA%CPTWFP(1,1))
       CALL SET_GPOINTER(WA%CPROJ_RED, WDES1%NPROD_RED,  NDIM*WDES1%NB_PAR, WA%CPROJ(1,1))
    ELSE
       WA%CW_RED   =>WA%CPTWFP
       WA%CPROJ_RED=>WA%CPROJ
    ENDIF
!$ACC ENTER DATA CREATE(WA%CW_RED,WA%CPROJ_RED) IF(OFFLOAD_ON.AND.ACC_IS_PRESENT(WA,1)) ASYNC(ACC_ASYNC_Q)

! optionally allocate space for the orbitals in real space
    IF (PRESENT(ALLOC_REAL)) THEN
       IF (ALLOC_REAL) THEN
          ALLOCATE(WA%CR(WDES1%GRID%MPLWV*WDES1%NRSPINORS,NDIMTOT))
!$ACC ENTER DATA CREATE(WA%CR) IF(OFFLOAD_ON.AND.ACC_IS_PRESENT(WA,1)) ASYNC(ACC_ASYNC_Q)
       ENDIF
    ENDIF

    WA%ISP=-1

    WA%WDES1=>WDES1
# 1258

  END SUBROUTINE NEWWAVA


!************************* SUBROUTINE NEWWAVA_PROJ *********************
!
!>  create storage for (1._q,0._q) wavefunction W array
!>  to store non local part only
!
!***********************************************************************

  SUBROUTINE NEWWAVA_PROJ(WA, WDES1, NDIM)

!! USE moffload

    USE prec
    USE wave_mpi
    IMPLICIT NONE
    TYPE (wavefuna), TARGET :: WA
    TYPE (wavedes1), TARGET :: WDES1
    INTEGER, OPTIONAL :: NDIM
    IF (PRESENT(NDIM)) THEN
       ALLOCATE(WA%CPROJ(WDES1%NPROD,NDIM))
    ELSE
       ALLOCATE(WA%CPROJ(WDES1%NPROD,WDES1%NBANDS))
    ENDIF
!$ACC ENTER DATA CREATE(WA%CPROJ) IF(OFFLOAD_ON.AND.ACC_IS_PRESENT(WA,1)) ASYNC(ACC_ASYNC_Q)
    WA%FIRST_DIM=0
    WA%ISP=-1
    NULLIFY(WA%CPTWFP)
    NULLIFY(WA%CW_RED)
    NULLIFY(WA%FERWE)
    NULLIFY(WA%AUX  )
    NULLIFY(WA%CELEN)
! set redistributed projection indices
    IF (WDES1%DO_REDIS) THEN
       IF (PRESENT(NDIM)) THEN
          CALL SET_GPOINTER(WA%CPROJ_RED, WDES1%NPROD_RED,  NDIM*WDES1%NB_PAR, WA%CPROJ(1,1))
       ELSE
          CALL SET_GPOINTER(WA%CPROJ_RED, WDES1%NPROD_RED,  WDES1%NB_TOT, WA%CPROJ(1,1))
       ENDIF
    ELSE
       WA%CPROJ_RED=>WA%CPROJ
    ENDIF
!$ACC ENTER DATA CREATE(WA%CPROJ_RED) IF(OFFLOAD_ON.AND.ACC_IS_PRESENT(WA,1)) ASYNC(ACC_ASYNC_Q)
    WA%WDES1=>WDES1
# 1311

  END SUBROUTINE NEWWAVA_PROJ


!************************* SUBROUTINE NEWWAV_ARRAY *********************
!
!>  create storage for an array of single orbital structures W1(:)
!
!***********************************************************************

  SUBROUTINE NEWWAV_ARRAY(W1,WDES1,ALLOC_REAL,WA)

!! USE moffload

    IMPLICIT NONE

    TYPE (wavefun1) :: W1(:)
    TYPE (wavedes1) :: WDES1
    TYPE (wavefuna), OPTIONAL :: WA

    LOGICAL :: ALLOC_REAL

! local variables
    INTEGER :: NDIM,I

    NDIM=SIZE(W1)

    IF (PRESENT(WA)) THEN
!
! Allocate a wavefuna structure WA to create storage space for the pointer
! members of the array of single orbital structures W1(:).
! This guarantees that the respective orbital data arrays W1(:)%CPTWFP and
! W1(:)%CR will be laid out contiguously in memory.
!
!$ACC ENTER DATA CREATE(WA) IF(OFFLOAD_ON.AND.ACC_IS_PRESENT(W1(:),1)) ASYNC(ACC_ASYNC_Q)
! allocate a wavefuna structure WA of size NDIM
       CALL NEWWAVA(WA, WDES1, NDIM, ALLOC_REAL=ALLOC_REAL)
! connect the W1 members to WA
       DO I=1,NDIM
! connect pointers
          W1(I) = ELEMENT(WA, I)
!$ACC ENTER DATA COPYIN(W1(I)%CPTWFP,W1(I)%CPROJ) IF(OFFLOAD_ON.AND.ACC_IS_PRESENT(W1(I),1)) ASYNC(ACC_ASYNC_Q)
!$ACC ENTER DATA COPYIN(W1(I)%CR) IF(OFFLOAD_ON.AND.ACC_IS_PRESENT(W1(I),1).AND.ALLOC_REAL) ASYNC(ACC_ASYNC_Q)
!!!$ACC ENTER DATA CREATE(W1(I)%CPTWFP,W1(I)%CPROJ) IF(OFFLOAD_ON.AND.ACC_IS_PRESENT(W1(I),1)) ASYNC(ACC_ASYNC_Q)
!!!$ACC ENTER DATA CREATE(W1(I)%CR) IF(OFFLOAD_ON.AND.ACC_IS_PRESENT(W1(I),1).AND.ALLOC_REAL) ASYNC(ACC_ASYNC_Q)

! set other members of W1 (see NEWWAV)
          W1(I)%FERWE= 0
          W1(I)%CELEN= 0
          W1(I)%NB   =-1
          W1(I)%ISP  =-1
          W1(I)%LDO  = .TRUE.

# 1369

       ENDDO
    ELSE
!
! In case WA is not present we fall back to the use of NEWWAV.
!
       DO I=1,NDIM
          CALL NEWWAV(W1(I), WDES1, ALLOC_REAL)
       ENDDO
    ENDIF
  END SUBROUTINE NEWWAV_ARRAY


!************************* SUBROUTINE DELWAV_ARRAY *********************
!
!>  destroy storage for an array of single orbital structures W1(:)
!>  that was possibly created using a wavefuna structure
!
!***********************************************************************

  SUBROUTINE DELWAV_ARRAY(W1,DEALLOC_REAL,WA)

!! USE moffload

    IMPLICIT NONE

    TYPE (wavefun1) :: W1(:)
    TYPE (wavefuna), OPTIONAL :: WA

    LOGICAL :: DEALLOC_REAL

! local variables
    INTEGER :: NDIM,I

    NDIM=SIZE(W1)

    IF (PRESENT(WA)) THEN
       IF (NDIM/=SIZE(WA%CPTWFP,2)) CALL vtutor%bug("DELWAV_ARRAY: partial deallocation is not supported","wave_high.F",1406)
! deallocate WA structure
       CALL DELWAVA(WA)
!$ACC EXIT DATA DELETE(WA) IF(OFFLOAD_ON.AND.ACC_IS_PRESENT(WA,1)) ASYNC(ACC_ASYNC_Q)
! nullify W1 pointer members
       DO I=1,NDIM
!$ACC EXIT DATA DELETE(W1(I)%CPTWFP,W1(I)%CPROJ,W1(I)%WDES1) IF(OFFLOAD_ON.AND.ACC_IS_PRESENT(W1(I),1)) ASYNC(ACC_ASYNC_Q)
          NULLIFY(W1(I)%CPTWFP,W1(I)%CPROJ,W1(I)%WDES1)
          IF (DEALLOC_REAL) THEN
!$ACC EXIT DATA DELETE(W1(I)%CR) IF(OFFLOAD_ON.AND.ACC_IS_PRESENT(W1(I),1)) ASYNC(ACC_ASYNC_Q)
             NULLIFY(W1(I)%CR)
          ENDIF
       ENDDO
    ELSE
       DO I=1,NDIM
          CALL DELWAV(W1(I), DEALLOC_REAL)
       ENDDO
    ENDIF
  END SUBROUTINE DELWAV_ARRAY


!************************* FUNCTION WA_CONNECT2_W1 *********************
!
!>  Connect the pointer members of a wavefuna structure to the members
!>  of an array of single orbital structures
!>
!>  Mind: the latter has to have been created using NEWWAV_ARRAY, to
!>  ensure that things like W1(:)%CPTWFP are contiguous in memory
!>
!>  WA_CONNECT2_W1 will try to make sure that the members of WA and W1(:)
!>  are properly mapped
!
!***********************************************************************

  SUBROUTINE WA_CONNECT2_W1(W1,WA)

!! USE moffload

    USE iso_c_binding
    IMPLICIT NONE
    TYPE (wavefun1) :: W1(:)
    TYPE (wavefuna) :: WA

! local variables
    INTEGER :: I

    CALL c_f_pointer(c_loc(W1(1)%CPTWFP(1)),WA%CPTWFP,[SIZE(W1(1)%CPTWFP),SIZE(W1)])
    CALL c_f_pointer(c_loc(W1(1)%CPROJ (1)),WA%CPROJ ,[SIZE(W1(1)%CPROJ ),SIZE(W1)])

    IF (ASSOCIATED(W1(1)%CR)) CALL c_f_pointer(c_loc(W1(1)%CR(1)),WA%CR,[SIZE(W1(1)%CR),SIZE(W1)])

    WA%WDES1 => W1(1)%WDES1

! test consistency
    DO I=1,SIZE(W1)
       IF (.NOT.c_associated(c_loc(WA%CPTWFP(1,I)),c_loc(W1(I)%CPTWFP(1)))) &
          CALL vtutor%bug("WA_CONNECT2_W1: target member CPTWFP appears to be non-contiguous","wave_high.F",1462)
       IF (.NOT.c_associated(c_loc(WA%CPROJ(1,I)),c_loc(W1(I)%CPROJ(1)))) &
          CALL vtutor%bug("WA_CONNECT2_W1: target member CPROJ appears to be non-contiguous", "wave_high.F",1464)
        IF (ASSOCIATED(WA%CR) .AND. .NOT.c_associated(c_loc(WA%CR(1,I)),c_loc(W1(I)%CR(1)))) &
          CALL vtutor%bug("WA_CONNECT2_W1: target member CR appears to be non-contiguous",    "wave_high.F",1466)
    ENDDO

!$ACC ENTER DATA COPYIN(WA) IF(OFFLOAD_ON.AND.ACC_IS_PRESENT(W1,1)) ASYNC(ACC_ASYNC_Q)
!$ACC ENTER DATA COPYIN(WA%CPTWFP,WA%CPROJ,WA%WDES1) IF(OFFLOAD_ON.AND.ACC_IS_PRESENT(W1,1)) ASYNC(ACC_ASYNC_Q)
!$ACC ENTER DATA COPYIN(WA%CR) IF(OFFLOAD_ON.AND.ACC_IS_PRESENT(W1,1).AND.ASSOCIATED(WA%CR)) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE WA_CONNECT2_W1


!************************* FUNCTION WA_CONNECT2_AR *********************
!
!>  Allocate a wavefun1 array and connect its pointer members to
!>  two 2d-arrays, CPTWFP and CPROJ, and a wavedes1 structure.
!>  Mind the leading dimensions of CPTWFP and CPROJ must be larger than
!>  or equal to WDES1%NRPLWV and WDES1%NPROD, respectively.
!>  The second dimensions of CPTWFP and CPROJ must be equal.
!
!***********************************************************************

  SUBROUTINE WA_CONNECT2_AR(CPTWFP,CPROJ,WDES1,WA)

!! USE moffload

    IMPLICIT NONE
    COMPLEX(q), TARGET :: CPTWFP(:,:)
    COMPLEX(q),       TARGET :: CPROJ(:,:)

    TYPE (wavedes1), TARGET :: WDES1
    TYPE (wavefuna) :: WA

! local
    INTEGER :: NDIM2, I

# 1502


    IF (SIZE(CPTWFP,2)/=SIZE(CPROJ,2)) &
       CALL vtutor%bug("WA_CONNECT2_AR: array sizes inconsistent","wave_high.F",1505)

    IF (SIZE(CPTWFP,1)<WDES1%NRPLWV) &
       CALL vtutor%bug("WA_CONNECT2_AR: first dimension of CW too small "// &
          str(SIZE(CPTWFP,1))//" "//str(WDES1%NRPLWV),"wave_high.F",1509)

    IF (SIZE(CPROJ,1)<WDES1%NPROD) &
       CALL vtutor%bug("WA_CONNECT2_AR: first dimension of CPROJ too small "// &
          str(SIZE(CPROJ,1))//" "//str(WDES1%NPROD),"wave_high.F",1513)

    WA%CPTWFP=>CPTWFP
    WA%CPROJ =>CPROJ
    WA%WDES1 =>WDES1

!$ACC ENTER DATA CREATE(WA) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
!$ACC ENTER DATA COPYIN(WA%CPTWFP,WA%CPROJ,WA%WDES1) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE WA_CONNECT2_AR


!************************* SUBROUTINE WA_NULLIFY ***********************
!
!>  nullify the pointer members of a wavefuna structure
!
!***********************************************************************

  SUBROUTINE WA_NULLIFY(WA)

!! USE moffload

    IMPLICIT NONE
    TYPE (wavefuna) :: WA
!$ACC EXIT DATA DELETE(WA%CR) IF(OFFLOAD_ON.AND.ACC_IS_PRESENT(WA,1).AND.ASSOCIATED(WA%CR)) ASYNC(ACC_ASYNC_Q)
!$ACC EXIT DATA DELETE(WA%CPTWFP,WA%CPROJ,WA%WDES1) IF(OFFLOAD_ON.AND.ACC_IS_PRESENT(WA,1)) ASYNC(ACC_ASYNC_Q)
!$ACC EXIT DATA DELETE(WA) IF(OFFLOAD_ON.AND.ACC_IS_PRESENT(WA,1))  ASYNC(ACC_ASYNC_Q)
    NULLIFY(WA%CPTWFP,WA%CPROJ,WA%WDES1); IF (ASSOCIATED(WA%CR)) NULLIFY(WA%CR)
  END SUBROUTINE WA_NULLIFY


!************************* FUNCTION W1_CONNECT2_AR *********************
!
!>  Allocate a wavefun1 array and connect its pointer members to
!>  two 2d-arrays, CPTWFP and CPROJ, and a wavedes1 structure.
!>
!>  Mind the leading dimensions of CPTWFP and CPROJ must be larger than
!>  or equal to WDES1%NRPLWV and WDES1%NPROD, respectively.
!>  The second dimensions of CPTWFP and CPROJ must be equal.
!
!***********************************************************************

  SUBROUTINE W1_CONNECT2_AR(CPTWFP,CPROJ,WDES1,W1)

!! USE moffload

    IMPLICIT NONE
    COMPLEX(q), TARGET :: CPTWFP(:,:)
    COMPLEX(q),       TARGET :: CPROJ(:,:)

    TYPE (wavedes1), TARGET  :: WDES1
    TYPE (wavefun1), POINTER :: W1(:)

! local
    INTEGER :: NDIM2, I

# 1571


    IF (SIZE(CPTWFP,2)/=SIZE(CPROJ,2)) &
       CALL vtutor%bug("W1_CONNECT2_AR: array sizes inconsistent","wave_high.F",1574)

    IF (SIZE(CPTWFP,1)<WDES1%NRPLWV) &
       CALL vtutor%bug("W1_CONNECT2_AR: first dimension of CW too small "// &
          str(SIZE(CPTWFP,1))//" "//str(WDES1%NRPLWV),"wave_high.F",1578)

    IF (SIZE(CPROJ,1)<WDES1%NPROD) &
       CALL vtutor%bug("W1_CONNECT2_AR: first dimension of CPROJ too small "// &
          str(SIZE(CPROJ,1))//" "//str(WDES1%NPROD),"wave_high.F",1582)

    NDIM2=SIZE(CPTWFP,2); ALLOCATE(W1(NDIM2))
!$ACC ENTER DATA CREATE(W1(:)) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
    DO I=1,NDIM2
       W1(I)%CPTWFP => CPTWFP(:,I)
       W1(I)%CPROJ  => CPROJ(:,I)
       W1(I)%WDES1  => WDES1
       W1(I)%LDO    =  .TRUE.
!$ACC ENTER DATA COPYIN(W1(I)%CPTWFP,W1(I)%CPROJ,W1(I)%WDES1) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
!$ACC UPDATE DEVICE(W1(I)%LDO) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
    ENDDO
  END SUBROUTINE W1_CONNECT2_AR


!************************* SUBROUTINE W1_NULLIFY ***********************
!
!>  nullify the pointer members of a wavefun1 array
!
!***********************************************************************

  SUBROUTINE W1_NULLIFY(W1)

!! USE moffload

    IMPLICIT NONE
    TYPE (wavefun1), POINTER :: W1(:)
! local
    INTEGER :: I
    DO I=1,SIZE(W1)
!$ACC EXIT DATA DELETE(W1(I)%CPTWFP,W1(I)%CPROJ,W1(I)%WDES1) IF(OFFLOAD_ON.AND.ACC_IS_PRESENT(W1(I),1)) ASYNC(ACC_ASYNC_Q)
       NULLIFY(W1(I)%CPTWFP,W1(I)%CPROJ,W1(I)%WDES1)
    ENDDO
!$ACC EXIT DATA DELETE(W1(:)) IF(OFFLOAD_ON.AND.ACC_IS_PRESENT(W1,1)) ASYNC(ACC_ASYNC_Q)
    DEALLOCATE(W1)
  END SUBROUTINE W1_NULLIFY


!************************* SUBROUTINE FFTWAV_WA ************************
!
!>  batched FFT (q -> r) of a group of orbitals
!
!***********************************************************************

  SUBROUTINE FFTWAV_WA(WA)
    IMPLICIT NONE
    TYPE (wavefuna) :: WA
! local variables
    INTEGER :: ISPINOR

    IF (WA%WDES1%NK==0) &
       CALL vtutor%bug("FFTWAV_WA: NK is set to zero: WDES not set up properly", "wave_high.F", 1633)

    IF (.NOT.ASSOCIATED(WA%CR)) &
       CALL vtutor%bug("FFTWAV_WA: can not perform FFT: CR is not allocated", "wave_high.F", 1636)

    IF (WA%WDES1%LUSEINV) THEN
       DO ISPINOR=0,WA%WDES1%NRSPINORS-1
          CALL FFTWAV_USEINV_MU(WA%WDES1%NGVECTOR,SIZE(WA%CPTWFP,2),WA%WDES1%NINDPW(1), &
                                WA%WDES1%NINDPW_INV(1),WA%WDES1%FFTSCA(1,2), &
                                WA%CR(1+ISPINOR*WA%WDES1%GRID%MPLWV,1),SIZE(WA%CR,1), &
                                WA%CPTWFP(1+ISPINOR*WA%WDES1%NGVECTOR,1),SIZE(WA%CPTWFP,1), &
                                WA%WDES1%GRID)
       ENDDO
    ELSE
       DO ISPINOR=0,WA%WDES1%NRSPINORS-1
          CALL FFTWAV_MU(WA%WDES1%NGVECTOR,SIZE(WA%CPTWFP,2),WA%WDES1%NINDPW(1), &
                         WA%CR(1+ISPINOR*WA%WDES1%GRID%MPLWV,1),SIZE(WA%CR,1), &
                         WA%CPTWFP(1+ISPINOR*WA%WDES1%NGVECTOR,1),SIZE(WA%CPTWFP,1), &
                         WA%WDES1%GRID)
       ENDDO
    ENDIF
  END SUBROUTINE FFTWAV_WA


!************************* SUBROUTINE FFTEXT_WA ************************
!
!>  batched FFT (r -> q) of a group of orbitals
!
!***********************************************************************

  SUBROUTINE FFTEXT_WA(WA,LADD)
    IMPLICIT NONE
    TYPE (wavefuna) :: WA
    LOGICAL :: LADD
! local variables
    INTEGER :: ISPINOR

    IF (WA%WDES1%NK==0) &
       CALL vtutor%bug("FFTEXT_WA: NK is set to zero: WDES not set up properly", "wave_high.F", 1671)

    IF (.NOT.ASSOCIATED(WA%CR)) &
       CALL vtutor%bug("FFTEXT_WA: can not perform FFT: CR is not allocated", "wave_high.F", 1674)

    IF (WA%WDES1%LUSEINV) THEN
       DO ISPINOR=0,WA%WDES1%NRSPINORS-1
          CALL FFTEXT_USEINV_MU(WA%WDES1%NGVECTOR,SIZE(WA%CPTWFP,2),WA%WDES1%NINDPW(1),WA%WDES1%FFTSCA(1,1), &
                                WA%CR(1+ISPINOR*WA%WDES1%GRID%MPLWV,1),SIZE(WA%CR,1), &
                                WA%CPTWFP(1+ISPINOR*WA%WDES1%NGVECTOR,1),SIZE(WA%CPTWFP,1), &
                                WA%WDES1%GRID,LADD)
       ENDDO
    ELSE
       DO ISPINOR=0,WA%WDES1%NRSPINORS-1
          CALL FFTEXT_MU(WA%WDES1%NGVECTOR,SIZE(WA%CPTWFP,2),WA%WDES1%NINDPW(1), &
                         WA%CR(1+ISPINOR*WA%WDES1%GRID%MPLWV,1),SIZE(WA%CR,1), &
                         WA%CPTWFP(1+ISPINOR*WA%WDES1%NGVECTOR,1),SIZE(WA%CPTWFP,1), &
                         WA%WDES1%GRID,LADD)
       ENDDO
    ENDIF
  END SUBROUTINE FFTEXT_WA


!************************* SUBROUTINE SETWAVA **************************
!
!>  set (1._q,0._q) single wavefunction array (WA) from an array of wavefunctions
!
!***********************************************************************

  SUBROUTINE SETWAVA(W, WA, WDES1, ISP)
    USE prec
    IMPLICIT NONE
    INTEGER NB,ISP
    TYPE (wavespin) W
    TYPE (wavefuna) WA
    TYPE (wavedes1), TARGET :: WDES1
    INTEGER NK

    NK=WDES1%NK

    WA%CPTWFP=>W%CPTWFP(:,:,NK,ISP)
    WA%CPROJ =>W%CPROJ(:,:,NK,ISP)
    WA%FERWE =>W%FERWE(:,NK,ISP)
    WA%AUX   =>W%AUX  (:,NK,ISP)
    WA%CELEN =>W%CELEN(:,NK,ISP)
    WA%WDES1 =>WDES1
! (1._q,0._q) dimensional indexing assumed
    WA%FIRST_DIM=0
! remember spin index
    WA%ISP =ISP
! set redistributed wavefunction indices
    IF (WDES1%DO_REDIS) THEN
       CALL SET_WPOINTER(WA%CW_RED,    WA%WDES1%NRPLWV_RED, W%WDES%NB_TOT, WA%CPTWFP(1,1))
       CALL SET_GPOINTER(WA%CPROJ_RED, WA%WDES1%NPROD_RED,  W%WDES%NB_TOT, WA%CPROJ(1,1))
    ELSE
       WA%CW_RED=>WA%CPTWFP
       WA%CPROJ_RED=>WA%CPROJ
    ENDIF

  END SUBROUTINE SETWAVA


!************************* SUBROUTINE DELWAVA **************************
!
!>  destroy storage for a wavefunctionarray WA
!
!***********************************************************************

      SUBROUTINE DELWAVA(WA)

!! USE moffload

      USE prec
      IMPLICIT NONE
      TYPE (wavefuna) WA

      IF (.NOT. ASSOCIATED(WA%CPTWFP) .OR. .NOT. ASSOCIATED(WA%FERWE) .OR. &
          .NOT. ASSOCIATED(WA%CELEN)  .OR. .NOT. ASSOCIATED(WA%AUX)) THEN
         CALL vtutor%bug("internal error in DELWAVA: not all enities are associated, try DELWAVA_PROJ", "wave_high.F", 1749)
      ENDIF

!$ACC EXIT DATA DELETE(WA%CPTWFP,WA%CPROJ,WA%CW_RED,WA%CPROJ_RED) IF(OFFLOAD_ON.AND.ACC_IS_PRESENT(WA,1)) ASYNC(ACC_ASYNC_Q)
      DEALLOCATE(WA%CPTWFP,WA%CPROJ,WA%CELEN, WA%FERWE, WA%AUX)
      NULLIFY(WA%CW_RED,WA%CPROJ_RED)

      IF (ASSOCIATED(WA%CR)) THEN
!$ACC EXIT DATA DELETE(WA%CR) IF(OFFLOAD_ON.AND.ACC_IS_PRESENT(WA,1)) ASYNC(ACC_ASYNC_Q)
         DEALLOCATE(WA%CR)
      ENDIF

# 1767

      NULLIFY(WA%WDES1)
!!$ACC EXIT DATA DELETE(WA) IF(OFFLOAD_ON.AND.ACC_IS_PRESENT(WA,1)) ASYNC(ACC_ASYNC_Q)
      END SUBROUTINE DELWAVA


!************************* SUBROUTINE DELWAVA_PROJ *********************
!
!>  destroy storage for a wavefunctionarray WA
!>  wavefunction character  only
!
!***********************************************************************

      SUBROUTINE DELWAVA_PROJ(WA)

!! USE moffload

      USE prec
      IMPLICIT NONE
      TYPE (wavefuna) WA

      IF (ASSOCIATED(WA%CPTWFP) .OR. ASSOCIATED(WA%FERWE) .OR. &
          ASSOCIATED(WA%CELEN) .OR. ASSOCIATED(WA%AUX)) THEN
         CALL vtutor%bug("internal error in DELWAVA_PROJ: enities are associated, try DELWAVA", "wave_high.F", 1790)
      ENDIF
!$ACC EXIT DATA DELETE(WA%CPROJ,WA%CPROJ_RED) IF(OFFLOAD_ON.AND.ACC_IS_PRESENT(WA,1)) ASYNC(ACC_ASYNC_Q)
      DEALLOCATE(WA%CPROJ)
      NULLIFY(WA%CPROJ_RED)
# 1801

      NULLIFY(WA%WDES1)
!!$ACC EXIT DATA DELETE(WA) IF(OFFLOAD_ON.AND.ACC_IS_PRESENT(WA,1)) ASYNC(ACC_ASYNC_Q)
    END SUBROUTINE DELWAVA_PROJ


!************************* FUNCTION PCW ********************************
!
!> index the wavefunction array or character array in WA
!
!***********************************************************************

  FUNCTION PCW( WA, N1, N2)
    INTEGER N1, N2
    TYPE (wavefuna) WA
    COMPLEX(q), POINTER :: PCW(:)

    IF (N1<=0 .OR. N1> WA%FIRST_DIM) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1) // " " // &
          str(WA%FIRST_DIM), "wave_high.F", 1820)
    ENDIF
    IF (N1+ (N2-1)*WA%FIRST_DIM<=0 .OR. N1+ (N2-1)*WA%FIRST_DIM> SIZE(WA%CPTWFP,2)) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1+(N2-1)) // " " // &
          str(SIZE(WA%CPTWFP,2)), "wave_high.F", 1824)
    ENDIF

    PCW=>WA%CPTWFP(:,N1+ (N2-1)*WA%FIRST_DIM)
  END FUNCTION PCW

  FUNCTION PCPROJ( WA, N1, N2)
    INTEGER N1, N2
    TYPE (wavefuna) WA
    COMPLEX(q), POINTER :: PCPROJ(:)

    IF (N1<=0 .OR. N1> WA%FIRST_DIM) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1) // " " // &
          str(WA%FIRST_DIM), "wave_high.F", 1837)
    ENDIF
    IF (N1+ (N2-1)*WA%FIRST_DIM<=0 .OR. N1+ (N2-1)*WA%FIRST_DIM> SIZE(WA%CPTWFP,2)) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1+(N2-1)) // " " // &
          str(SIZE(WA%CPTWFP,2)), "wave_high.F", 1841)
    ENDIF

    PCPROJ=>WA%CPROJ(:,N1+ (N2-1)*WA%FIRST_DIM)
  END FUNCTION PCPROJ


!************************* SUBROUTINE ORTHON ***************************
!
!> orthogonalize a wavefunction W1 to all other bands
!> including the current band
!>
!> the subroutine uses BLAS 3 calls,
!
!***********************************************************************

  SUBROUTINE ORTHON(NK, W, W1, CQIJ, ISP)
    IMPLICIT NONE

    INTEGER NK
    TYPE (wavespin)   W
    TYPE (wavefun1)   W1
    LOGICAL LOVERL
    COMPLEX(q) CQIJ(:,:,:,:)
    INTEGER ISP
! local
    COMPLEX(q) :: CPRO(W%WDES%NBANDS),CWORK(W%WDES%NPRO)
    REAL(q) :: WFMAG
    INTEGER :: I


    IF (W%WDES%COMM_KIN%NCPU /= W%WDES%COMM_INB%NCPU) THEN
       CALL vtutor%bug("internal error: ORTHON does not support band-par.", "wave_high.F", 1873)
    ENDIF

    IF (W1%WDES1%LOVERL) THEN
       CALL OVERL1(W1%WDES1, SIZE(CQIJ,1),CQIJ(1,1,1,ISP),CQIJ(1,1,1,ISP), 0.0_q, W1%CPROJ(1),CWORK(1))
    ENDIF

    CALL ZGEMV( 'C' ,  W%WDES%NPLWKP(NK) , W%WDES%NBANDS ,(1._q,0._q) , W%CPTWFP(1,1,NK,ISP), &
                        W%WDES%NRPLWV, W1%CPTWFP(1) , 1 , (0._q,0._q) ,  CPRO(1), 1)

    IF (W1%WDES1%LOVERL) THEN
       IF (W%WDES%NPRO /= 0) &
            CALL ZGEMV( 'C' ,  W%WDES%NPRO , W%WDES%NBANDS ,(1._q,0._q) , W%CPROJ(1,1,NK,ISP) , &
            W%WDES%NPROD, CWORK(1), 1 , (1._q,0._q) ,  CPRO(1), 1)
    ENDIF

    CALL M_sum_z(W%WDES%COMM_KIN, CPRO(1), W%WDES%NBANDS)

    CALL ZGEMM( 'N', 'N' ,  W%WDES%NPLWKP(NK) , 1 , W%WDES%NBANDS , -(1._q,0._q) , &
         W%CPTWFP(1,1,NK,ISP),  W%WDES%NRPLWV , CPRO(1) , W%WDES%NBANDS , &
         (1._q,0._q) , W1%CPTWFP(1) ,  W%WDES%NRPLWV )

    IF (W%WDES%NPRO /= 0) &
         CALL ZGEMM( 'N', 'N' ,  W%WDES%NPRO , 1 , W%WDES%NBANDS  , -(1._q,0._q) , &
         W%CPROJ(1,1,NK,ISP) ,  W%WDES%NPROD , CPRO(1) , W%WDES%NBANDS , &
         (1._q,0._q) , W1%CPROJ(1) ,  W%WDES%NPROD  )

  END SUBROUTINE ORTHON

  SUBROUTINE ORTHON1P(NK, W, W1, CQIJ, ISP, NB)
    IMPLICIT NONE

    INTEGER NK
    TYPE (wavespin)   W
    TYPE (wavefun1)   W1(:)
    LOGICAL LOVERL
    COMPLEX(q) CQIJ(:,:,:,:)
    INTEGER ISP,NB(:)
! local
    COMPLEX(q) :: CPRO(W%WDES%NBANDS,SIZE(W1)),CWORK(W%WDES%NPRO,SIZE(W1))
    COMPLEX(q),ALLOCATABLE :: CW1(:,:),CPROJ1(:,:)
    REAL(q) :: WFMAG
    INTEGER :: I,BW

    BW=128


    IF (W%WDES%COMM_KIN%NCPU /= W%WDES%COMM_INB%NCPU) THEN
       CALL vtutor%bug("internal error: ORTHON does not support band-par.", "wave_high.F", 1921)
    ENDIF


    ALLOCATE(CW1(SIZE(W1(1)%CPTWFP,1),SIZE(W1)))
    DO I=1,SIZE(W1)
       CW1(:,I)=W1(I)%CPTWFP
    ENDDO

    DO I=1,SIZE(W1)
       IF (W1(I)%WDES1%LOVERL)THEN
          CALL OVERL1(W1(I)%WDES1, SIZE(CQIJ,1),CQIJ(1,1,1,ISP),CQIJ(1,1,1,ISP), 0.0_q, W1(I)%CPROJ(1),CWORK(1,I))
       ELSE
          CWORK(:,I)=0.0_q
       ENDIF
    ENDDO

    CALL ZGEMM( 'C' , 'N', W%WDES%NBANDS, size(w1),  W%WDES%NPLWKP(NK), &
           (1._q,0._q) , W%CPTWFP(1,1,NK,ISP),  W%WDES%NRPLWV, CW1(1,1) ,  W%WDES%NRPLWV, (0._q,0._q) ,  CPRO(1,1), W%WDES%NBANDS)


    IF (W%WDES%NPRO /= 0) &
         CALL ZGEMM( 'C' , 'N', W%WDES%NBANDS , size(w1), W%WDES%NPRO , (1._q,0._q) ,&
           W%CPROJ(1,1,NK,ISP) , W%WDES%NPROD, CWORK(1,1), W%WDES%NPRO, (1._q,0._q),  CPRO(1,1), W%WDES%NBANDS)

    CALL M_sum_z(W%WDES%COMM_KIN, CPRO(1,1), W%WDES%NBANDS*SIZE(W1))

    DO I=1,SIZE(W1)
       CPRO(NB(I),I)=0.d0
       IF(NB(I)-BW>=1)CPRO(1:NB(I)-BW,I)=0.d0
       IF(NB(I)+BW<=SIZE(CPRO,1))CPRO(NB(I)+BW:SIZE(CPRO,1),I)=0.d0
    ENDDO

    CALL ZGEMM( 'N', 'N' ,  W%WDES%NPLWKP(NK) , SIZE(W1) , W%WDES%NBANDS , -(1._q,0._q) , &
         W%CPTWFP(1,1,NK,ISP),  W%WDES%NRPLWV , CPRO(1,1) , W%WDES%NBANDS , &
         (1._q,0._q) , CW1(1,1) ,  W%WDES%NRPLWV )

    IF (W%WDES%NPRO /= 0)THEN
       DO I=1,SIZE(W1)
          CWORK(:,i)=W1(i)%CPROJ(:)
       ENDDO
       CALL ZGEMM( 'N', 'N' ,  W%WDES%NPRO , SIZE(W1) , W%WDES%NBANDS  , -(1._q,0._q) , &
            W%CPROJ(1,1,NK,ISP) ,  W%WDES%NPROD , CPRO(1,1) , W%WDES%NBANDS , &
            (1._q,0._q) , CWORK(1,1) ,  W%WDES%NPROD  )
    ENDIF

    DO I=1,SIZE(W1)
       W1(I)%CPTWFP=CW1(:,I)
       W1(I)%CPROJ=CWORK(:,I)
    ENDDO

    DEALLOCATE(CW1)

  END SUBROUTINE ORTHON1P


!************************* SUBROUTINE CNORMN ***************************
!
!> this subroutine normalises a wavefunction
!>
!> subroutine is not important for performance
!
!***********************************************************************

  SUBROUTINE CNORMN(W, CQIJ, ISP, WSCAL)

!! USE moffload

    IMPLICIT NONE
    TYPE (wavefun1) W

    COMPLEX(q) CQIJ(:,:,:,:)
    INTEGER ISP
    REAL(q) WSCAL
! local
    COMPLEX(q)    CP
    REAL(q) WFMAG
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NT, NIS, NI, LMMAXC, NPRO2, NPRO2_

    COMPLEX(q), EXTERNAL ::  ZDOTC
!$ACC ROUTINE(ECCP_NL) VECTOR

    WFMAG=ZDOTC(W%WDES1%NPL,W%CPTWFP(1),1,W%CPTWFP(1),1)
!=======================================================================
! if necessary caclulate <w| P |w>
!=======================================================================
    IF (W%WDES1%LOVERL) THEN
       CP  =0

!$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(W,CQIJ) REDUCTION(+:CP) &
!$ACC PRIVATE(NPRO,NPRO_,NPRO2,NPRO2_,NT,LMMAXC) DEFAULT(none) 
       spinor: DO ISPINOR=0,W%WDES1%NRSPINORS-1
          DO ISPINOR_=0,W%WDES1%NRSPINORS-1

        NPRO =ISPINOR *(W%WDES1%NPRO/2)
        NPRO_=ISPINOR_*(W%WDES1%NPRO/2)

             DO NI=1,W%WDES1%NIONS
!!           NPRO =ISPINOR *(W%WDES1%NPRO/2)
!!           NPRO_=ISPINOR_*(W%WDES1%NPRO/2)
                NT=W%WDES1%ITYP(NI)
                LMMAXC=W%WDES1%LMMAX(NT)
                IF (LMMAXC==0) CYCLE

                NPRO2 =W%WDES1%LMBASE(NI)+NPRO
                NPRO2_=W%WDES1%LMBASE(NI)+NPRO_
                CALL ECCP_NL(SIZE(CQIJ,1),LMMAXC,CQIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),W%CPROJ(NPRO2_+1),W%CPROJ(NPRO2+1),CP)
             ENDDO
          ENDDO
       ENDDO spinor

!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
       WFMAG=WFMAG+CP
    ENDIF

    CALL M_sum_d(W%WDES1%COMM_INB, WFMAG, 1)

!-----check that it is non-(0._q,0._q)
    IF(WFMAG<=0) THEN
!=======================================================================
! if it is smaller (0._q,0._q) write a warning
!=======================================================================

       IF (W%WDES1%COMM_INB%NODE_ME == W%WDES1%COMM_INB%IONODE) THEN

          WRITE(*,*)'WARNING: CNORMN: search vector ill defined'

       ENDIF

       WSCAL= -1._q/SQRT(-WFMAG)
    ELSE
       WSCAL= 1._q/SQRT(WFMAG)
    ENDIF
    CALL ZDSCAL( W%WDES1%NPL ,WSCAL,W%CPTWFP(1),1)
    CALL ZDSCAL( W%WDES1%NPRO,WSCAL,W%CPROJ(1),1)

  END SUBROUTINE CNORMN


!************************* SUBROUTINE CNORMN_REAL **********************
!
!> performs operations on real space part of wavefunction
!> after a call to CPROJCN and CNORMN
!
!***********************************************************************

  SUBROUTINE CNORMN_REAL(W, W1, ISP, WSCAL, CSCPD )

!! USE moffload_struct_def

    IMPLICIT NONE
    TYPE (wavefun1)    W, W1

    INTEGER   ISP
    REAL(q) WSCAL
    COMPLEX(q) :: CSCPD
! local
    INTEGER ISPINOR, K, KK

!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(W,W%WDES1,W%WDES1%GRID,W%WDES1%GRID%RL,W1) PRIVATE(KK) 
    DO ISPINOR=0,W%WDES1%NRSPINORS-1
       DO K=1,W%WDES1%GRID%RL%NP
          KK=K+ISPINOR*W%WDES1%GRID%MPLWV
          W%CR(KK)=(W%CR(KK)-CSCPD*W1%CR(KK))*WSCAL
       ENDDO
    ENDDO

  END SUBROUTINE CNORMN_REAL


!************************* SUBROUTINE CNORMA ***************************
!
!> this subroutine calculates the norm of a wavefunction
!
!> @details @ref openmp :
!> a ZDOTC call is replaced by an explicit OMP PARALLEL loop
!> (for some reason BLAS1 call thread really badly).
!> The nested loop over \"types\" + \"ions-of-type\" is replaced by
!> a loop over \"all ions\", that is distributed over all available
!> threads.
!
!***********************************************************************

  SUBROUTINE CNORMA(W, CQIJ, ISP, WSCAL)

!! USE moffload

    IMPLICIT NONE

    TYPE (wavefun1) W
    COMPLEX(q) :: CQIJ(:,:,:,:)
    INTEGER :: ISP
    REAL(q) :: WSCAL
! local
    COMPLEX(q)    :: CP, CTMP
    REAL(q) :: WFMAG
    INTEGER :: ISPINOR, ISPINOR_, NPRO, NPRO_, NT, NI, LMMAXC, LMDIM

    REAL(q), EXTERNAL ::  DDOT

!$  INTEGER I
!$ACC ROUTINE(ECCP_NL) VECTOR

    


!      WFMAG=ZDOTC(W%WDES1%NPL,W%CPTWFP(1),1,W%CPTWFP(1),1)
       WFMAG=DDOT(W%WDES1%NPL*2,W%CPTWFP(1),1,W%CPTWFP(1),1)
# 2136

!=======================================================================
! if necessary caclulate <w| P |w>
!=======================================================================
    IF (W%WDES1%LOVERL) THEN

       CP=0
       LMDIM=SIZE(CQIJ,1)

       

# 2151

!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) DEFAULT(NONE) &
!$OMP PRIVATE(ISPINOR,ISPINOR_,NI,NT,LMMAXC,NPRO,NPRO_,CTMP) &
!$OMP SHARED(W,CQIJ,LMDIM,ISP) REDUCTION(+:CP)

       spinor: DO ISPINOR=0,W%WDES1%NRSPINORS-1
          DO ISPINOR_=0,W%WDES1%NRSPINORS-1
             DO NI=1,W%WDES1%NIONS
                NT=W%WDES1%ITYP(NI)
                LMMAXC=W%WDES1%LMMAX(NT)
                IF (LMMAXC==0) CYCLE
                NPRO =ISPINOR *(W%WDES1%NPRO/2)+W%WDES1%LMBASE(NI)
                NPRO_=ISPINOR_*(W%WDES1%NPRO/2)+W%WDES1%LMBASE(NI)

                CTMP=0; CALL ECCP_NL(LMDIM,LMMAXC,CQIJ(1,1,NI,ISP+ISPINOR_+2*ISPINOR),W%CPROJ(NPRO_+1),W%CPROJ(NPRO+1),CTMP)
                CP=CP+CTMP
             ENDDO
          ENDDO
       ENDDO spinor
# 2173

!$OMP END PARALLEL DO

       

       WFMAG=WFMAG+CP
    ENDIF
!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
    CALL M_sum_d(W%WDES1%COMM_INB, WFMAG, 1)

!-----check that it is non-(0._q,0._q)
    IF(WFMAG<=0) THEN
!=======================================================================
! if it is smaller (0._q,0._q) write a warning
!=======================================================================

       IF (W%WDES1%COMM_INB%NODE_ME == W%WDES1%COMM_INB%IONODE) THEN

          WRITE(*,*)'WARNING: CNORMN: search vector ill defined'

       ENDIF

       WSCAL= -1._q/SQRT(-WFMAG)
    ELSE
       WSCAL= 1._q/SQRT(WFMAG)
    ENDIF

    

  END SUBROUTINE CNORMA


!************************* SUBROUTINE CINPROD **************************
!
!> this subroutine calculates the inproduct between two wavefunctions
!> ~~~
!>  <W2 | S | W1>      = W2(G)*  W1(G) +  \sum_i W2_i*Q_ij W1_j
!> ~~~
!***********************************************************************

  SUBROUTINE CINPROD(W1,W2,CQIJ,CWFMAG)
    IMPLICIT NONE

    TYPE (wavefun1)    W1
    TYPE (wavefun1)    W2

    COMPLEX(q)   CQIJ(:,:,:,:)
    COMPLEX(q) CWFMAG
! local
    COMPLEX(q)      CP
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NT, NIS, NI, LMMAXC
    COMPLEX(q), EXTERNAL :: ZDOTC
    REAL(q), EXTERNAL ::  DDOT

    CWFMAG=ZDOTC(W1%WDES1%NPL,W2%CPTWFP(1),1,W1%CPTWFP(1),1)
    CP  =0
!=======================================================================
! if necessary caclulate <w| P |w>
!=======================================================================
    IF (W1%WDES1%LOVERL) THEN
       NPRO=0

       spinor: DO ISPINOR=0,W1%WDES1%NRSPINORS-1
          DO ISPINOR_=0,W1%WDES1%NRSPINORS-1

             NPRO =ISPINOR *(W1%WDES1%NPRO/2)
             NPRO_=ISPINOR_*(W1%WDES1%NPRO/2)

             NIS =1
             DO NT=1,W1%WDES1%NTYP
                LMMAXC=W1%WDES1%LMMAX(NT)
                IF (LMMAXC==0) GOTO 230

                DO NI=NIS,W1%WDES1%NITYP(NT)+NIS-1
                   CALL ECCP_NL(SIZE(CQIJ,1),LMMAXC,CQIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),W1%CPROJ(NPRO_+1),W2%CPROJ(NPRO+1),CP)
                   NPRO = LMMAXC+NPRO
                   NPRO_= LMMAXC+NPRO_
                ENDDO
230             NIS = NIS+W1%WDES1%NITYP(NT)
             ENDDO
          ENDDO
       ENDDO spinor

       CWFMAG=CWFMAG+CP
    ENDIF

    CALL M_sum_d(W1%WDES1%COMM_INB, CWFMAG, 2)
  END SUBROUTINE CINPROD


!************************* SUBROUTINE PROJCN ***************************
!
!> this subroutine projects out from (1._q,0._q) wavefunction
!> CF another wavefunction  CPRO
!>
!> subroutine is not important for performance
!
!***********************************************************************

  SUBROUTINE PROJCN(W1, W2, CQIJ, ISP, CSCPD)

!! USE moffload

    IMPLICIT NONE

    TYPE (wavefun1)    W1,W2
    COMPLEX(q)    CADD
    COMPLEX(q) CQIJ(:,:,:,:)
    INTEGER :: ISP
    COMPLEX(q) :: CSCPD
! local
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NT, NIS, NI, LMMAXC
    COMPLEX(q), EXTERNAL ::  ZDOTC
!$ACC ROUTINE(ECCP_NL) VECTOR

    CSCPD= (ZDOTC(W1%WDES1%NPL,W2%CPTWFP(1),1,W1%CPTWFP(1),1))
!=======================================================================
! if necessary caclulate <p| P |w>
!=======================================================================
    IF (W1%WDES1%LOVERL) THEN
       CADD=0
!$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(W1,W1%WDES1,W2,CQIJ) REDUCTION(+:CADD) &
!$ACC PRIVATE(NPRO,NPRO_,NT,LMMAXC) 
       spinor: DO ISPINOR=0,W1%WDES1%NRSPINORS-1
          DO ISPINOR_=0,W1%WDES1%NRSPINORS-1
             DO NI=1,W1%WDES1%NIONS
                   NT=W1%WDES1%ITYP(NI)
                   LMMAXC=W1%WDES1%LMMAX(NT)
                   IF (LMMAXC==0) CYCLE

                   NPRO =ISPINOR *(W1%WDES1%NPRO/2)+W1%WDES1%LMBASE(NI)
                   NPRO_=ISPINOR_*(W1%WDES1%NPRO/2)+W1%WDES1%LMBASE(NI)
                   CALL ECCP_NL(SIZE(CQIJ,1),LMMAXC,CQIJ(1,1,NI,ISP+ISPINOR_+2*ISPINOR),W1%CPROJ(NPRO_+1),W2%CPROJ(NPRO+1),CADD)
             ENDDO
          ENDDO
       ENDDO spinor

!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
       CSCPD=(CSCPD+CADD)
    ENDIF
!=======================================================================
! performe orthogonalisations
!=======================================================================
    CALL M_sum_z(W1%WDES1%COMM_INB, CSCPD, 1)

    CALL ZAXPY(W1%WDES1%NPL ,-CSCPD,W2%CPTWFP(1)   ,1,W1%CPTWFP(1)   ,1)
    IF (W1%WDES1%LGAMMA) THEN
       CALL DAXPY(W1%WDES1%NPRO,-CSCPD,W2%CPROJ(1),1,W1%CPROJ(1),1)
    ELSE
       CALL ZAXPY(W1%WDES1%NPRO,-CSCPD,W2%CPROJ(1),1,W1%CPROJ(1),1)
    ENDIF

    RETURN
  END SUBROUTINE PROJCN


!************************ SUBROUTINE W1_GATHER ************************
!
!> This subroutine gathers a set of wavefunctions starting
!> from band NB1 until NB2 to all nodes
!
!**********************************************************************

  SUBROUTINE W1_GATHER( W, NB1, NB2, ISP, W1)
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1(:)  ! array into which the merge is performed

! local
    INTEGER :: NN, N, NLOC, NCPU

    NCPU=W%WDES%NB_PAR

    DO N=NB1,NB2
       NN=(N-NB1)*NCPU+W%WDES%NB_LOW
       CALL W1_COPY( ELEMENT( W, W1(NN)%WDES1, N, ISP), W1(NN) )
       CALL FFTWAV_W1( W1(NN))
    ENDDO

    NLOC=(NB2-NB1+1)*W%WDES%NB_PAR

! distribute W1 to all nodes

    IF (W%WDES%DO_REDIS) THEN
       DO NN=1,NLOC
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CPTWFP(1), &
               SIZE(W1(NN)%CPTWFP),MOD(NN-1,NCPU)+1)
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CR(1), &
               W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,MOD(NN-1,NCPU)+1)

          IF (W%WDES%LOVERL) THEN

             CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CPROJ(1), &
                  W%WDES%NPROD,MOD(NN-1,NCPU)+1)
# 2372

          ENDIF
       ENDDO


      CALL M_barrier( W%WDES%COMM_INTER )

    ENDIF

  END SUBROUTINE W1_GATHER


!************************ SUBROUTINE W1_GATHER_N **********************
!
!> This subroutine gathers a set of wavefunctions starting
!> from band NB1 until NB2 to all nodes
!>
!> compared to the routine only up to NLOC bands are collected
!
!**********************************************************************

  SUBROUTINE W1_GATHER_N( W, NB1, NB2, ISP, W1, NLOC)
!> wavefunction
    TYPE (wavespin) W
!> starting band
    INTEGER:: NB1
!> final band
    INTEGER:: NB2
!> spin
    INTEGER:: ISP
!> array into which the merge is performed
    TYPE (wavefun1):: W1(:)
!> total number of bands to be collected
    INTEGER:: NLOC

! local
    INTEGER :: NN, N, NCPU

    

    NCPU=W%WDES%NB_PAR

    DO N=NB1,NB2
       NN=(N-NB1)*NCPU+W%WDES%NB_LOW
       IF (NN>NLOC) EXIT
       CALL W1_COPY( ELEMENT( W, W1(NN)%WDES1, N, ISP), W1(NN) )
       CALL FFTWAV_W1( W1(NN))
    ENDDO

! distribute W1 to all nodes

    IF (W%WDES%DO_REDIS) THEN
       DO NN=1,NLOC
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CPTWFP(1), &
               SIZE(W1(NN)%CPTWFP),MOD(NN-1,NCPU)+1)
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CR(1), &
               W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,MOD(NN-1,NCPU)+1)

          IF (W%WDES%LOVERL) THEN

             CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CPROJ(1), &
                  W%WDES%NPROD,MOD(NN-1,NCPU)+1)
# 2437

          ENDIF
       ENDDO


      CALL M_barrier( W%WDES%COMM_INTER )


    ENDIF

    

  END SUBROUTINE W1_GATHER_N


!************************ SUBROUTINE W1_GATHER_GLB ********************
!
!> This subroutine gathers a set of wavefunctions starting
!> from band NB1 until NB2 to all nodes
!>
!> compared to the previous version the global instead of local band
!> indices are supplied
!
!**********************************************************************


  SUBROUTINE W1_GATHER_GLB( W, NB1, NB2, ISP, W1)
# 2466

    IMPLICIT NONE
!> wavefunction
    TYPE (wavespin) W
!> starting band
    INTEGER:: NB1
!> final band
    INTEGER:: NB2
!> spin
    INTEGER:: ISP
!> array into which the merge is performed
    TYPE (wavefun1):: W1(:)

! local
    INTEGER :: N_INTO_TOT, N, NCPU
    INTEGER :: ierror
# 2484


    

    NCPU=W%WDES%NB_PAR
    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          CALL W1_COPY( ELEMENT( W, W1(N_INTO_TOT-NB1+1)%WDES1, N, ISP), W1(N_INTO_TOT-NB1+1) )
          CALL FFTWAV_W1( W1(N_INTO_TOT-NB1+1))
       ENDIF
    ENDDO

! distribute W1 to all nodes

    IF (W%WDES%DO_REDIS) THEN
# 2502

       DO N_INTO_TOT=NB1,NB2
          N=N_INTO_TOT-NB1+1
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N)%CPTWFP(1), &
               SIZE(W1(N)%CPTWFP),MOD(N_INTO_TOT-1,NCPU)+1)
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N)%CR(1), &
               W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,MOD(N_INTO_TOT-1,NCPU)+1)

          IF (W%WDES%LOVERL) THEN

              CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N)%CPROJ(1), &
                   W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
# 2517

          ENDIF
       ENDDO
# 2522



      CALL M_barrier( W%WDES%COMM_INTER )


    ENDIF

    

  END SUBROUTINE W1_GATHER_GLB

!************************ SUBROUTINE W1_GATHER_GLB ********************
!
! This subroutine gathers a set of wavefunctions starting
! from band NB1 until NB2 to all nodes
! compared to the previous version the global instead of local band
! indices are supplied
!
!**********************************************************************


  SUBROUTINE W1_GATHER_GLB_( W, NB1, NB2, ISP, W1)
    IMPLICIT NONE
!> wavefunction
    TYPE (wavespin) W
!> starting band
    INTEGER:: NB1
!> final band
    INTEGER:: NB2
!> spin
    INTEGER:: ISP
!> array into which the merge is performed
    TYPE (wavefun1):: W1(:)

! local
    INTEGER :: N_INTO_TOT, N, NB_LOCAL

    

    DO N=NB1,NB2
       N_INTO_TOT=N-NB1+1
       IF (MOD(N-1,W%WDES%NB_PAR)+1==W%WDES%NB_LOW) THEN
          NB_LOCAL=1+(N-1)/W%WDES%NB_PAR
          CALL W1_COPY(ELEMENT(W,W1(N_INTO_TOT)%WDES1,NB_LOCAL,ISP),W1(N_INTO_TOT))
          CALL FFTWAV_W1(W1(N_INTO_TOT))
       ENDIF
    ENDDO

! distribute W1 to all nodes

    IF (W%WDES%DO_REDIS) THEN
       DO N=NB1,NB2
          N_INTO_TOT=N-NB1+1
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N_INTO_TOT)%CPTWFP(1), &
               SIZE(W1(N_INTO_TOT)%CPTWFP),MOD(N-1,W%WDES%NB_PAR)+1)
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N_INTO_TOT)%CR(1), &
               W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,MOD(N-1,W%WDES%NB_PAR)+1)

          IF (W%WDES%LOVERL) THEN

             CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N_INTO_TOT)%CPROJ(1), &
                  W%WDES%NPROD,MOD(N-1,W%WDES%NB_PAR)+1)
# 2588

          ENDIF
       ENDDO

       CALL M_barrier(W%WDES%COMM_INTER)

    ENDIF

    

  END SUBROUTINE W1_GATHER_GLB_


!************************ SUBROUTINE W1_GATHER_GLB_NOFFT ********************
!> identical to W1_GATHER_GLB but does not FFT to real space
!**********************************************************************

  SUBROUTINE W1_GATHER_GLB_NOFFT( W, NB1, NB2, ISP, W1)
    IMPLICIT NONE
!> wavefunction
    TYPE (wavespin) W
!> starting band
    INTEGER:: NB1
!> final band
    INTEGER:: NB2
!> spin
    INTEGER:: ISP
!> array into which the merge is performed
    TYPE (wavefun1):: W1(:)

! local
    INTEGER :: N_INTO_TOT, N, NCPU
    INTEGER :: ierror

    

    NCPU=W%WDES%NB_PAR
    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          CALL W1_COPY( ELEMENT( W, W1(N_INTO_TOT-NB1+1)%WDES1, N, ISP), W1(N_INTO_TOT-NB1+1) )
       ENDIF
    ENDDO

! distribute W1 to all nodes

    IF (W%WDES%DO_REDIS) THEN
       DO N_INTO_TOT=NB1,NB2
          N=N_INTO_TOT-NB1+1
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N)%CPTWFP(1), &
               SIZE(W1(N)%CPTWFP),MOD(N_INTO_TOT-1,NCPU)+1)

          IF (W%WDES%LOVERL) CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N)%CPROJ(1), &
               W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
# 2645

       ENDDO


      CALL M_barrier( W%WDES%COMM_INTER )


    ENDIF

    

  END SUBROUTINE W1_GATHER_GLB_NOFFT


!************************ SUBROUTINE W1_IGATHER_GLB *******************
!
!> Gathers a set of wavefunctions starting from band NB1 until NB2
!> (global band indices) to all nodes using non-blocking bcast_from.
!>
!> In case the code is compiled with-Dshmem_bcast_buffer communication
!> will be between ranks within COMM_inter_node, i.e., only those ranks
!> within COMM_INTER that are NOT on the same physical node will talk to
!> eachother. 1 communication between ranks within COMM_INTER that
!> reside on the same node is not necessary since they access a common
!> shared memory segment.
!>
!> In case the code is NOT compiled with-Dshmem_bcast_buffer, the
!> communication will be between all ranks in COMM_INTER.
!
!**********************************************************************
# 2795

  SUBROUTINE W1_IGATHER_GLB( W, NB1, NB2, ISP, W1)
    IMPLICIT NONE
!> wavefunction
    TYPE (wavespin) W
!> starting band
    INTEGER:: NB1
!> final band
    INTEGER:: NB2
!> spin
    INTEGER:: ISP
!> array into which the merge is performed
    TYPE (wavefun1):: W1(:)

! local
    INTEGER :: NI, N, NB_LOCAL
    INTEGER :: ierror

    INTEGER :: nrequests,requests(3*(NB2-NB1+1))

    

    DO N=NB1,NB2
       NI=N-NB1+1
       IF (MOD(N-1,W%WDES%NB_PAR)+1==W%WDES%NB_LOW) THEN
          NB_LOCAL=1+(N-1)/W%WDES%NB_PAR
          CALL W1_COPY(ELEMENT(W,W1(NI)%WDES1,NB_LOCAL,ISP),W1(NI))
          CALL FFTWAV_W1(W1(NI))
       ENDIF
    ENDDO

! distribute W1 to all nodes

    IF (W%WDES%COMM_INTER%NCPU>1) THEN

       nrequests=0
       DO N=NB1,NB2
          NI=N-NB1+1

          nrequests=nrequests+1
          CALL M_ibcast_z_from(W%WDES%COMM_INTER,W1(NI)%CPTWFP(1), &
               SIZE(W1(NI)%CPTWFP),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))
          nrequests=nrequests+1
          CALL M_ibcast_z_from(W%WDES%COMM_INTER,W1(NI)%CR(1), &
               SIZE(W1(NI)%CR),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))

          IF (W%WDES%LOVERL) THEN
             nrequests=nrequests+1

             CALL M_ibcast_z_from(W%WDES%COMM_INTER,W1(NI)%CPROJ(1), &
                  SIZE(W1(NI)%CPROJ),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))
# 2849

          ENDIF
       ENDDO

       CALL M_waitall(nrequests,requests(1))

    ENDIF

    

  END SUBROUTINE W1_IGATHER_GLB


!************************ SUBROUTINE W1_REDUCE_GLB ********************
!
!> Take the sum of W1(i)%CPTWFP, W1(i)%CPROJ, for i = 1, NB2-NB1+1 over all
!> ranks in COMM_inter_node, using reduce_to.
!>
!> The result is reduced onto the 1-rank in COMM_inter_node that would
!> normally own (part of) the band with the global index N = NB1, ...,NB2.
!
!**********************************************************************

  SUBROUTINE W1_REDUCE_GLB(WDES1,W1,NB1,NB2)
    IMPLICIT NONE
    TYPE (wavedes1):: WDES1
!> array into which the merge is performed
    TYPE (wavefun1):: W1(:)
!> starting band
    INTEGER:: NB1
!> final band
    INTEGER:: NB2

! local
    INTEGER :: N_INTO_TOT, N
    INTEGER :: ierror

    INTEGER :: IDO(NB2-NB1+1)

    

    DO N=NB1,NB2
       N_INTO_TOT=N-NB1+1
       IF (MOD(N-1,WDES1%NB_PAR)+1==WDES1%NB_LOW) THEN
          IDO(N_INTO_TOT)=WDES1%COMM_inter_node%NODE_ME
       ELSE
          IDO(N_INTO_TOT)=0
       ENDIF
    ENDDO

    CALL M_sum_i(WDES1%COMM_inter_node,IDO,NB2-NB1+1)

! reduce W1
    IF (WDES1%DO_REDIS) THEN

       DO N=NB1,NB2
          N_INTO_TOT=N-NB1+1
          IF (IDO(N_INTO_TOT)>0) THEN
             CALL M_reduce_z_to(WDES1%COMM_inter_node,W1(N_INTO_TOT)%CPTWFP(1), &
                  SIZE(W1(N_INTO_TOT)%CPTWFP),IDO(N_INTO_TOT))

             IF (WDES1%LOVERL) THEN

                CALL M_reduce_z_to(WDES1%COMM_inter_node,W1(N_INTO_TOT)%CPROJ(1), &
                     SIZE(W1(N_INTO_TOT)%CPROJ),IDO(N_INTO_TOT))
# 2917

             ENDIF
          ENDIF
       ENDDO

    ENDIF

    

  END SUBROUTINE W1_REDUCE_GLB


!************************ SUBROUTINE W1_IREDUCE_GLB *******************
!
!> Take the sum of W1(i)%CPTWFP, W1(i)%CPROJ, for i = 1, NB2-NB1+1 over all
!> ranks in COMM_INTER, using non-blocking reduce_to.
!>
!> The result is reduced onto the 1-rank in COMM_INTER that would
!> normally own (part of) the band with the global index N = NB1, ...,NB2.
!
!**********************************************************************

  SUBROUTINE W1_IREDUCE_GLB(WDES1,W1,NB1,NB2)
    IMPLICIT NONE
    TYPE (wavedes1) :: WDES1
!> array into which the merge is performed
    TYPE (wavefun1):: W1(:)
!> starting band
    INTEGER:: NB1
!> final band
    INTEGER:: NB2

! local
    INTEGER :: NI, N
    INTEGER :: requests(2*(NB2-NB1+1))
    INTEGER :: nrequests

    

! reduce W1

    IF (WDES1%DO_REDIS) THEN

       nrequests=0
       DO N=NB1,NB2
          NI=N-NB1+1

          nrequests=nrequests+1
          CALL M_ireduce_z_to(WDES1%COMM_INTER,W1(NI)%CPTWFP(1), &
               SIZE(W1(NI)%CPTWFP),MOD(N-1,WDES1%NB_PAR)+1,requests(nrequests))

          IF (WDES1%LOVERL) THEN
             nrequests=nrequests+1

             CALL M_ireduce_z_to(WDES1%COMM_INTER,W1(NI)%CPROJ(1), &
                  SIZE(W1(NI)%CPROJ),MOD(N-1,WDES1%NB_PAR)+1,requests(nrequests))
# 2976

          ENDIF
       ENDDO

       CALL M_waitall(nrequests,requests(1))

    ENDIF

    

  END SUBROUTINE W1_IREDUCE_GLB


!************************ SUBROUTINE W1_GATHER_GLB_NOCR ***************
!
!> This subroutine gathers a set of wavefunctions starting
!> from band NB1 until NB2 to all nodes
!>
!> compared to the previous version the global instead of local band
!> indices are supplied and the real space part is not
!> communicated
!
!**********************************************************************

  SUBROUTINE W1_GATHER_GLB_NOCR( W, NB1, NB2, ISP, W1)
    IMPLICIT NONE
!> wavefunction
    TYPE (wavespin) W
!> starting band
    INTEGER:: NB1
!> final band
    INTEGER:: NB2
!> spin
    INTEGER:: ISP
!> array into which the merge is performed
    TYPE (wavefun1):: W1(:)

! local
    INTEGER :: N_INTO_TOT, N, NCPU

    

    NCPU=W%WDES%NB_PAR
    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          CALL W1_COPY( ELEMENT( W, W1(N_INTO_TOT-NB1+1)%WDES1, N, ISP), W1(N_INTO_TOT-NB1+1) )
       ENDIF
    ENDDO

! distribute W1 to all nodes

    IF (W%WDES%DO_REDIS) THEN
       DO N_INTO_TOT=NB1,NB2
          N=N_INTO_TOT-NB1+1
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N)%CPTWFP(1), &
               SIZE(W1(N)%CPTWFP),MOD(N_INTO_TOT-1,NCPU)+1)

          IF (W%WDES%LOVERL) THEN

             CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N)%CPROJ(1), &
                  W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
# 3041

          ENDIF
       ENDDO
    ENDIF

    

  END SUBROUTINE W1_GATHER_GLB_NOCR

# 3196


!************************ SUBROUTINE W1_GATHER_GLB_ALLK ***************
!
!> Gathers a set of wavefunctions including band NB1 to NB2 over
!> all k-points to all nodes
!
!**********************************************************************

  SUBROUTINE W1_GATHER_GLB_ALLK( W, NB1, NB2, ISP, WF)
    IMPLICIT NONE
!> wavefunction
    TYPE (wavespin) W
!> starting band
    INTEGER:: NB1
!> final band
    INTEGER:: NB2
!> spin
    INTEGER:: ISP
!> array into which the merge is performed
    TYPE (wavefun1):: WF(:,:)
! local
    TYPE (wavefun1):: WAUX
    TYPE(wavedes1), TARGET :: WDESAUX
    INTEGER :: N_INTO_TOT,N,NCPU,IK,NKPTS


    CALL SETWDES(W%WDES,WDESAUX,0)
    CALL NEWWAV(WAUX,WDESAUX,.TRUE.)
    NCPU=W%WDES%NB_PAR
    NKPTS=W%WDES%NKPTS

! W1 contains bands NB1 to NB2 for all k-points
    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          DO IK=1,NKPTS
! set wave descriptor to current k-point for FFT
             CALL SETWDES(W%WDES,WDESAUX,IK)
             CALL W1_COPY_NOCR( ELEMENT( W, WDESAUX, N, ISP), WAUX)
             CALL FFTWAV_W1(WAUX)
             CALL ZCOPY( WDESAUX%GRID%MPLWV*WDESAUX%NRSPINORS, WAUX%CR(1), 1, WF(N_INTO_TOT-NB1+1,IK)%CR(1), 1)
             IF (W%WDES%LGAMMA) THEN
                CALL DCOPY( W%WDES%NPROD, WAUX%CPROJ(1), 1,  WF(N_INTO_TOT-NB1+1,IK)%CPROJ(1), 1)
             ELSE
                CALL ZCOPY( W%WDES%NPROD, WAUX%CPROJ(1), 1,  WF(N_INTO_TOT-NB1+1,IK)%CPROJ(1), 1)
             ENDIF
          ENDDO
       ENDIF
    ENDDO
    CALL DELWAV(WAUX,.TRUE.)
! distribute W1 to all nodes

    IF (W%WDES%DO_REDIS) THEN
       DO N_INTO_TOT=NB1,NB2
          N=N_INTO_TOT-NB1+1
          DO IK=1,NKPTS
             CALL M_bcast_z_from(W%WDES%COMM_INTER,WF(N,IK)%CR(1), &
                  W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,MOD(N_INTO_TOT-1,NCPU)+1)

             IF (W%WDES%LOVERL) THEN

                CALL M_bcast_z_from(W%WDES%COMM_INTER,WF(N,IK)%CPROJ(1), &
                     W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
# 3263

             ENDIF
          ENDDO
       ENDDO


       CALL M_barrier( W%WDES%COMM_INTER )


    ENDIF

  END SUBROUTINE W1_GATHER_GLB_ALLK


!************************ SUBROUTINE W1_GATHER_ARRAY ******************
!
!> This subroutine gathers a set of orbitals
!> from band NB1 until NB2 to all nodes
!>
!> compared to previous versions the set in collected
!> into work arrays
!
!**********************************************************************

  SUBROUTINE W1_GATHER_ARRAY( W, NB1, NB2, ISP, W1, CR, CPROJ)
!> wavefunction
    TYPE (wavespin) W
!> starting band
    INTEGER:: NB1
!> final band
    INTEGER:: NB2
!> spin
    INTEGER:: ISP
!> array for FFT
    TYPE (wavefun1):: W1
!> collected real space orbitals
    COMPLEX(q)    :: CR(:,:)
!> collected projected orbitals
    COMPLEX(q)    :: CPROJ(:,:)

! local
    INTEGER :: NN, N, NLOC, NCPU

    IF (W%WDES%NRSPINORS/=1) THEN
       CALL vtutor%bug("internal error in W1_GATHER_ARRAY: at the moment spinors are not supported", "wave_high.F", 3307)
    ENDIF

    NCPU=W%WDES%NB_PAR

    DO N=NB1,NB2
       NN=(N-NB1)*NCPU+W%WDES%NB_LOW
       CALL W1_COPY( ELEMENT( W, W1%WDES1, N, ISP), W1 )
       CALL FFTWAV_W1( W1)

       IF (NN> SIZE(CR,2) .OR. NN >SIZE(CPROJ,2)) THEN
          CALL vtutor%bug("internal error in W1_GATHER_ARRAY: bound exceed " // str(NN) // " " // &
             str(SIZE(CR,2)) // " " // str(SIZE(CPROJ,2)), "wave_high.F", 3319)
       ENDIF
       IF (SIZE(CR,1)/=W1%WDES1%GRID%NPLWV) THEN
          WRITE(*,*) 'internal error in W1_GATHER_ARRAY: size mismatch ',SIZE(CR,1),W1%WDES1%GRID%NPLWV
       ENDIF

       CR(1:W1%WDES1%GRID%RL%NP, NN)=W1%CR(1:W1%WDES1%GRID%RL%NP)
! pad CR with zeros (just in case)
       CR(W1%WDES1%GRID%RL%NP+1:SIZE(CR,1), NN)=0
       CPROJ(:,NN)=W1%CPROJ(:)

! distribute W1 to all nodes

       IF (W%WDES%DO_REDIS) THEN
        DO NN=1,NCPU

          CALL M_bcast_z_from(W%WDES%COMM_INTER,CR(1,NN+(N-NB1)*NCPU), &
               SIZE(CR,1),MOD(NN-1,NCPU)+1)
# 3340

          IF (W%WDES%LOVERL) THEN

             CALL M_bcast_z_from(W%WDES%COMM_INTER,CPROJ(1,NN+(N-NB1)*NCPU), &
                  W%WDES%NPROD,MOD(NN-1,NCPU)+1)
# 3348

            ENDIF
         ENDDO
       ENDIF

      ENDDO

  END SUBROUTINE W1_GATHER_ARRAY


!************************ SUBROUTINE W1_GATHER_ARRAY_RECIPROCAL *******
!
!> This subroutine gathers a set of orbitals
!> from band NB1 until NB2 to all nodes
!>
!> this version collects the place wave coefficients and the
!> CPROJ coefficients
!
!**********************************************************************

  SUBROUTINE W1_GATHER_ARRAY_RECIPROCAL( W, NB1, NB2, ISP, W1, CG, CPROJ)

!! USE moffload_struct_def

!> wavefunction
    TYPE (wavespin) W
!> starting band
    INTEGER:: NB1
!> final band
    INTEGER:: NB2
!> spin
    INTEGER:: ISP
!> array for FFT
    TYPE (wavefun1):: W1
!> collected real space orbitals
    COMPLEX(q):: CG(:,:)
!> collected projected orbitals
    COMPLEX(q)       :: CPROJ(:,:)

    COMPLEX(q) :: CWBUFF(SIZE(CG,1))       ! send buffer
    COMPLEX(q)       :: CPROJBUFF(SIZE(CPROJ,1)) ! send buffer for projectors

! local
    INTEGER :: NN, N, NCPU

    

!$ACC ENTER DATA CREATE(CWBUFF,CPROJBUFF) 

    IF (W%WDES%NRSPINORS/=1) THEN
       CALL vtutor%bug("W1_GATHER_ARRAY_RECIPROCAL: at the moment spinors are not supported", "wave_high.F", 3398)
    ENDIF

    IF (SIZE(CPROJ,1) > SIZE(W1%CPROJ,1)) THEN
       CALL vtutor%bug("W1_GATHER_ARRAY_RECIPROCAL: CPROJ size inconsistent " // &
          str(SIZE(CPROJ,1)) // " " // str(SIZE(W1%CPROJ,1)), "wave_high.F", 3403)
    ENDIF

    NN=(NB2-NB1)*W%WDES%NB_PAR+W%WDES%NB_LOW
    IF (NN>SIZE(CG,2) .OR. NN>SIZE(CPROJ,2)) THEN
       CALL vtutor%bug("W1_GATHER_ARRAY_RECIPROCAL: bound exceeded " // str(NN) &
          // " " // str(SIZE(CG,2)) // " " // str(SIZE(CPROJ,2)), "wave_high.F", 3409)
    ENDIF

    NCPU=W%WDES%NB_PAR

! NB1 and NB2 are local indices for collecting the bands
    DO N=NB1,NB2
!NB_LOW is the off-set of the local node
!NN is index into CG array
       NN=(N-NB1)*NCPU+W%WDES%NB_LOW
!this copies orbital N (reciprocal and real (if allocated) part to W1)
!these are local indices
       CALL W1_COPY( ELEMENT( W, W1%WDES1, N, ISP), W1 )

! copy data over to return arrays: CG and CPROJ
!$ACC KERNELS PRESENT(CG,CPROJ,W1) 
       CG(1:W1%WDES1%NPL ,NN)=W1%CPTWFP(1:W1%WDES1%NPL)
       CG(W1%WDES1%NPL+1:,NN)=0.0_q                 ! pad with (0._q,0._q)

       CPROJ(:,NN)=W1%CPROJ(1:SIZE(CPROJ,1))
!$ACC END KERNELS

! distribute W1 to all nodes, CG(1,(N-NB)*NCPU)
       IF (W%WDES%DO_REDIS) THEN
! MPI_IN_PLACE was tested but found to be much slower (code has been removed)
!$ACC KERNELS PRESENT(CWBUFF,CPROJBUFF,CG,CPROJ) 
          CPROJBUFF(:)=CPROJ(:,NN)
          CWBUFF(:)   =CG(:,NN)
!$ACC END KERNELS
          CALL M_allgathero_z(W%WDES%COMM_INTER,SIZE(CWBUFF),CWBUFF(1),CG(1,(N-NB1)*NCPU+1))

          IF (W%WDES%LOVERL) THEN
# 3443

             CALL M_allgathero_z(W%WDES%COMM_INTER,SIZE(CPROJBUFF),CPROJBUFF(1),CPROJ(1,(N-NB1)*NCPU+1))

          ENDIF
       ENDIF

    ENDDO

!$ACC EXIT DATA DELETE(CWBUFF,CPROJBUFF) 

    

  END SUBROUTINE W1_GATHER_ARRAY_RECIPROCAL

!
!> A version of W1_GATHER_ARRAY_RECIPROCAL using non-blocking allgather operations
!
  SUBROUTINE W1_IGATHER_ARRAY_RECIPROCAL( W, NB1, NB2, ISP, W1, CG, CPROJ)
!> wavefunction
    TYPE (wavespin) W
!> starting band
    INTEGER:: NB1
!> final band
    INTEGER:: NB2
!> spin
    INTEGER:: ISP
!> array for FFT
    TYPE (wavefun1):: W1
!> collected real space orbitals
    COMPLEX(q):: CG(:,:)
!> collected projected orbitals
    COMPLEX(q)       :: CPROJ(:,:)

    COMPLEX(q) :: CWBUFF(SIZE(CG,1))       ! send buffer
    COMPLEX(q)       :: CPROJBUFF(SIZE(CPROJ,1)) ! send buffer for projectors

! local
    INTEGER :: NN, N, NCPU

    INTEGER :: nrequests, requests(2*(NB2-NB1+1))

    

    IF (W%WDES%NRSPINORS/=1) THEN
       CALL vtutor%bug("internal error in W1_IGATHER_ARRAY_RECIPROCAL: at the moment spinors are not supported", "wave_high.F", 3487)
    ENDIF

    IF (SIZE(CPROJ,1) > SIZE(W1%CPROJ,1)) THEN
       CALL vtutor%bug("W1_IGATHER_ARRAY_RECIPROCAL: CPROJ size inconsistent " // &
          str(SIZE(CPROJ,1)) // " " // str(SIZE(W1%CPROJ,1)), "wave_high.F", 3492)
    ENDIF

    NN=(NB2-NB1)*W%WDES%NB_PAR+W%WDES%NB_LOW
    IF (NN>SIZE(CG,2) .OR. NN>SIZE(CPROJ,2)) THEN
       CALL vtutor%bug("W1_IGATHER_ARRAY_RECIPROCAL: bound exceeded " // str(NN) &
          // " " // str(SIZE(CG,2)) // " " // str(SIZE(CPROJ,2)), "wave_high.F", 3498)
    ENDIF

    NCPU=W%WDES%NB_PAR

    nrequests=0

! NB1 and NB2 are local indices for collecting the bands
    DO N=NB1,NB2
!NB_LOW is the off-set of the local node
!NN is index into CG array
       NN=(N-NB1)*NCPU+W%WDES%NB_LOW
!this copies orbital N (reciprocal and real (if allocated) part to W1)
!these are local indices
       CALL W1_COPY( ELEMENT( W, W1%WDES1, N, ISP), W1 )

! copy data over to return arrays: CG and CPROJ
       CG(1:W1%WDES1%NPL ,NN)=W1%CPTWFP(1:W1%WDES1%NPL)
       CG(W1%WDES1%NPL+1:,NN)=0.0_q                 ! pad with (0._q,0._q)

       CPROJ(:,NN)=W1%CPROJ(1:SIZE(CPROJ,1))

! distribute W1 to all nodes, CG(1,(N-NB)*NCPU)
       IF (W%WDES%DO_REDIS) THEN
! MPI_IN_PLACE was tested but found to be much slower (code has been removed)
          CPROJBUFF(:)=CPROJ(:,NN)
          CWBUFF(:)   =CG(:,NN)

          nrequests=nrequests+1
          CALL M_iallgathero_z(W%WDES%COMM_INTER,SIZE(CWBUFF),CWBUFF(2),CG(1,(N-NB1)*NCPU+1),requests(nrequests))

          IF (W%WDES%LOVERL) THEN
             nrequests=nrequests+1
# 3533

             CALL M_iallgathero_z(W%WDES%COMM_INTER,SIZE(CPROJBUFF),CPROJBUFF(1),CPROJ(1,(N-NB1)*NCPU+1),requests(nrequests))

          ENDIF
       ENDIF

    ENDDO

    CALL M_waitall(nrequests,requests(1))

    

  END SUBROUTINE W1_IGATHER_ARRAY_RECIPROCAL

!
!> old version, uses bcast and is much slower
!
  SUBROUTINE W1_GATHER_ARRAY_RECIPROCAL_OLD( W, NB1, NB2, ISP, W1, CG, CPROJ)
    TYPE (wavespin) W        !< wavefunction
    INTEGER :: NB1           !< starting band
    INTEGER :: NB2           !< final band
    INTEGER :: ISP           !< spin
    TYPE (wavefun1):: W1     !< array for FFT
    COMPLEX(q) :: CG(:,:)    !< collected real space orbitals
    COMPLEX(q)       :: CPROJ(:,:) !< collected projected orbitals

! local
    INTEGER :: NN, N, NCPU

    

    IF (W%WDES%NRSPINORS/=1) THEN
       CALL vtutor%bug("internal error in W1_GATHER_ARRAY_RECIPROCAL: at the moment spinors are not supported", "wave_high.F", 3565)
    ENDIF

    IF (SIZE(CPROJ,1) > SIZE(W1%CPROJ,1)) THEN
       CALL vtutor%bug("W1_GATHER_ARRAY_RECIPROCAL: CPROJ size inconsistent " // &
          str(SIZE(CPROJ,1)) // " " // str(SIZE(W1%CPROJ,1)), "wave_high.F", 3570)
    ENDIF

    NN=(NB2-NB1)*W%WDES%NB_PAR+W%WDES%NB_LOW
    IF (NN>SIZE(CG,2) .OR. NN>SIZE(CPROJ,2)) THEN
       CALL vtutor%bug("W1_GATHER_ARRAY_RECIPROCAL: bound exceeded " // str(NN) &
          // " " // str(SIZE(CG,2)) // " " // str(SIZE(CPROJ,2)), "wave_high.F", 3576)
    ENDIF

    NCPU=W%WDES%NB_PAR

    DO N=NB1,NB2
!NB_LOW will be the off-set of the local node
       NN=(N-NB1)*NCPU+W%WDES%NB_LOW
!this copies orbital N (reciprocal and real (if allocated) part to W1)
!these are local indices
       CALL W1_COPY( ELEMENT( W, W1%WDES1, N, ISP), W1 )

!copy the local data to the position with global index in the strip
       CG(1:W1%WDES1%NPL, NN)=W1%CPTWFP(1:W1%WDES1%NPL)
       CG(W1%WDES1%NPL+1: SIZE(CG,1), NN)=0

!the first index CPROJ size can be smaller than the second (1._q,0._q) (less cpus/tau than total number)
!but this copies only the data for valid indices in the first array (or?)
       CPROJ(:,NN)=W1%CPROJ(1:SIZE(CPROJ,1))

! distribute CG and CPROJ to all nodes

       IF (W%WDES%DO_REDIS) THEN
          DO NN=1,NCPU
!so broadcast the data from node holding index NN to all other nodes
!I wonder how efficient this is, might be all right 10000 coefficients is 160kB of data,
!a bit on the low side but fine
             CALL M_bcast_z_from(W%WDES%COMM_INTER,CG(1,NN+(N-NB1)*NCPU), &
                  SIZE(CG,1),MOD(NN-1,NCPU)+1)

             IF (W%WDES%LOVERL) THEN
! OK, the point is that the W%WDES%NPROD array can be longer than the size of CPROJ
! since both are a multiple of NCPU but once the total number and in the second case per tau

                CALL M_bcast_z_from(W%WDES%COMM_INTER,CPROJ(1,NN+(N-NB1)*NCPU), &
                     SIZE(CPROJ,1),MOD(NN-1,NCPU)+1)
# 3615

             ENDIF
          ENDDO
       ENDIF

    ENDDO

    

  END SUBROUTINE W1_GATHER_ARRAY_RECIPROCAL_OLD


!************************ SUBROUTINE W1_GATHER_DISTR ******************
!
!> This subroutine gathers a set of wavefunctions starting
!> from band NB1 until NB2 to all nodes
!>
!> compared to the previous version the global instead of local band
!> indices are supplied
!
!**********************************************************************

  SUBROUTINE W1_GATHER_DISTR( W, NB1, NB2, ISP, W1)
    IMPLICIT NONE
!> wavefunction
    TYPE (wavespin) W
!> starting band
    INTEGER:: NB1
!> final band
    INTEGER:: NB2
!> spin
    INTEGER:: ISP
!> array into which the merge is performed
    TYPE (wavefun1):: W1(:)
! local variables
    INTEGER :: N_INTO_TOT, N, NCPU
    INTEGER :: NBMIN, NBMAX
    TYPE (wavefun1):: WTMP

    NCPU=W%WDES%NB_PAR

    CALL NEWWAV(WTMP, W1(1)%WDES1, .FALSE.)

! establish global band index interval
    NBMIN=NB1
    NBMAX=NB2

    DO N=1,NCPU
       CALL M_bcast_i_from(W%WDES%COMM_INTER, NBMIN, 1, n)
       CALL M_bcast_i_from(W%WDES%COMM_INTER, NBMAX, 1, n)
       IF (NB1<NBMIN) NBMIN=NB1
       IF (NB2>NBMAX) NBMAX=NB2
    ENDDO


    DO N_INTO_TOT=NBMIN,NBMAX

! if band resides on this node copy it to WTMP
       IF (MOD(N_INTO_TOT-W%WDES%NB_LOW,NCPU)==0) THEN
          N=(N_INTO_TOT-W%WDES%NB_LOW)/NCPU+1
          CALL W1_COPY( ELEMENT( W, WTMP%WDES1, N, ISP), WTMP )
       ENDIF


! broadcast WTMP from the node where it resides
       CALL M_bcast_z_from(W%WDES%COMM_INTER,WTMP%CPTWFP(1), &
            SIZE(WTMP%CPTWFP),MOD(N_INTO_TOT-1,NCPU)+1)

       IF (W%WDES%LOVERL) THEN

          CALL M_bcast_z_from(W%WDES%COMM_INTER,WTMP%CPROJ(1), &
               W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
# 3690

       ENDIF

! if this band is targeted to reside on this node
! copy WTMP to W1(N_INTO_TOT-NB1+1)
       IF (N_INTO_TOT>=NB1 .AND. N_INTO_TOT<=NB2) THEN
          N=N_INTO_TOT-NB1+1
          CALL W1_COPY(WTMP,W1(N))
       ENDIF
    ENDDO

! FFT to real space
    DO N=1,(NB2-NB1)+1
       CALL FFTWAV_W1(W1(N))
    ENDDO

    CALL DELWAV(WTMP, .FALSE.)

  END SUBROUTINE W1_GATHER_DISTR

!************************ SUBROUTINE W1_GATHER_KSEL_NOFFT *************
!
! This subroutine gathers a set of wavefunctions starting
! from band NB1 until NB2 and distributes the data over k in a round
! robin fashion
! the data distribution is based on an index (k-point index) supplied as
! the last argument
!
!**********************************************************************


  SUBROUTINE W1_GATHER_KSEL_NOFFT( W, NB1, NB2, ISP, W1, NODE, NK)
# 3724

    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1(:)  ! array into which the merge is performed
    INTEGER :: NODE          ! node to recieve k-point
    INTEGER :: NK  ! k-points intex

! local
    INTEGER :: NN, N, NLOC, NCPU
    INTEGER :: NCP
# 3738



    

    NCPU=W%WDES%NB_PAR
    MPLWV = SIZE(W%CPTWFP(:,1,NK,ISP))
    NPROD = W%WDES%NPROD

# 3749

    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          NCP = N
       ENDIF

       IF (MOD(NODE-1,NCPU)+1 ==  W%WDES%NB_LOW) THEN
! receive from all other nodes and local copy
          DO NN=1,NCPU
             N_INTO_TOT=(N-1)*NCPU+NN
             IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
                IF (NN /= W%WDES%NB_LOW) THEN
                   CALL M_recv_z(W%WDES%COMM_INTER, NN , &
                        W1(N_INTO_TOT-NB1+1)%CPTWFP(1), SIZE(W1(N_INTO_TOT-NB1+1)%CPTWFP))

                   IF (W%WDES%LOVERL) THEN

                      CALL M_recv_z(W%WDES%COMM_INTER, NN, &
                           W1(N_INTO_TOT-NB1+1)%CPROJ(1), W%WDES%NPROD)
# 3772

                   ENDIF
                ENDIF
             ENDIF
          ENDDO
       ELSE
          IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
             CALL M_send_z(W%WDES%COMM_INTER, MOD(NODE-1,NCPU)+1, &
                      W%CPTWFP(:,NCP,NK,ISP),SIZE(W%CPTWFP(:,NCP,NK,ISP)))

             IF (W%WDES%LOVERL) THEN

                CALL M_send_z(W%WDES%COMM_INTER, MOD(NODE-1,NCPU)+1, &
                     W%CPROJ(:,NCP,NK,ISP), W%WDES%NPROD)
# 3789

             ENDIF
          ENDIF
       ENDIF

    ENDDO
# 3797


    IF (MOD(NODE-1,NCPU)+1 ==  W%WDES%NB_LOW) THEN
!$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) PRIVATE(N_INTO_TOT,NCP) PRESENT(W1,W) 
       DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
     N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
     IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
        NCP = N
        DO NN = 1, MPLWV
           W1(N_INTO_TOT-NB1+1)%CPTWFP(NN) = W%CPTWFP(NN,NCP,NK,ISP)
        ENDDO
     ENDIF
!!     DO NN = 1, MPLWV
!!        N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
!!        IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
!!           NCP = N
!!           W1(N_INTO_TOT-NB1+1)%CPTWFP(NN) = W%CPTWFP(NN,NCP,NK,ISP)
!!        ENDIF
!!     ENDDO
       ENDDO
!$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) PRIVATE(N_INTO_TOT,NCP) PRESENT(W1,W) 
       DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
     N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
     IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
        NCP = N
        DO NN = 1, NPROD
           W1(N_INTO_TOT-NB1+1)%CPROJ(NN) = W%CPROJ(NN,NCP,NK,ISP)
        ENDDO
     ENDIF
!!     DO NN = 1, NPROD
!!        N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
!!        IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
!!           NCP = N
!!           W1(N_INTO_TOT-NB1+1)%CPROJ(NN) = W%CPROJ(NN,NCP,NK,ISP)
!!        ENDIF
!!     ENDDO
       ENDDO
    ENDIF

    

  END SUBROUTINE W1_GATHER_KSEL_NOFFT

!************************ SUBROUTINE W1_GATHER_KSEL *******************
!
!> This subroutine gathers a set of wavefunctions starting
!> from band NB1 until NB2 and distributes the data over k in a round
!> robin fashion
!>
!> the data distribution is based on an index (k-point index) supplied as
!> the last argument
!
!**********************************************************************


  SUBROUTINE W1_GATHER_KSEL( W, NB1, NB2, ISP, W1, NK)
!> wavefunction
    TYPE (wavespin) W
!> starting band
    INTEGER:: NB1
!> final band
    INTEGER:: NB2
!> spin
    INTEGER:: ISP
!> array into which the merge is performed
    TYPE (wavefun1):: W1(:)
!> k-points intex
    INTEGER:: NK

! local
    TYPE (wavefun1):: WTMP   ! temporary for FFT
    INTEGER :: NN, N, NLOC, NCPU
# 3871



    

    NCPU=W%WDES%NB_PAR

!$ACC ENTER DATA CREATE(WTMP) 
    CALL NEWWAV(WTMP, W1(1)%WDES1, .TRUE.)

    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          CALL W1_COPY( ELEMENT( W, W1(1)%WDES1, N, ISP), WTMP )
          CALL FFTWAV_W1( WTMP)
       ENDIF

# 3890

       IF (MOD(NK-1,NCPU)+1 ==  W%WDES%NB_LOW) THEN
          
! receive from all other nodes and local copy
          DO NN=1,NCPU
             N_INTO_TOT=(N-1)*NCPU+NN
             IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
                IF (NN==W%WDES%NB_LOW) THEN
                   CALL W1_COPY( WTMP, W1(N_INTO_TOT-NB1+1) )
                ELSE
!                   WRITE(*,*) W%WDES%NB_LOW, 'receive', N_INTO_TOT-NB1+1, 'from' , NN
                   CALL M_recv_z(W%WDES%COMM_INTER, NN , &
                        W1(N_INTO_TOT-NB1+1)%CR(1), W%WDES%GRID%MPLWV*W%WDES%NRSPINORS)

                   IF (W%WDES%LOVERL) THEN

                      CALL M_recv_z(W%WDES%COMM_INTER, NN, &
                           W1(N_INTO_TOT-NB1+1)%CPROJ(1), W%WDES%NPROD)
# 3911

                   ENDIF
                ENDIF
             ENDIF
          ENDDO
          
       ELSE
          IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
            
!             WRITE(*,*) W%WDES%NB_LOW,'send', N_INTO_TOT, 'to' , MOD(NK-1,NCPU)+1
             CALL M_send_z(W%WDES%COMM_INTER, MOD(NK-1,NCPU)+1, &
                  WTMP%CR(1), W%WDES%GRID%MPLWV*W%WDES%NRSPINORS)

             IF (W%WDES%LOVERL) THEN

                CALL M_send_z(W%WDES%COMM_INTER, MOD(NK-1,NCPU)+1, &
                     WTMP%CPROJ(1), W%WDES%NPROD)
# 3931

             ENDIF
            
          ENDIF
       ENDIF
# 3938


    ENDDO

    CALL DELWAV(WTMP, .TRUE.)
!$ACC EXIT DATA DELETE(WTMP) 

    

  END SUBROUTINE W1_GATHER_KSEL



!************************ SUBROUTINE W1_GATHER_KNODESEL ***************
!
!> This subroutine is similar to #W1_GATHER_KSEL,
!> but at variance with
!> it, (1._q,0._q) explicitly select the node on which you gather the wave
!> fucntion (instead of distributing k-points in round robin fashion).
!
!**********************************************************************

  SUBROUTINE W1_GATHER_KNODESEL( W, NB1, NB2, ISP, W1, NODE)
!> wavefunction
    TYPE (wavespin) W
!> starting band
    INTEGER:: NB1
!> final band
    INTEGER:: NB2
!> spin
    INTEGER:: ISP
!> array into which the merge is performed
    TYPE (wavefun1):: W1(:)
!> which node receives
    INTEGER:: NODE

! local
    TYPE (wavefun1):: WTMP   ! temporary for FFT
    INTEGER :: NN, N, NLOC, NCPU

    NCPU=W%WDES%NB_PAR

    CALL NEWWAV(WTMP, W1(1)%WDES1, .TRUE.)

    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          CALL W1_COPY( ELEMENT( W, W1(1)%WDES1, N, ISP), WTMP )
          CALL FFTWAV_W1( WTMP)
       ENDIF

       IF (NODE==W%WDES%NB_LOW) THEN
! receive from all other nodes and local copy
          DO NN=1,NCPU
             N_INTO_TOT=(N-1)*NCPU+NN
             IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
                IF (NN==W%WDES%NB_LOW) THEN
                   CALL W1_COPY( WTMP, W1(N_INTO_TOT-NB1+1) )
                ELSE
                   CALL M_recv_z(W%WDES%COMM_INTER, NN , &
                        W1(N_INTO_TOT-NB1+1)%CR(1), W%WDES%GRID%MPLWV*W%WDES%NRSPINORS)

                   IF (W%WDES%LOVERL) THEN

                      CALL M_recv_z(W%WDES%COMM_INTER, NN, &
                           W1(N_INTO_TOT-NB1+1)%CPROJ(1), W%WDES%NPROD)
# 4007

                   ENDIF
                ENDIF
             ENDIF
          ENDDO
       ELSE
          IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
             CALL M_send_z(W%WDES%COMM_INTER,NODE, &
                  WTMP%CR(1), W%WDES%GRID%MPLWV*W%WDES%NRSPINORS)

             IF (W%WDES%LOVERL) THEN

                CALL M_send_z(W%WDES%COMM_INTER,NODE, &
                     WTMP%CPROJ(1), W%WDES%NPROD)
# 4024

             ENDIF
          ENDIF
       ENDIF

    ENDDO
    CALL DELWAV(WTMP, .TRUE.)
  END SUBROUTINE W1_GATHER_KNODESEL


!************************ SUBROUTINE W1_GATHER_W1 *********************
!
!> This subroutine gathers from a wavefunction array
!> instead of W
!
!**********************************************************************

  SUBROUTINE W1_GATHER_W1( W, NB2, W1_ORIG, W1)
!> wavefunction
    TYPE (wavespin):: W
!> wavefunction array to be merged
    TYPE (wavefun1):: W1_ORIG(:)
!> array into which the merge is performed
    TYPE (wavefun1):: W1(:)
!> final band
    INTEGER:: NB2

! local
    INTEGER :: NN, N, NLOC, NCPU

    NCPU=W%WDES%NB_PAR
    IF (SIZE(W1_ORIG)< NB2) THEN
       CALL vtutor%bug("internal error in W1_GATHER_W1: W1_ORIG is not sufficiently large", "wave_high.F", 4056)
    ENDIF

    DO N=1,NB2
       NN=(N-1)*NCPU+W%WDES%NB_LOW
       CALL W1_COPY_NOCR( W1_ORIG(N), W1(NN) )
       CALL FFTWAV_W1( W1(NN))
    ENDDO

    NLOC=NB2*W%WDES%NB_PAR

! distribute W1 to all nodes

    IF (W%WDES%DO_REDIS) THEN
       DO NN=1,NLOC
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CR(1), &
               W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,MOD(NN-1,NCPU)+1)

          IF (W%WDES%LOVERL) THEN

             CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CPROJ(1), &
                  W%WDES%NPROD,MOD(NN-1,NCPU)+1)
# 4081

          ENDIF
       ENDDO


      CALL M_barrier( W%WDES%COMM_INTER )


    ENDIF


  END SUBROUTINE W1_GATHER_W1


!************************ SUBROUTINE W1_GATHER_STRIP ******************
!
!> This subroutine gathers from a wavefunction array instead of W.
!> Collected are the bands W1_ORIG(N): NB1 <= NB_LOW+(N-1)*NB_PAR <= NB2,
!> into W1(1:NB2-NB1+1)
!
!**********************************************************************

  SUBROUTINE W1_GATHER_STRIP( W, NB1, NB2, W1_ORIG, W1)
!> wavefunction
    TYPE (wavespin):: W
!> wavefunction array to be merged
    TYPE (wavefun1):: W1_ORIG(:)
!> array into which the merge is performed
    TYPE (wavefun1):: W1(:)
!> first band
    INTEGER:: NB1
!> final band
    INTEGER:: NB2

! local
    INTEGER :: NN, N, NCPU

    

    NCPU=W%WDES%NB_PAR
    IF (SIZE(W1_ORIG)*NCPU<NB2) THEN
       CALL vtutor%bug("internal error in W1_GATHER_STRIP: W1_ORIG is not sufficiently large", "wave_high.F", 4122)
    ENDIF

    IF (SIZE(W1)<NB2-NB1+1) THEN
       CALL vtutor%bug("internal error in W1_GATHER_STRIP: W1 is not sufficiently large", "wave_high.F", 4126)
    ENDIF

    DO N=1,SIZE(W1_ORIG)
       NN=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NN<NB1.OR.NN>NB2) CYCLE
       CALL W1_COPY_NOCR( W1_ORIG(N), W1(NN-NB1+1) )
       CALL FFTWAV_W1( W1(NN-NB1+1))
    ENDDO

! distribute W1 to all nodes

    IF (W%WDES%DO_REDIS) THEN
       DO NN=NB1,NB2
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN-NB1+1)%CR(1), &
               W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,MOD(NN-1,NCPU)+1)

          IF (W%WDES%LOVERL) THEN

             CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN-NB1+1)%CPROJ(1), &
                  W%WDES%NPROD,MOD(NN-1,NCPU)+1)
# 4150

          ENDIF
       ENDDO


      CALL M_barrier( W%WDES%COMM_INTER )


    ENDIF

    

  END SUBROUTINE W1_GATHER_STRIP


!************************ SUBROUTINE W1_IGATHER_STRIP *****************
!
!> Gathers a set of wavefunctions starting from band NB1 until NB2
!> (global band indices) to all nodes using non-blocking bcast_from.
!> This subroutine gathers from a wavefunction array W1_ORIG instead
!> of W.
!>
!> In case the code is compiled with-Dshmem_bcast_buffer communication
!> will be between ranks within COMM_inter_node, i.e., only those ranks
!> within COMM_INTER that are NOT on the same physical node will talk to
!> eachother. 1 communication between ranks within COMM_INTER that
!> reside on the same node is not necessary since they access a common
!> shared memory segment.
!>
!> In case the code is NOT compiled with-Dshmem_bcast_buffer, the
!> communication will be between all ranks in COMM_INTER.
!
!**********************************************************************
# 4302

  SUBROUTINE W1_IGATHER_STRIP( W, NB1, NB2, W1_ORIG, W1)
!> wavefunction
    TYPE (wavespin):: W
!> wavefunction array to be merged
    TYPE (wavefun1):: W1_ORIG(:)
!> array into which the merge is performed
    TYPE (wavefun1):: W1(:)
!> first band
    INTEGER:: NB1
!> final band
    INTEGER:: NB2

! local
    INTEGER :: NI,N,NB_LOCAL
    INTEGER :: ierror

    INTEGER :: nrequests,requests(3*(NB2-NB1+1))


    

    IF (SIZE(W1_ORIG)*W%WDES%NB_PAR<NB2) THEN
       CALL vtutor%bug("internal error in W1_IGATHER_STRIP: W1_ORIG is not sufficiently large", "wave_high.F", 4325)
    ENDIF

    IF (SIZE(W1)<NB2-NB1+1) THEN
       CALL vtutor%bug("internal error in W1_IGATHER_STRIP: W1 is not sufficiently large", "wave_high.F", 4329)
    ENDIF

    DO N=NB1,NB2
       NI=N-NB1+1
       IF (MOD(N-1,W%WDES%NB_PAR)+1==W%WDES%NB_LOW) THEN
          NB_LOCAL=1+(N-1)/W%WDES%NB_PAR
          CALL W1_COPY_NOCR(W1_ORIG(NB_LOCAL),W1(NI))
          CALL FFTWAV_W1(W1(NI))
       ENDIF
    ENDDO

! distribute W1 to all nodes

    IF (W%WDES%COMM_INTER%NCPU>1) THEN

       nrequests=0
       DO N=NB1,NB2
          NI=N-NB1+1

          nrequests=nrequests+1
          CALL M_ibcast_z_from(W%WDES%COMM_INTER,W1(NI)%CPTWFP(1), &
               SIZE(W1(NI)%CPTWFP),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))
          nrequests=nrequests+1
          CALL M_ibcast_z_from(W%WDES%COMM_INTER,W1(NI)%CR(1), &
               SIZE(W1(NI)%CR),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))

          IF (W%WDES%LOVERL) THEN
             nrequests=nrequests+1

             CALL M_ibcast_z_from(W%WDES%COMM_INTER,W1(NI)%CPROJ(1), &
                  SIZE(W1(NI)%CPROJ),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))
# 4364

          ENDIF
       ENDDO

       CALL M_waitall(nrequests,requests(1))

    ENDIF

    

  END SUBROUTINE W1_IGATHER_STRIP


!************************* SUBROUTINE WVREAL_PRECISE *******************
!
!> this subroutine forces the wavefunction to be real at the Gamma-point
!> it is required for the gamma point only mode
!> to avoid that small non real components develop
!> this version is exact and works through an FFT to real space
!> and then forcing the wavefunction to be real
!> the routine is required only if subspace rotations are performed
!> since the subspace rotation routine can not force wavefunctions
!> to become real
!>
!> if LORBITALREAL is .TRUE., the routine also tries to make the
!> orbitals real (in real space) by calculating the orbital
!> ~~~
!>   phi(r) = u_k(r) e^ikr
!> ~~~
!> taking the real part, storing the coefficients back
!> this is possible at special k-points such as points at the BZ boundary
!> ~~~
!> e.g. k=(+-0.5, +-0.5, +-0.5)
!>      k=(+-0.5, +-0.5,  0)
!>      k=(+-0.5,  0  ,  0)
!> ~~~
!>and permutations thereof
!>
!> the routine will not work properly at any other k-points
!> (and is in fact by-passed in this case)
!> ideally no subspace rotation should be called after calling this
!> routine (though orthogonalization is fine)
!
!***********************************************************************

  SUBROUTINE WVREAL_PRECISE(W)

!! USE moffload

    USE constant
    USE random_seeded, ONLY: RANE
    IMPLICIT NONE
    TYPE (wavespin) W
! local
    INTEGER :: NK, ISP, NB, ISPINOR, N
    TYPE (wavedes1)    WDES1          ! descriptor for (1._q,0._q) k-point
    TYPE (wavefun1)    W1             ! current wavefunction
    LOGICAL :: LPHASE
    COMPLEX(q) :: CPHASE(W%WDES%GRID%MPLWV)
    COMPLEX(q) :: AVERAGE_PHASE

! for spinors it is not (always) possible to make wave functions real-valued
! (for that reason we need to use the "_ncl" of VASP even for Gamma-only ...)
     IF (W%WDES%LNONCOLLINEAR) RETURN

# 4464

    IF (W%WDES%LORBITALREAL) THEN
!$ACC ENTER DATA CREATE(W1,WDES1) IF(OFFLOAD_ON)
    CALL SETWDES(W%WDES,WDES1,0)
    CALL NEWWAV(W1, WDES1, .TRUE.)
!
! force orbitals to be real for non-gamma point only method
!
    DO NK  =1,W%WDES%NKPTS

       IF (MOD(NK-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1) CYCLE

       IF (ABS(MOD(W%WDES%VKPT(1,NK)*2+100,1._q))< 1E-6 .AND. &
           ABS(MOD(W%WDES%VKPT(2,NK)*2+100,1._q))< 1E-6 .AND. &
           ABS(MOD(W%WDES%VKPT(3,NK)*2+100,1._q))< 1E-6) THEN
       CALL SETPHASE_WVREAL(W%WDES%VKPT(:,NK), W%WDES%GRID, CPHASE, LPHASE)
!$ACC ENTER DATA COPYIN(CPHASE) IF(OFFLOAD_ON .AND. LPHASE) ASYNC(ACC_ASYNC_Q)

       CALL SETWDES(W%WDES,WDES1,NK)
       DO ISP=1,W%WDES%ISPIN
          DO NB=1,W%WDES%NBANDS
             CALL W1_COPY(ELEMENT(W, WDES1, NB, ISP), W1)
             CALL FFTWAV_W1(W1)
! chose random phase in 1st and 4th quadrant, avoid proximity to +-i
             AVERAGE_PHASE=EXP(CITPI*(RANE()-0.5)/2.1)
             CALL M_bcast_z( WDES1%COMM_INB, AVERAGE_PHASE, 1)

             DO ISPINOR =0,WDES1%NRSPINORS-1
                IF (LPHASE) THEN
!$ACC PARALLEL LOOP PRESENT(W1,W1%CR,WDES1,WDES1%GRID,WDES1%GRID%MPLWV,CPHASE) &
!$ACC 
                   DO N=1,WDES1%GRID%RL%NP
                      W1%CR(N+ISPINOR*WDES1%GRID%MPLWV)=CONJG(CPHASE(N))*REAL(W1%CR(N+ISPINOR*WDES1%GRID%MPLWV)*CPHASE(N)*AVERAGE_PHASE,q) &
                      *(1.0_q/WDES1%GRID%NPLWV)
                   ENDDO
                ELSE
!$ACC PARALLEL LOOP PRESENT(W1,W1%CR,WDES1,WDES1%GRID,WDES1%GRID%MPLWV) 
                   DO N=1,WDES1%GRID%RL%NP
                      W1%CR(N+ISPINOR*WDES1%GRID%MPLWV)=REAL(W1%CR(N+ISPINOR*WDES1%GRID%MPLWV)*AVERAGE_PHASE,q) &
                      *(1.0_q/WDES1%GRID%NPLWV)
                   ENDDO
                ENDIF

                CALL FFTEXT(WDES1%NGVECTOR, WDES1%NINDPW(1), W1%CR(1+ISPINOR*WDES1%GRID%MPLWV),W%CPTWFP(1+ISPINOR*WDES1%NGVECTOR,NB,NK,ISP),WDES1%GRID,.FALSE.)
             ENDDO
          ENDDO
       ENDDO
!$ACC EXIT DATA DELETE(CPHASE) IF(OFFLOAD_ON .AND. LPHASE) ASYNC(ACC_ASYNC_Q)
       ENDIF
    ENDDO

    CALL DELWAV(W1, .TRUE.)
# 4519

    ENDIF


    RETURN
  END SUBROUTINE WVREAL_PRECISE


  SUBROUTINE SETPHASE_WVREAL(VKPT, GRID, CPHASE, LPHASE)
    USE constant

    IMPLICIT NONE

    REAL(q) :: VKPT(3)
    TYPE (grid_3d) GRID
    COMPLEX(q) :: CPHASE(GRID%MPLWV)
    LOGICAL LPHASE
! local
    REAL(q),PARAMETER :: TINY=1E-6_q
    REAL(q) F1, F2, F3
    INTEGER NC, N, IND
    COMPLEX(q) C, CD, CSUM
    CSUM=0
# 4545


    IF (ABS(VKPT(1))>TINY .OR. ABS(VKPT(2))>TINY .OR. ABS(VKPT(3))>TINY) THEN
       LPHASE=.TRUE.
       F1=TPI/GRID%NGX*VKPT(1)
       F2=TPI/GRID%NGY*VKPT(2)
       F3=TPI/GRID%NGZ*VKPT(3)

       IF (GRID%RL%NFAST==3) THEN
          CD=EXP(CMPLX(0,F3,q))
          IND=0
          DO NC=1,GRID%RL%NCOL
             C=EXP(CMPLX(0,F1*(GRID%RL%I2(NC)-1)+F2*(GRID%RL%I3(NC)-1),q))
             DO N=1,GRID%RL%NROW
                IND=IND+1
                CPHASE(IND)=C
                C=C*CD
             ENDDO
          ENDDO
       ELSE
          CD=EXP(CMPLX(0,F1,q))
          IND=0
          DO NC=1,GRID%RL%NCOL
             C=EXP(CMPLX(0,F2*(GRID%RL%I2(NC)-1)+F3*(GRID%RL%I3(NC)-1),q))
             DO N=1,GRID%RL%NROW
                IND=IND+1
                CPHASE(IND)=C
                CSUM=CSUM+C
                C=C*CD
             ENDDO
          ENDDO
       ENDIF
       LPHASE=.TRUE.
    ELSE
       LPHASE=.FALSE.
    ENDIF
  END SUBROUTINE SETPHASE_WVREAL

END MODULE wave_high


!***********************************************************************
!
!> assign a W1 structure an initial value from a wavefunction
!> array
!>
!> the plane wave coefficients and the wave function character (optional)
!> must be supplied
!> for performance reason no runtime checking on anything is performed
!
!***********************************************************************

  SUBROUTINE ARRAY_TO_W1( W1, C, CPROJ)
    USE prec
    USE wave
    IMPLICIT NONE
    TYPE (wavefun1)    W1
    COMPLEX(q):: C(*)
    COMPLEX(q), OPTIONAL :: CPROJ(*)

# 4607

    CALL ZCOPY( W1%WDES1%NRPLWV, C, 1, W1%CPTWFP(1), 1)

    IF (PRESENT(CPROJ)) THEN
       IF (W1%WDES1%LGAMMA) THEN
          CALL DCOPY( W1%WDES1%NPROD, CPROJ(1), 1, W1%CPROJ(1), 1)
       ELSE
          CALL ZCOPY( W1%WDES1%NPROD, CPROJ(1), 1, W1%CPROJ(1), 1)
       ENDIF
    ENDIF

  END SUBROUTINE ARRAY_TO_W1


!***********************************************************************
!
!> assign a W1 structure to an array (reverse of the previous operation)
!
!***********************************************************************

  SUBROUTINE W1_TO_ARRAY( W1, C,  CPROJ)
    USE prec
    USE wave
    IMPLICIT NONE
    TYPE (wavefun1)    W1
    COMPLEX(q) :: C(*)
    COMPLEX(q), OPTIONAL :: CPROJ(*)

# 4637

    CALL ZCOPY( W1%WDES1%NRPLWV, W1%CPTWFP(1), 1, C, 1)

    IF (PRESENT(CPROJ)) THEN
       IF (W1%WDES1%LGAMMA) THEN
          CALL DCOPY( W1%WDES1%NPROD,  W1%CPROJ(1), 1, CPROJ(1), 1)
       ELSE
          CALL ZCOPY( W1%WDES1%NPROD,  W1%CPROJ(1), 1, CPROJ(1), 1)
       ENDIF
    ENDIF

  END SUBROUTINE W1_TO_ARRAY


!************************* SUBROUTINE ECCP_NL **************************
!
!> this subroutine calculates the expectation value of < c|H|cp>
!> where c and cp are two wavefunctions;
!>
!> non local part only for (1._q,0._q) ion. I have put this in a
!> separate routine because optimization
!> is than easier
!
!***********************************************************************

  SUBROUTINE ECCP_NL(LMDIM,LMMAXC,CDIJ,CPROJ1,CPROJ2,CNL)
!$ACC ROUTINE VECTOR
    USE prec
    IMPLICIT NONE
    COMPLEX(q)      CNL
    INTEGER LMDIM, LMMAXC
    COMPLEX(q) CDIJ(LMDIM,LMDIM)
    COMPLEX(q) CPROJ1(LMMAXC),CPROJ2(LMMAXC)
! local
    INTEGER L, LP

!   

!DIR$ IVDEP
!OCL NOVREL
!$ACC LOOP VECTOR COLLAPSE(2) REDUCTION(+:CNL)
    DO L=1,LMMAXC
       DO LP=1,LMMAXC
          CNL=CNL+CDIJ(LP,L)*CPROJ1(LP)*CONJG(CPROJ2(L))
       ENDDO
    ENDDO

!   

  END SUBROUTINE ECCP_NL


!************************* SUBROUTINE OVERL ***************************
!
!>
!> calculate the result of the overlap-operator acting onto a set of
!> wave function characters; F77 low level routine
!> ~~~
!>  CRESUL^n_N, nlm = sum_n'l'm' D_N, n'l'm',nlm CPROF^n_N, n'l'm'
!>  CRESUL^n_N, i   = sum_j D_N, j, i CPROF^n_N, j
!> ~~~
!>
!> n is the band index, N the ion index, nlm = i the index for the
!> (1._q,0._q) centre partial waves
!>
!> @details @ref openmp :
!> the loop over the bands owned locally by a particular 1-rank
!> is distributed over all available OpenMP threads.
!
!**********************************************************************

  SUBROUTINE OVERL(WDES1, LOVERL, LMDIM, CQIJ, CPROF, CRESUL)
!! USE moffload
    USE wave

!$ACC ROUTINE(SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX) SEQ

    IMPLICIT NONE

    TYPE (wavedes1) WDES1
    LOGICAL LOVERL
    INTEGER LMDIM
    COMPLEX(q) CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
    COMPLEX(q) CRESUL(WDES1%NPROD,WDES1%NBANDS),CPROF(WDES1%NPROD,WDES1%NBANDS)
! local
    COMPLEX(q) CTMP
    INTEGER NB,NP,ISPINOR,ISPINOR_,NPRO,NPRO_,NT,NI,LMMAXC,L,LP

    

    IF (LOVERL) THEN
# 4732

 !$OMP PARALLEL DO DEFAULT(NONE) &
 !$OMP SHARED(WDES1,CRESUL,CQIJ,CPROF) &
 !$OMP PRIVATE(NB,NP,ISPINOR,ISPINOR_,NPRO,NPRO_,NI,NT,LMMAXC,L,LP,CTMP)
!$ACC PARALLEL LOOP PRESENT(WDES1,CQIJ,CPROF,CRESUL) COLLAPSE(4) GANG  &
!$ACC PRIVATE(NP,ISPINOR,ISPINOR_,NI,NT,LMMAXC,NPRO,NPRO_,L,LP,CTMP)
       bands: DO NB=1,WDES1%NBANDS
     DO NP=1,WDES1%NPRO
        CRESUL(NP,NB)=0
     ENDDO

          spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
          DO ISPINOR_=0,WDES1%NRSPINORS-1
             DO NI=1,WDES1%NIONS
                NT=WDES1%ITYP(NI)
                LMMAXC=WDES1%LMMAX(NT)
                IF (LMMAXC==0) CYCLE
                NPRO =WDES1%LMBASE(NI)+ISPINOR *(WDES1%NPRO/2)
                NPRO_=WDES1%LMBASE(NI)+ISPINOR_*(WDES1%NPRO/2)

!$ACC LOOP
!NEC$ select_vector
                DO L =1,LMMAXC

                   CTMP=0
!DIR$ IVDEP
!OCL NOVREC
!$ACC LOOP VECTOR REDUCTION(+:CTMP)
                   DO LP=1,LMMAXC
                      CTMP=CTMP+CQIJ(LP,L,NI,1+ISPINOR_+2*ISPINOR)*CPROF(LP+NPRO_,NB)
                   ENDDO

!$ACC ATOMIC UPDATE
                   CRESUL(L+NPRO,NB)=CRESUL(L+NPRO,NB)+CTMP
# 4768

# 4774

                ENDDO
             ENDDO
          ENDDO
          ENDDO spinor

       ENDDO bands
 !$OMP END PARALLEL DO
    ENDIF

    

  END SUBROUTINE OVERL


!************************* SUBROUTINE OVERL1 **************************
!
!
!> calculate the result of the overlap-operator acting onto (1._q,0._q)
!> wave function character;  F77 low level routine
!> ~~~
!>  CRESUL_N, nlm = sum_n'l'm' (D_N, n'l'm',nlm-e Q_N, n'l'm',nlm) CPROF_N, n'l'm'
!>  CRESUL_N, i   = sum_j (D_N, j, i-e Q_N, j, i) CPROF_N, j
!> ~~~
!>
!> N the ion index, nlm = i the index for the (1._q,0._q) centre partial waves
!>
!> @details @ref openmp :
!> the nested loops over \"types\" + \"ions-of-type\" are replaced
!> by a single loop over \"all ions\" that is distributed over
!> all available threads.
!
!**********************************************************************


!
! scheduled for removal to be replaced by OVERL1_
!
  SUBROUTINE OVERL1(WDES1, LMDIM, CDIJ, CQIJ, EVALUE, CPROF, CRESUL)

!! USE moffload_struct_def

    USE wave
    IMPLICIT NONE

    TYPE (wavedes1) WDES1
    INTEGER LMDIM
    COMPLEX(q) CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS), &
         CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
    REAL(q) :: EVALUE
    COMPLEX(q) CRESUL(WDES1%NPRO),CPROF(WDES1%NPRO)
! local
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NPRO2, NPRO2_, NT, NI, LMMAXC, L, LP
    COMPLEX(q) CTMP

    

!$ACC KERNELS PRESENT(CRESUL) 
    CRESUL=0
!$ACC END KERNELS

!$ACC ENTER DATA COPYIN(EVALUE) 

    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
    DO ISPINOR_=0,WDES1%NRSPINORS-1

       NPRO =ISPINOR *(WDES1%NPRO/2)
       NPRO_=ISPINOR_*(WDES1%NPRO/2)
# 4845

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
!$OMP SHARED(NPRO,NPRO_,WDES1,EVALUE,ISPINOR,ISPINOR_,CDIJ,CQIJ,CPROF,CRESUL) &
!$OMP PRIVATE(NI,NT,LMMAXC,NPRO2,NPRO2_,L,LP,CTMP)

!NEC$ novector
       DO NI=1,WDES1%NIONS
          NT=WDES1%ITYP(NI)
          LMMAXC=WDES1%LMMAX(NT)
          IF (LMMAXC==0) CYCLE
          NPRO2 =WDES1%LMBASE(NI)+NPRO
          NPRO2_=WDES1%LMBASE(NI)+NPRO_
          IF (EVALUE==0) THEN
!$ACC LOOP VECTOR PRIVATE(CTMP)
!NEC$ select_vector
             DO L =1,LMMAXC

                CTMP=0
!DIR$ IVDEP
!OCL NOVREC
                DO LP=1,LMMAXC
                   CTMP=CTMP + &
                        CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)*CPROF(LP+NPRO2_)
                ENDDO
                CRESUL(L+NPRO2)=CRESUL(L+NPRO2)+CTMP
# 4876

             ENDDO
          ELSE
!$ACC LOOP VECTOR PRIVATE(CTMP)
!NEC$ select_vector
             DO L =1,LMMAXC

                CTMP=0
!DIR$ IVDEP
!OCL NOVREC
                DO LP=1,LMMAXC
                   CTMP=CTMP+ &
                       (CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)- &
                        EVALUE*CQIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)) * CPROF(LP+NPRO2_)
                ENDDO
                CRESUL(L+NPRO2)=CRESUL(L+NPRO2)+CTMP
# 4899

             ENDDO
          ENDIF
       ENDDO

!$OMP END PARALLEL DO

    ENDDO
    ENDDO spinor

!$ACC EXIT DATA DELETE(EVALUE) 

    

  END SUBROUTINE OVERL1

# 4991


!
! identical to previous version but with complex EVALUE
!
!> @details @ref openmp :
!> the nested loops over \"types\" + \"ions-of-type\" are replaced
!> by a single loop over \"all ions\" that is distributed over
!> all available threads.
!
  SUBROUTINE OVERL1_C(WDES1, LMDIM, CDIJ, CQIJ, EVALUE, CPROF,CRESUL)

!! USE moffload_struct_def

    USE wave
    IMPLICIT NONE

    TYPE (wavedes1) WDES1
    INTEGER LMDIM
    COMPLEX(q) CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS), &
         CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
    COMPLEX(q) EVALUE
    COMPLEX(q) CRESUL(WDES1%NPRO),CPROF(WDES1%NPRO)
! local
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NPRO2, NPRO2_, NT, NI, LMMAXC, L, LP
    COMPLEX(q) CTMP

    

!$ACC KERNELS PRESENT(CRESUL) 
    CRESUL=0
!$ACC END KERNELS

!$ACC ENTER DATA COPYIN(EVALUE) 

    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
    DO ISPINOR_=0,WDES1%NRSPINORS-1

       NPRO =ISPINOR *(WDES1%NPRO/2)
       NPRO_=ISPINOR_*(WDES1%NPRO/2)
# 5034

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
!$OMP SHARED(NPRO,NPRO_,WDES1,EVALUE,ISPINOR,ISPINOR_,CDIJ,CQIJ,CPROF,CRESUL) &
!$OMP PRIVATE(NI,NT,LMMAXC,NPRO2,NPRO2_,L,LP,CTMP)

       DO NI=1,WDES1%NIONS
          NT=WDES1%ITYP(NI)
          LMMAXC=WDES1%LMMAX(NT)
          IF (LMMAXC==0) CYCLE
          NPRO2 =WDES1%LMBASE(NI)+NPRO
          NPRO2_=WDES1%LMBASE(NI)+NPRO_
          DO L =1,LMMAXC
             CTMP=0
!DIR$ IVDEP
!OCL NOVREC
!$ACC LOOP VECTOR REDUCTION(+:CTMP)
             DO LP=1,LMMAXC
                CTMP=CTMP+(CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)- &
                     EVALUE*CQIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)) * CPROF(LP+NPRO2_)
             ENDDO
             CRESUL(L+NPRO2)=CRESUL(L+NPRO2)+CTMP
          ENDDO
       ENDDO

!$OMP END PARALLEL DO

    ENDDO
    ENDDO spinor

!$ACC EXIT DATA DELETE(EVALUE) 

    

  END SUBROUTINE OVERL1_C

!
!> indentical to #OVERL1 but with CDIJ and CQIJ always complex
!
  SUBROUTINE OVERL1_CCDIJ(WDES1, LMDIM, CDIJ, CQIJ, EVALUE, CPROF,CRESUL)
    USE wave
    IMPLICIT NONE

    TYPE (wavedes1) WDES1
    INTEGER LMDIM
    COMPLEX(q) CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS), &
         CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
    REAL(q) :: EVALUE
    COMPLEX(q) CRESUL(WDES1%NPRO),CPROF(WDES1%NPRO)
! local
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NT, NIS, NI, LMMAXC, L, LP

    CRESUL=0
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
    DO ISPINOR_=0,WDES1%NRSPINORS-1

       NPRO =ISPINOR *(WDES1%NPRO/2)
       NPRO_=ISPINOR_*(WDES1%NPRO/2)

       NIS =1
       DO NT=1,WDES1%NTYP
          LMMAXC=WDES1%LMMAX(NT)
          IF (LMMAXC/=0) THEN
             DO NI=NIS,WDES1%NITYP(NT)+NIS-1
                IF (EVALUE==0) THEN
                   DO L =1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
                   DO LP=1,LMMAXC
                      CRESUL(L+NPRO)=CRESUL(L+NPRO)+ &
                           CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)*CPROF(LP+NPRO_)
                   ENDDO
                   ENDDO
                ELSE
                   DO L =1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
                   DO LP=1,LMMAXC
                      CRESUL(L+NPRO)=CRESUL(L+NPRO)+ &
                           (CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)- &
                           EVALUE*CQIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)) * CPROF(LP+NPRO_)
                   ENDDO
                   ENDDO
                ENDIF
                NPRO = LMMAXC+NPRO
                NPRO_= LMMAXC+NPRO_
             ENDDO
          ENDIF
          NIS = NIS+WDES1%NITYP(NT)
       ENDDO
    ENDDO
    ENDDO spinor

  END SUBROUTINE OVERL1_CCDIJ
