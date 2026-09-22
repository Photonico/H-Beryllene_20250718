# 1 "dfast.F"
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


# 2 "dfast.F" 2 
!***********************************************************************
!
!> Collection of a few high level BLAS 3 routines
!>
!> Calculate inproducts between sets of wave functions
!> and to perform transformation (sub space rotations)
!
!***********************************************************************
MODULE dfast
  USE prec
  USE mpimy
!
!> @details The overlap or Hamiltonmatrix between states is calculated blockwise
!> if the upper of lower triangle is required
!> NSTRIPD is the maximum value ever used, whereas NSTRIP_STANDARD
!> is the typical value used
  INTEGER, PARAMETER :: NSTRIPD=16
!> Typical value used, see #NSTRIPD
  INTEGER, SAVE :: NSTRIP_STANDARD
  INTEGER, SAVE :: NSTRIP_STANDARD_GLOBAL
!
!> When the wavefunctions are rotated (unitary transformations)
!> this is 1._q in blocks to save storage for the transformed wavefunctions
  INTEGER :: NBLK=256

  INTERFACE
!************************ SUBROUTINE ORTH1 *****************************
!> Calculates a stripe of columns of the overlap between two set of vectors
!>
!> Actual documentation of this subroutine is in ::ORTH1()
!***********************************************************************
     SUBROUTINE ORTH1(CSEL,CPTWFP,CFW,CPROJ,CPROW,NBANDS, &
          &  NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL)
       USE prec
       IMPLICIT COMPLEX(q) (C)
       IMPLICIT REAL(q) (A-B,D-H,O-Z)
       
       COMPLEX(q)      CPROJ
       COMPLEX(q)      CPROW
       COMPLEX(q)      COVL
       CHARACTER (LEN=*) CSEL
     END SUBROUTINE ORTH1
  END INTERFACE
  
  INTERFACE
!************************ SUBROUTINE ORTH2 *****************************
!> Calculates a stripe of columns of the overlap between two set of vectors
!>
!> Actual documentation of this subroutine is in ::ORTH2()
!***********************************************************************
     SUBROUTINE ORTH2(CPTWFP,CFW,CPROJ,CPROW,NBANDS, &
          &  NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL)
       USE prec
       IMPLICIT COMPLEX(q) (C)
       IMPLICIT REAL(q) (A-B,D-H,O-Z)

       COMPLEX(q)      CPROJ
       COMPLEX(q)      CPROW
       COMPLEX(q)      COVL
     END SUBROUTINE ORTH2
  END INTERFACE

  INTERFACE
!************************ SUBROUTINE ORTH3 *****************************
!> Similar to ::ORTH2 but allows for a different leading dimension NDIM
!> for the matrix to be set up (COVL)
!>
!> Actual documentation of this subroutine is in ::ORTH3()
!***********************************************************************
     SUBROUTINE ORTH3(CPTWFP,CFW,CPROJ,CPROW, &
          & NBANDS,NSTRIP,NDIM,NPOS, &
          & NPL,NPRO,NPLDIM,NPROD,COVL)
       USE prec
       IMPLICIT COMPLEX(q) (C)
       IMPLICIT REAL(q) (A-B,D-H,O-Z)

       COMPLEX(q)      CPROJ
       COMPLEX(q)      CPROW
       COMPLEX(q)      COVL
     END SUBROUTINE ORTH3
  END INTERFACE


  TYPE parallel_gemm
     INTEGER :: N          !< dimension of matrix
     TYPE (communic) :: COMM
     INTEGER :: NPROC      !<  number of cores
     INTEGER, ALLOCATABLE :: NCOL(:), OFFSET(:), NCTOT(:), OFFDATA(:)
  END TYPE

  CONTAINS
!***********************************************************************
!
!> Set default for the blocking
!
!***********************************************************************

    SUBROUTINE SET_NBLK_NSTRIP(WDES)
      USE wave
      TYPE (wavedes)  WDES
      
      INTEGER NCPU


      NCPU   =WDES%COMM_INTER%NCPU ! number of procs involved in band dis.
# 109

! use typically (1._q,0._q) quarter of plane wave coefficients but never more than 256
# 113

      IF (NBLK==-1) NBLK=MIN(256,MAX(32,(WDES%NRPLWV/256)*64))

      CALL M_max_i(WDES%COMM, NBLK, 1)
! set NSTRIP between [1 and NSTRIPD]
      NSTRIP_STANDARD=MAX(MIN(NSTRIPD,32/NCPU,WDES%NBANDS),1)
      NSTRIP_STANDARD_GLOBAL=NSTRIP_STANDARD*NCPU

    END SUBROUTINE SET_NBLK_NSTRIP


!************************ SUBROUTINE LINCOM ****************************
!
!> Build linear combinations of wavefunctions according to matrix CTRANS
!>
!> This subroutine performes implicitly a MATRIX x MATRIX multiplication,
!> it is needed for the unitary transformation of the wavefunctions or
!> for orthogonalisation routines and uses a blocked algorithm
!> to save storage.
!> ~~~
!> COUT_n,k =   sum_kp CIN_n,kp CTRANS kp,k
!> where n  =  1 ... NPL    (leading dimension of array = NPLDIM)
!>       k  =  1 ... NOUT
!>       kp =  1 ... NIN    (NIN must be greater or equal NOUT)
!> ~~~
!> @warning On exit the input array will be overwritten by output data!
!>
!> `MODE` determines the mode for the transformation
!> - `U`  `CTRANS`   upper triangle set
!> - `L`  `CTRANS`   lower triangle set
!> - `F`  `CTRANS`   all components set
!> - `C`  `CTRANS`   contains the conjugated transformation matrix
!> - `T`  `CTRANS`   contains the transposed transformation matrix
!>
!> - `A`  used only for Davidson
!> - `B`  used only for Davidson
!>
!> #LINCOM_BLAS has a BLAS like calling sequence and does not
!> allow for the Davidson like functionality
!
!***********************************************************************


    SUBROUTINE LINCOM(MODE,CF,CPROF,CTRANS,NIN,NOUT,NPL, &
       &           NPRO,NPLDIM,NPROD,LDTRAN,CFA,CPROFA)

!! USE moffload_struct_def

      USE prec
      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      CHARACTER (1) MODE
      DIMENSION   CF(NPLDIM,NIN),CFA(NPLDIM,NIN)
      COMPLEX(q)        CPROF(NPROD,NIN),CPROFA(NPROD,NIN)
      COMPLEX(q)        CTRANS(LDTRAN,NIN)
# 171


! work array
      COMPLEX(q),ALLOCATABLE ::   CBLOCK(:,:)

      

# 181

      ALLOCATE(CBLOCK(NBLK,LDTRAN))
! CBLOCK=(0.0_q, 0.0_q) ! No test failures reported.

      IF (NPL/=0) THEN
      CALL LINBAS(MODE,CF,CBLOCK,CTRANS,NIN,NOUT, NPL, &
     &             NPLDIM,LDTRAN, NBLK,CFA)
      ENDIF
      IF (NPRO/=0) THEN
      CALL LINBAS(MODE,CPROF,CBLOCK,CTRANS,NIN,NOUT, NPRO, &
     &             NPROD,LDTRAN, NBLK,CPROFA)
      ENDIF
      DEALLOCATE(CBLOCK)

# 197

      

      RETURN
    END SUBROUTINE LINCOM


    SUBROUTINE LINCOM_BLAS(MODE, NPL, NPRO, NIN, NOUT, CF, NPLDIM, CPROF, NPROD, CTRANS, LDTRAN)
      USE prec

      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      CHARACTER (1) MODE
      DIMENSION   CF(NPLDIM,NIN)
      COMPLEX(q)        CPROF(NPROD,NIN)
      COMPLEX(q)        CTRANS(LDTRAN,NIN)

! work array
      COMPLEX(q),ALLOCATABLE ::   CBLOCK(:,:)
      ALLOCATE(CBLOCK(NBLK,LDTRAN))
! CBLOCK=(0.0_q, 0.0_q) ! No test failures reported.

      CALL LINBAS_BLAS(MODE,  NPL, NIN, NOUT, CF,  NPLDIM, CTRANS, LDTRAN, CBLOCK,  NBLK)
      IF (NPRO/=0) THEN
      CALL LINBAS_BLAS(MODE, NPRO,   NIN, NOUT, CPROF, NPROD,  CTRANS, LDTRAN, CBLOCK,  NBLK)
      ENDIF
      DEALLOCATE(CBLOCK)

      RETURN
    END SUBROUTINE LINCOM_BLAS

!************************ SUBROUTINE SETUP_PARALLEL_GEMM ***************
!
!> Sets up a simple communicator structure to perform parallel DGEMM calls
!>
!> The source matrices need to be known on all nodes.
!> The final "destination" matrix is assembled using an
!> allgatherv.
!
!***********************************************************************

    SUBROUTINE SETUP_PARALLEL_GEMM(COMM, NB_TOT, PGEMM_HANDLE)
      USE mpimy
      IMPLICIT NONE
      TYPE(communic) :: COMM          !< used communicator
      INTEGER        :: NB_TOT        !< matrix dimension
      TYPE (parallel_gemm), POINTER ::  PGEMM_HANDLE
! local
      INTEGER        :: I


      ALLOCATE(PGEMM_HANDLE)

      PGEMM_HANDLE%NPROC =COMM%NCPU
      PGEMM_HANDLE%N     =NB_TOT
      PGEMM_HANDLE%COMM  =COMM

      ALLOCATE( PGEMM_HANDLE%NCOL(PGEMM_HANDLE%NPROC),  PGEMM_HANDLE%OFFSET(PGEMM_HANDLE%NPROC), & 
                PGEMM_HANDLE%NCTOT(PGEMM_HANDLE%NPROC), PGEMM_HANDLE%OFFDATA(PGEMM_HANDLE%NPROC) )

      CALL DISTRIBUTE_COLUMNS_GEMM(NB_TOT, PGEMM_HANDLE%NPROC, PGEMM_HANDLE%NCOL)

      PGEMM_HANDLE%OFFSET   = 0
      PGEMM_HANDLE%OFFDATA  = 0
      PGEMM_HANDLE%NCTOT    = 0

      PGEMM_HANDLE%NCTOT(1) = PGEMM_HANDLE%NCOL(1) * NB_TOT

      DO I=2,PGEMM_HANDLE%NPROC
         IF (PGEMM_HANDLE%NCOL(I)>0) THEN
            PGEMM_HANDLE%NCTOT(I)   = NB_TOT * PGEMM_HANDLE%NCOL(I)
            PGEMM_HANDLE%OFFSET(I)  = PGEMM_HANDLE%OFFSET(I-1) + PGEMM_HANDLE%NCOL(I-1)
            PGEMM_HANDLE%OFFDATA(I) = NB_TOT * PGEMM_HANDLE%OFFSET(I)
         ENDIF
      ENDDO


      CONTAINS

!
!> Calculates number of columns per core
!
      SUBROUTINE  DISTRIBUTE_COLUMNS_GEMM( N, NPROC, NCOL )

        INTEGER :: N_PER_PROC, REMAINDER, ME1, NPROC, N, NCOL(NPROC)
        INTEGER :: NCOL_MIN=16
        INTEGER :: NPROC_MAX

        NCOL = 0

! only NPROC_MAX cores take part in the communication
! otherwise communication gets too expensive
        NPROC_MAX=MIN(MAX(1,N/NCOL_MIN),NPROC)
        
        N_PER_PROC = N/NPROC_MAX
        REMAINDER   = N - NPROC_MAX * N_PER_PROC
        DO ME1 = 1,NPROC_MAX
           IF( ME1 <= REMAINDER )THEN
              NCOL(ME1) = N_PER_PROC + 1
           ELSE
              NCOL(ME1) = N_PER_PROC
           ENDIF
        ENDDO
        
        RETURN
      END SUBROUTINE DISTRIBUTE_COLUMNS_GEMM
      
     
    END SUBROUTINE SETUP_PARALLEL_GEMM


    SUBROUTINE RELEASE_PARALLEL_GEMM(PGEMM_HANDLE)
      IMPLICIT NONE

      TYPE (parallel_gemm), POINTER ::  PGEMM_HANDLE

      IF (ASSOCIATED(PGEMM_HANDLE)) THEN
         DEALLOCATE( PGEMM_HANDLE%NCOL,  PGEMM_HANDLE%OFFSET, & 
                     PGEMM_HANDLE%NCTOT, PGEMM_HANDLE%OFFDATA )
         DEALLOCATE(PGEMM_HANDLE)
      ENDIF

    END SUBROUTINE RELEASE_PARALLEL_GEMM

!************************ SUBROUTINE PARALLEL_GEMM *********************
!
!> Same interface as ZGEMM except for the handle which needs to be passed
!> along as well
!>
!> @note This version works in principle only for square matrices
!> and is not intended for general purpose use.
!> I think also the leading dimensions must be correct and identical.
!
!***********************************************************************

    SUBROUTINE PARALLEL_GGEMM( PGEMM_HANDLE, TRANSA , TRANSB, N1, N2, N3, ALPHA, A, LDA, &
         B, LDB, BETA, C, LDC)

!! USE moffload

      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      CHARACTER(LEN=1), INTENT(IN) :: TRANSA
      CHARACTER(LEN=1), INTENT(IN) :: TRANSB
      INTEGER N1, N2, N3
      INTEGER LDA, LDB, LDC
      COMPLEX(q)             :: ALPHA, BETA
      COMPLEX(q)             :: A(LDA,*),B(LDB,*),C(LDC,*)
      TYPE (parallel_gemm), POINTER ::  PGEMM_HANDLE
! local
      INTEGER NODE_ME

      COMPLEX(q) :: CTMP(LDC,LDC)

      IF (.NOT. ASSOCIATED(PGEMM_HANDLE)) &
         CALL vtutor%bug("PARALLEL_GGEMM: PGEMM_HANDLE not set up", "dfast.F", 355)

      IF (N1 /= N2 .OR. N2 /= N3 .OR. LDC /= N1) &
         CALL vtutor%bug("PARALLEL_GGEMM: not a square matrix " // str(N1) // " " &
            // str(N2) // " " // str(N3) // " " // str(LDC), "dfast.F", 359)

      IF (N1 /= PGEMM_HANDLE%N) &
         CALL vtutor%bug("PARALLEL_GGEMM: not correctly set up " // &
            str(N1 /= PGEMM_HANDLE%N), "dfast.F", 363)

      IF (TRANSB=='N' .OR. TRANSB=='n') THEN

         NODE_ME = PGEMM_HANDLE%COMM%NODE_ME

         IF (PGEMM_HANDLE%NCOL(NODE_ME)>0) &
            CALL ZGEMM( TRANSA, TRANSB, N1,  PGEMM_HANDLE%NCOL(NODE_ME), N3, ALPHA, A, LDA, &
                 B(1,1+PGEMM_HANDLE%OFFSET(NODE_ME)), LDB, BETA, C(1,1+PGEMM_HANDLE%OFFSET(NODE_ME)), LDC )

!$ACC ENTER DATA CREATE(CTMP) 
!$ACC KERNELS PRESENT(CTMP,C) 
         CTMP(1:LDC,1:LDC)=C(1:LDC,1:LDC)
!$ACC END KERNELS
# 380

         CALL M_allgatherv_z (PGEMM_HANDLE%COMM, CTMP(1,1+PGEMM_HANDLE%OFFSET(NODE_ME)),  &
                 C, PGEMM_HANDLE%NCTOT, PGEMM_HANDLE%OFFDATA)

!$ACC EXIT DATA DELETE(CTMP) 
      ELSE
         CALL vtutor%bug("PARALLEL_GGEMM: the second matrix needs to be stored 'N'", "dfast.F", 386)
      ENDIF
# 390


    END SUBROUTINE PARALLEL_GGEMM
END MODULE dfast


!************************ SUBROUTINE LINBAS ****************************
!
!> Build linear combinations of set of vectors according to matrix CTRANS
!>
!> This subroutine performes implicitly a MATRIX x MATRIX multiplication,
!> it is needed for the unitary transformation of the wavefunctions or
!> for orthogonalisation routines and uses a blocked algorithm
!> to save storage.
!>
!> LINBAS is only called from LINCOM.
!>
!> LINBAS_BLAS allows only in place transformation of a set of
!> of orbitals without adding another set .
!> Furthermore the calling sequence is more BLAS like
!
!***********************************************************************

    SUBROUTINE LINBAS(MODE,CF,CBLOCK,CTRANS,NIN,NOUT,NPL, &
     &           NPLDIM,LDTRAN,NBLK,CFA)

!! USE moffload_struct_def
!! USE moffload


      USE prec
      USE tutor, ONLY: vtutor

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

! Might be that ZGEMM performs much much faster than ZTRMM ... ??
# 429

      PARAMETER(IUSETR=1)


      CHARACTER (1) MODE
      LOGICAL     LTRI,LADD,LBOTH,LTRANS
      COMPLEX(q)     CF(NPLDIM,NIN),CFA(NPLDIM,NIN)
      COMPLEX(q)     CBLOCK(NBLK,LDTRAN)
      COMPLEX(q)     CTRANS(LDTRAN,NIN)

      IF (NOUT>NIN) THEN
         CALL vtutor%bug("LINBAS: wrong arguments, NOUT>NIN", "dfast.F", 440)
      ENDIF

      LTRI =(MODE=='U').OR.(MODE=='u').OR. &
     &      (MODE=='L').OR.(MODE=='l')
      LADD =(MODE=='A').OR.(MODE=='a')
      LBOTH=(MODE=='B').OR.(MODE=='b')
      LTRANS=(MODE=='T').OR.(MODE=='t').OR.(MODE=='C').OR.(MODE=='c')
      IF (LTRI.AND.(IUSETR==0)) THEN
         IF ((MODE=='L').OR.(MODE=='l')) THEN
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CTRANS) 
            DO N2=1,NIN
!DIR$ IVDEP
!OCL NOVREC
               DO N1=1,N2-1
                  CTRANS(N1,N2)= (0._q,0._q)
               ENDDO
            ENDDO
         ELSE
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CTRANS) 
            DO N2=1,NIN
!DIR$ IVDEP
!OCL NOVREC
               DO N1=N2+1,NIN
                  CTRANS(N1,N2)= (0._q,0._q)
               ENDDO
            ENDDO
         ENDIF
      ENDIF

! Try to get best load balance, maximum block size < NBLK ...
      NBLOCK=NBLK

      DO IBLOCK=0,NPL-1,NBLOCK
         ILENPL=MIN(NBLOCK,NPL-IBLOCK)
         IADDPL=MIN(IBLOCK,NPL-1)
         ILENPL=MAX(ILENPL,0)

         IF (LTRI.AND.(IUSETR/=0)) THEN
! 'Triangular update':
# 482

            CALL ZTRMM &

     &                ('R',MODE,'N','N', ILENPL,NOUT,(1._q,0._q), &
     &                 CTRANS,LDTRAN,CF(1+IADDPL,1), NPLDIM)
         ELSE
! 'Full update':
!$ACC ENTER DATA CREATE(CBLOCK) 
            IF (LBOTH.OR.(.NOT.LADD)) THEN
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CBLOCK,CF) 
               DO N1=1,NIN
                  DO M=1,ILENPL
                     CBLOCK(M,N1)=CF(M+IADDPL,N1)
                  ENDDO
               ENDDO
               IF (LTRANS) THEN
               CALL ZGEMM('N',MODE, ILENPL, NOUT, NIN, (1._q,0._q), &
     &               CBLOCK(1,1), NBLK, CTRANS(1,1), &
     &               LDTRAN, (0._q,0._q), CF(IADDPL+1,1), NPLDIM)
               ELSE
               CALL ZGEMM('N', 'N', ILENPL, NOUT, NIN, (1._q,0._q), &
     &               CBLOCK(1,1), NBLK, CTRANS(1,1), &
     &               LDTRAN, (0._q,0._q), CF(IADDPL+1,1), NPLDIM)
               ENDIF
            ENDIF
            IF (LBOTH.OR.LADD) THEN
               IADDT=0
               IF (LBOTH) IADDT=NIN
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CBLOCK,CF) 
               DO N1=1,NIN
                  DO M=1,ILENPL
                     CBLOCK(M,N1)=CFA(M+IADDPL,N1)
                  ENDDO
               ENDDO
               CALL ZGEMM('N', 'N', ILENPL, NOUT, NIN, (1._q,0._q), &
     &               CBLOCK(1,1), NBLK, CTRANS(1+IADDT,1), &
     &               LDTRAN, (1._q,0._q), CF(IADDPL+1,1), NPLDIM)
            ENDIF
!$ACC EXIT DATA DELETE(CBLOCK) 
         ENDIF

      ENDDO

      RETURN
    END SUBROUTINE


    SUBROUTINE LINBAS_BLAS(MODE, NPL, NIN, NOUT, CF,NPLDIM, CTRANS, LDTRAN, CBLOCK, NBLK)
      USE prec
      USE tutor, ONLY: vtutor

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

! Might be that ZGEMM performs much much faster than ZTRMM ... ??
# 539

      PARAMETER(IUSETR=1)


      CHARACTER (1) MODE
      LOGICAL     LTRI,LADD,LBOTH,LTRANS
      COMPLEX(q)     CF(NPLDIM,NIN)
      COMPLEX(q)     CBLOCK(NBLK,LDTRAN)
      COMPLEX(q)     CTRANS(LDTRAN,NIN)

      IF (NOUT>NIN) THEN
         CALL vtutor%bug("internal error in routine LINBAS: wrong arguments, NOUT>NIN", "dfast.F", 550)
      ENDIF

      LTRI =(MODE=='U').OR.(MODE=='u').OR. &
     &      (MODE=='L').OR.(MODE=='l')
      LADD =(MODE=='A').OR.(MODE=='a')
      LTRANS=(MODE=='T').OR.(MODE=='t').OR.(MODE=='C').OR.(MODE=='c')
      IF (LTRI.AND.(IUSETR==0)) THEN
         DO 4 N2=1,NIN
            IF ((MODE=='L').OR.(MODE=='l')) THEN
!DIR$ IVDEP
!OCL NOVREC
               DO N1=1,N2-1
               CTRANS(N1,N2)= (0._q,0._q)
               ENDDO
            ELSE
!DIR$ IVDEP
!OCL NOVREC
               DO N1=N2+1,NIN
               CTRANS(N1,N2)= (0._q,0._q)
               ENDDO
            ENDIF
    4    ENDDO
      ENDIF

! Try to get best load balance, maximum block size < NBLK ...
      NBLOCK=NBLK

      DO 70 IBLOCK=0,NPL-1,NBLOCK
         ILENPL=MIN(NBLOCK,NPL-IBLOCK)
         IADDPL=MIN(IBLOCK,NPL-1)
         ILENPL=MAX(ILENPL,0)

         IF (LTRI.AND.(IUSETR/=0)) THEN
! 'Triangular update':
# 587

            CALL ZTRMM &

     &                ('R',MODE,'N','N', ILENPL,NOUT,(1._q,0._q), &
     &                 CTRANS,LDTRAN,CF(1+IADDPL,1), NPLDIM)
         ELSE
! 'Full update':
            IF (LBOTH.OR.(.NOT.LADD)) THEN
               DO 30 N1=1,NIN
                  DO 10 M=1,ILENPL
                     CBLOCK(M,N1)=CF(M+IADDPL,N1)
   10             CONTINUE
   30          CONTINUE
               IF (LTRANS) THEN
               CALL ZGEMM('N',MODE, ILENPL, NOUT, NIN, (1._q,0._q), &
     &               CBLOCK(1,1), NBLK, CTRANS(1,1), &
     &               LDTRAN, (0._q,0._q), CF(IADDPL+1,1), NPLDIM)
               ELSE
               CALL ZGEMM('N', 'N', ILENPL, NOUT, NIN, (1._q,0._q), &
     &               CBLOCK(1,1), NBLK, CTRANS(1,1), &
     &               LDTRAN, (0._q,0._q), CF(IADDPL+1,1), NPLDIM)
               ENDIF
            ENDIF
         ENDIF

   70 CONTINUE

      RETURN
    END SUBROUTINE

!************************ SUBROUTINE ORTH1 *****************************
!> Calculates a stripe of columns of the overlap between two set of vectors
!>
!> ~~~
!>   O(I,J) =  <CPTWFP(I) |  CFW(J) > +  <CPROJ(I) | CPROW(J) >
!>       J=NPOS ,..., NPOS+NSTRIP-1
!>       I=1 ,..., NBANDS
!> ~~~
!> Usually CWF might be either equal CPTWFP, and `CPROW = Q | CPROJ>`
!> or it might hold
!> ~~~
!> H_kinetic + H_local CPTWFP, and CPROW = D | CPROJ>
!> ~~~
!> ORTH1 determines only the lower part of the resultant matrix
!>  and is suitable if the resulting matrix is Hermitian.
!>
!> ::ORTH2 calculates the entire matrix:
!> ~~~
!>                ....    CFW(NPOS) CFW(NPOS+1) ... CFW(NPOS+NSTRIP-1) ....
!>
!> C(1)*             X    O2         O2       O2       O2              X
!> ...               X    O2         O2       O2       O2              X
!> C(NPOS)*          X    C          C        C        C               X
!> ...               X    C          C        C        C               X
!> C(NPOS+NSTRIP-1)* X    C          C        C        C               X
!> ...               X    C          C        C        C               X
!> ~~~
!>
!> In a single call the matrix elements marked with C are calculated
!> by ORTH1, ::ORTH2 additionally calculates those marked with O2.
!> Those marked with X are not updated in a single call.
!>
!> Explicit interface to this subroutine is in dfast::ORTH1 (beginning of
!> dfast module).
!***********************************************************************

    SUBROUTINE ORTH1(CSEL,CPTWFP,CFW,CPROJ,CPROW,NBANDS, &
     &  NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL)
      USE prec
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      DIMENSION CPTWFP(NPLDIM,NBANDS)
      DIMENSION CFW(NPLDIM,NSTRIP)
      COMPLEX(q)      CPROJ(NPROD,NBANDS)
      COMPLEX(q)      CPROW(NPROD,NSTRIP)
      COMPLEX(q)      COVL(NBANDS,NBANDS)
      CHARACTER (LEN=*) CSEL

      

      IF (NSTRIP+NPOS-1 > NBANDS) THEN
        CALL vtutor%bug("internal error in ORTH1: dim= " // str(NSTRIP+NPOS) // " " // str(NBANDS), "dfast.F", 671)
      ENDIF
!
! update of lower triangular part
!
    IF (CSEL(1:1) == 'L' .OR. CSEL(1:1) == 'l') THEN
      IF (NPL/=0) THEN
! for historic reasons (1._q,0._q) can "unroll" the loop over-plane
! wave coefficients (not used at present)
      NBLOCK= NPL

      DO NPOSPL=1, NPL-NBLOCK,NBLOCK
      CALL ZGEMM('C','N',NBANDS-NPOS+1,NSTRIP,NBLOCK,(1._q,0._q), &
              CPTWFP(NPOSPL,NPOS), NPLDIM,CFW(NPOSPL,1), &
               NPLDIM,(1._q,0._q),COVL(NPOS,NPOS),NBANDS)
      ENDDO
      CALL ZGEMM('C','N',NBANDS-NPOS+1,NSTRIP, NPL-NPOSPL+1,(1._q,0._q), &
              CPTWFP(NPOSPL,NPOS), NPLDIM,CFW(NPOSPL,1), &
               NPLDIM,(1._q,0._q),COVL(NPOS,NPOS),NBANDS)
      ENDIF

      IF (NPRO/=0) THEN
! for historic reasons (1._q,0._q) can "unroll" the loop non-local
! wave coefficients (not used at present)
      NBLOCK=NPRO
      DO NPOSPR=1,NPRO-NBLOCK,NBLOCK
      CALL ZGEMM('C','N',NBANDS-NPOS+1,NSTRIP,NBLOCK,(1._q,0._q), &
              CPROJ(NPOSPR,NPOS), NPROD,CPROW(NPOSPR,1), &
              NPROD,(1._q,0._q),COVL(NPOS,NPOS),NBANDS)
      ENDDO
      CALL ZGEMM('C','N',NBANDS-NPOS+1,NSTRIP,NPRO-NPOSPR+1,(1._q,0._q), &
              CPROJ(NPOSPR,NPOS), NPROD,CPROW(NPOSPR,1), &
              NPROD,(1._q,0._q),COVL(NPOS,NPOS),NBANDS)
      ENDIF
!
! update of upper triangular part
!
    ELSE IF (CSEL(1:1) == 'U' .OR. CSEL(1:1) == 'u') THEN
      IF (NPL/=0) THEN
      NBLOCK= NPL

      DO NPOSPL=1, NPL-NBLOCK,NBLOCK
      CALL ZGEMM('C','N',NSTRIP+NPOS-1,NSTRIP,NBLOCK,(1._q,0._q), &
              CPTWFP(NPOSPL,1), NPLDIM,CFW(NPOSPL,1), &
               NPLDIM,(1._q,0._q),COVL(1,NPOS),NBANDS)
      ENDDO
      CALL ZGEMM('C','N',NSTRIP+NPOS-1,NSTRIP, NPL-NPOSPL+1,(1._q,0._q), &
              CPTWFP(NPOSPL,1), NPLDIM,CFW(NPOSPL,1), &
               NPLDIM,(1._q,0._q),COVL(1,NPOS),NBANDS)
      ENDIF

      IF (NPRO/=0) THEN

      NBLOCK=NPRO
      DO NPOSPR=1,NPRO-NBLOCK,NBLOCK
      CALL ZGEMM('C','N',NSTRIP+NPOS-1,NSTRIP,NBLOCK,(1._q,0._q), &
              CPROJ(NPOSPR,1), NPROD,CPROW(NPOSPR,1), &
              NPROD,(1._q,0._q),COVL(1,NPOS),NBANDS)
      ENDDO
      CALL ZGEMM('C','N',NSTRIP+NPOS-1,NSTRIP,NPRO-NPOSPR+1,(1._q,0._q), &
              CPROJ(NPOSPR,1), NPROD,CPROW(NPOSPR,1), &
              NPROD,(1._q,0._q),COVL(1,NPOS),NBANDS)
      ENDIF

    ELSE
      WRITE(*,*)'internal error in ORTH1: CSEL=',CSEL
    ENDIF

      

      RETURN
    END SUBROUTINE ORTH1



!************************ SUBROUTINE ORTH2 *****************************
!> Calculates a stripe of columns of the overlap between two set of vectors
!>
!> ~~~
!>   O(I,J) =  <CPTWFP(I) |  CFa(J) > +  <CPROJ(I) | CPROW(J) >
!>       J=NPOS ,..., NPOS+NSTRIP-1
!>       I=1 ,..., NBANDS
!> ~~~
!> For general case (i.e. resulting matrix must not be Hermitian)
!> see comments in ::ORTH1.
!>
!> Explicit interface to this subroutine is in dfast::ORTH2 (beginning of
!> dfast module).
!***********************************************************************

    SUBROUTINE ORTH2(CPTWFP,CFW,CPROJ,CPROW,NBANDS, &
     &  NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL)
      USE prec

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      DIMENSION CPTWFP(NPLDIM,NBANDS)
      DIMENSION CFW(NPLDIM,NSTRIP)
      COMPLEX(q)      CPROJ(NPROD,NBANDS)
      COMPLEX(q)      CPROW(NPROD,NSTRIP)
      COMPLEX(q)      COVL(NBANDS,NBANDS)

      IF (NPL/=0) THEN
! here external blocking can be 1._q, but pretty useless
      NBLOCK= NPL
      DO NPOSPL=1, NPL-NBLOCK,NBLOCK
      CALL ZGEMM('C','N',NBANDS,NSTRIP,NBLOCK,(1._q,0._q), &
     &         CPTWFP(NPOSPL,1), NPLDIM,CFW(NPOSPL,1), &
     &          NPLDIM,(1._q,0._q),COVL(1,NPOS),NBANDS)
      ENDDO
      CALL ZGEMM('C','N',NBANDS,NSTRIP, NPL-NPOSPL+1,(1._q,0._q), &
     &         CPTWFP(NPOSPL,1), NPLDIM,CFW(NPOSPL,1), &
     &          NPLDIM,(1._q,0._q),COVL(1,NPOS),NBANDS)
      ENDIF

      IF (NPRO/=0) THEN
! here external blocking can be 1._q, but pretty useless
      NBLOCK=NPRO
      DO NPOSPR=1,NPRO-NBLOCK,NBLOCK
      CALL ZGEMM('C','N',NBANDS,NSTRIP,NBLOCK,(1._q,0._q), &
     &         CPROJ(NPOSPR,1), NPROD,CPROW(NPOSPR,1), &
     &          NPROD,(1._q,0._q),COVL(1,NPOS),NBANDS)
      ENDDO
      CALL ZGEMM('C','N',NBANDS,NSTRIP,NPRO-NPOSPR+1,(1._q,0._q), &
     &         CPROJ(NPOSPR,1), NPROD,CPROW(NPOSPR,1), &
     &          NPROD,(1._q,0._q),COVL(1,NPOS),NBANDS)
      ENDIF

      RETURN
    END SUBROUTINE ORTH2

!************************ SUBROUTINE ORTH3 *****************************
!> Similar to ::ORTH2 but allows for a different leading dimension NDIM
!> for the matrix to be set up (COVL)
!>
!> Explicit interface to this subroutine is in dfast::ORTH3 (beginning of
!> dfast module).
!***********************************************************************

    SUBROUTINE ORTH3(CPTWFP,CFW,CPROJ,CPROW, &
         & NBANDS,NSTRIP,NDIM,NPOS, &
         & NPL,NPRO,NPLDIM,NPROD,COVL)
      USE prec

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      DIMENSION CPTWFP(NPLDIM,NBANDS)
      DIMENSION CFW(NPLDIM,NSTRIP)
      COMPLEX(q)      CPROJ(NPROD,NBANDS)
      COMPLEX(q)      CPROW(NPROD,NSTRIP)
      COMPLEX(q)      COVL(NDIM,*)

      IF (NPL/=0) THEN
! here external blocking can be 1._q, but pretty useless
         NBLOCK= NPL
         DO NPOSPL=1, NPL-NBLOCK,NBLOCK
            CALL ZGEMM('C','N',NBANDS,NSTRIP,NBLOCK,(1._q,0._q), &
                 &         CPTWFP(NPOSPL,1), NPLDIM,CFW(NPOSPL,1), &
                 &          NPLDIM,(1._q,0._q),COVL(1,NPOS),NDIM)
         ENDDO
         CALL ZGEMM('C','N',NBANDS,NSTRIP, NPL-NPOSPL+1,(1._q,0._q), &
              &         CPTWFP(NPOSPL,1), NPLDIM,CFW(NPOSPL,1), &
              &          NPLDIM,(1._q,0._q),COVL(1,NPOS),NDIM)
      ENDIF
      
      IF (NPRO/=0) THEN
! here external blocking can be 1._q, but pretty useless
         NBLOCK=NPRO
         DO NPOSPR=1,NPRO-NBLOCK,NBLOCK
            CALL ZGEMM('C','N',NBANDS,NSTRIP,NBLOCK,(1._q,0._q), &
                 &         CPROJ(NPOSPR,1), NPROD,CPROW(NPOSPR,1), &
                 &          NPROD,(1._q,0._q),COVL(1,NPOS),NDIM)
         ENDDO
         CALL ZGEMM('C','N',NBANDS,NSTRIP,NPRO-NPOSPR+1,(1._q,0._q), &
              &         CPROJ(NPOSPR,1), NPROD,CPROW(NPOSPR,1), &
              &          NPROD,(1._q,0._q),COVL(1,NPOS),NDIM)
      ENDIF
      
      RETURN
    END SUBROUTINE ORTH3

!************************ SUBROUTINE ORTH1_NOSUBINDEX ******************
!> Variant of of ::ORTH1
!>
!> ~~~
!>   O(I,J) =  <CPTWFP(I) |  CFa(J) > +  <CPROJ(I) | CPROW(J) >
!>       J=NPOS ,..., NPOS+NSTRIP-1
!>       I=1 ,..., NBANDS
!> ~~~
!> This routine is for use in 1 aware versions of subrot.F and
!> choleski2.F.
!> The ::ORTH1 version accesses COVL at the storage position NPOS.
!> This version accesses COVL at storage position 1
!> regardless of NPOS.
!***********************************************************************

    SUBROUTINE ORTH1_NOSUBINDEX(CSEL,CPTWFP,CFW,CPROJ,CPROW,NBANDS, &
     &  NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL)
      USE prec
      USE string, ONLY: str
      USE tutor, ONLY: vtutor

      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      DIMENSION CPTWFP(NPLDIM,NBANDS)
      DIMENSION CFW(NPLDIM,NSTRIP)
      COMPLEX(q)      CPROJ(NPROD,NBANDS)
      COMPLEX(q)      CPROW(NPROD,NSTRIP)
      COMPLEX(q)      COVL(NBANDS,NBANDS)
      CHARACTER*(*) CSEL

      

      IF (NSTRIP+NPOS-1 > NBANDS) THEN
        CALL vtutor%bug("internal error in ORTH1_NOSUBINDEX: dim= " // str(NSTRIP+NPOS) // " " // &
           str(NBANDS), "dfast.F", 889)
      ENDIF
!
! update of lower triangular part
!
    IF (CSEL(1:1) == 'L' .OR. CSEL(1:1) == 'l') THEN
      IF (NPL/=0) THEN
      NBLOCK= NPL

      DO NPOSPL=1, NPL-NBLOCK,NBLOCK
      CALL ZGEMM('C','N',NBANDS-NPOS+1,NSTRIP,NBLOCK,(1._q,0._q), &
              CPTWFP(NPOSPL,NPOS), NPLDIM,CFW(NPOSPL,1), &
               NPLDIM,(1._q,0._q),COVL(NPOS,1),NBANDS)
      ENDDO
      CALL ZGEMM('C','N',NBANDS-NPOS+1,NSTRIP, NPL-NPOSPL+1,(1._q,0._q), &
              CPTWFP(NPOSPL,NPOS), NPLDIM,CFW(NPOSPL,1), &
               NPLDIM,(1._q,0._q),COVL(NPOS,1),NBANDS)
      ENDIF

      IF (NPRO/=0) THEN

      NBLOCK=NPRO
      DO NPOSPR=1,NPRO-NBLOCK,NBLOCK
      CALL ZGEMM('C','N',NBANDS-NPOS+1,NSTRIP,NBLOCK,(1._q,0._q), &
              CPROJ(NPOSPR,NPOS), NPROD,CPROW(NPOSPR,1), &
              NPROD,(1._q,0._q),COVL(NPOS,1),NBANDS)
      ENDDO
      CALL ZGEMM('C','N',NBANDS-NPOS+1,NSTRIP,NPRO-NPOSPR+1,(1._q,0._q), &
              CPROJ(NPOSPR,NPOS), NPROD,CPROW(NPOSPR,1), &
              NPROD,(1._q,0._q),COVL(NPOS,1),NBANDS)
      ENDIF
!
! update of upper triangular part
!
    ELSE IF (CSEL(1:1) == 'U' .OR. CSEL(1:1) == 'u') THEN
      IF (NPL/=0) THEN
      NBLOCK= NPL

      DO NPOSPL=1, NPL-NBLOCK,NBLOCK
      CALL ZGEMM('C','N',NSTRIP+NPOS-1,NSTRIP,NBLOCK,(1._q,0._q), &
              CPTWFP(NPOSPL,1), NPLDIM,CFW(NPOSPL,1), &
               NPLDIM,(1._q,0._q),COVL(1,1),NBANDS)
      ENDDO
      CALL ZGEMM('C','N',NSTRIP+NPOS-1,NSTRIP, NPL-NPOSPL+1,(1._q,0._q), &
              CPTWFP(NPOSPL,1), NPLDIM,CFW(NPOSPL,1), &
               NPLDIM,(1._q,0._q),COVL(1,1),NBANDS)
      ENDIF

      IF (NPRO/=0) THEN

      NBLOCK=NPRO
      DO NPOSPR=1,NPRO-NBLOCK,NBLOCK
      CALL ZGEMM('C','N',NSTRIP+NPOS-1,NSTRIP,NBLOCK,(1._q,0._q), &
              CPROJ(NPOSPR,1), NPROD,CPROW(NPOSPR,1), &
              NPROD,(1._q,0._q),COVL(1,1),NBANDS)
      ENDDO
      CALL ZGEMM('C','N',NSTRIP+NPOS-1,NSTRIP,NPRO-NPOSPR+1,(1._q,0._q), &
              CPROJ(NPOSPR,1), NPROD,CPROW(NPOSPR,1), &
              NPROD,(1._q,0._q),COVL(1,1),NBANDS)
      ENDIF

    ELSE
      WRITE(*,*)'internal error in ORTH1_NOSUBINDEX: CSEL=',CSEL
    ENDIF

    

    RETURN
    END SUBROUTINE ORTH1_NOSUBINDEX

!************************ SUBROUTINE ORTH2_NOSUBINDEX ******************
!> Variant of ::ORTH2
!>
!> ~~~
!>   O(I,J) =  <CPTWFP(I) |  CFa(J) > +  <CPROJ(I) | CPROW(J) >
!>       J=NPOS ,..., NPOS+NSTRIP-1
!>       I=1 ,..., NBANDS
!> ~~~
!> This routine is for use in 1 aware versions of subrot.F and
!> choleski2.F.
!> The ::ORTH2 version accesses COVL at the storage position NPOS.
!> This version accesses COVL at storage position 1
!> regardless of NPOS.
!***********************************************************************

    SUBROUTINE ORTH2_NOSUBINDEX(CPTWFP,CFW,CPROJ,CPROW,NBANDS, &
     &  NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL)
      USE prec

      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      DIMENSION CPTWFP(NPLDIM,NBANDS)
      DIMENSION CFW(NPLDIM,NSTRIP)
      COMPLEX(q)      CPROJ(NPROD,NBANDS)
      COMPLEX(q)      CPROW(NPROD,NSTRIP)
      COMPLEX(q)      COVL(NBANDS,NBANDS)
!
! update of upper triangular part
!
      IF (NPL/=0) THEN
      NBLOCK= NPL

      DO NPOSPL=1, NPL-NBLOCK,NBLOCK
      CALL ZGEMM('C','N',NBANDS,NSTRIP,NBLOCK,(1._q,0._q), &
              CPTWFP(NPOSPL,1), NPLDIM,CFW(NPOSPL,1), &
               NPLDIM,(1._q,0._q),COVL(1,1),NBANDS)
      ENDDO
      CALL ZGEMM('C','N',NBANDS,NSTRIP, NPL-NPOSPL+1,(1._q,0._q), &
              CPTWFP(NPOSPL,1), NPLDIM,CFW(NPOSPL,1), &
               NPLDIM,(1._q,0._q),COVL(1,1),NBANDS)
      ENDIF

      IF (NPRO/=0) THEN

      NBLOCK=NPRO
      DO NPOSPR=1,NPRO-NBLOCK,NBLOCK
      CALL ZGEMM('C','N',NBANDS,NSTRIP,NBLOCK,(1._q,0._q), &
              CPROJ(NPOSPR,1), NPROD,CPROW(NPOSPR,1), &
              NPROD,(1._q,0._q),COVL(1,1),NBANDS)
      ENDDO
      CALL ZGEMM('C','N',NBANDS,NSTRIP,NPRO-NPOSPR+1,(1._q,0._q), &
              CPROJ(NPOSPR,1), NPROD,CPROW(NPOSPR,1), &
              NPROD,(1._q,0._q),COVL(1,1),NBANDS)
      ENDIF

    RETURN
    END SUBROUTINE ORTH2_NOSUBINDEX

!************************ SUBROUTINE LINCOM_SLICE **********************
!
!> Calculates a linear combination
!>
!> Operates on a block of the matrix CF calculates a linear combination:
!> ~~~
!> COUT_n,k =   sum_kp CIN_n,kp CH_ kp,k
!> ~~~
!> This is the version required for 1 aware operation.
!
!***********************************************************************


    SUBROUTINE LINCOM_SLICE(MODE,CF,CPROF,CTRANS,NIN,NPOS,NOUT,NPL, &
     &           NPRO,NPLDIM,NPROD,LDTRAN,CF_RESULT,CPROF_RESULT)
      USE prec

      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      CHARACTER*1 MODE
      DIMENSION   CF(NPLDIM,NIN)
      DIMENSION   CF_RESULT(NPLDIM,NIN)
      COMPLEX(q)        CPROF(NPROD,NIN)
      COMPLEX(q)        CPROF_RESULT(NPROD,NIN)
      COMPLEX(q)        CTRANS(LDTRAN,NOUT)
      INTEGER     NIN,NPOS,NOUT,NPL,NPRO,NPLDIM,NPROD,LDTRAN

      

      CALL LINBAS_SLICE(MODE,CF,CF_RESULT,CTRANS,NIN,NPOS,NOUT, NPL, &
     &             NPLDIM,LDTRAN)
      IF (NPRO/=0) THEN
      CALL LINBAS_SLICE(MODE,CPROF,CPROF_RESULT,CTRANS,NIN,NPOS,NOUT, NPRO, &
     &             NPROD,LDTRAN)
      ENDIF

      

      RETURN
    END SUBROUTINE LINCOM_SLICE
      

    SUBROUTINE LINBAS_SLICE(MODE,CF,CF_RESULT,CTRANS,NIN,NPOS,NOUT,NPL, &
     &           NPLDIM,LDTRAN)

!! USE moffload_struct_def

      USE prec
      USE tutor, ONLY: vtutor

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      INTEGER NIN,NPOS,NOUT,NPL,NPLDIM,LDTRAN

      CHARACTER*1 MODE
      LOGICAL     LTRI
      COMPLEX(q)     CF(NPLDIM,NIN)
      COMPLEX(q)     CF_RESULT(NPLDIM,NIN)
      COMPLEX(q)     CTRANS(LDTRAN,NOUT)

      IF (NOUT>NIN) THEN
         CALL vtutor%bug("LINBAS_SLICE: wrong arguments, NOUT>NIN", "dfast.F", 1081)
      ENDIF

      LTRI =(MODE=='U').OR.(MODE=='u').OR. &
     &      (MODE=='L').OR.(MODE=='l')
! this part makes sure that the upper or lower triangle is cleared
! if we use DGEMM instead of ZTRMM
      IF (LTRI) THEN
         IF ((MODE=='L').OR.(MODE=='l')) THEN
            CALL vtutor%bug("LINBAS_SLICE: lower triangle is not yet supported", "dfast.F", 1090)
            DO N2=1,NOUT
!DIR$ IVDEP
!OCL NOVREC
               DO N1=1,N2+NPOS-2
                  CTRANS(N1,N2)= (0._q,0._q)
               ENDDO
            ENDDO
         ELSE
!$ACC PARALLEL LOOP PRESENT(CTRANS) PRIVATE(N1) 
            DO N2=1,NOUT
!DIR$ IVDEP
!OCL NOVREC
               DO N1=N2+NPOS,NIN
                  CTRANS(N1,N2)= (0._q,0._q)
               ENDDO
            ENDDO
         ENDIF
      ENDIF

      ILENPL=NPL

      IF (LTRI) THEN
! 'Triangular update':
         CALL ZGEMM('N', 'N', ILENPL, NOUT, NOUT+NPOS-1, (1._q,0._q), &
     &        CF(1,1), NPLDIM, CTRANS(1,1), &
     &        LDTRAN, (0._q,0._q), CF_RESULT(1,1), NPLDIM)
      ELSE
! 'Full update':
         CALL ZGEMM('N', 'N', ILENPL, NOUT, NIN, (1._q,0._q), &
     &        CF(1,1), NPLDIM, CTRANS(1,1), &
     &        LDTRAN, (0._q,0._q), CF_RESULT(1,1), NPLDIM)
      ENDIF

      RETURN
    END SUBROUTINE LINBAS_SLICE
