# 1 "hamil.F"
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


# 2 "hamil.F" 2 
!***********************************************************************
!
!> this module implements most of the low and high level
!> routines to calculated the action of a local Hamiltonian
!> onto the wavefunctions
!
!***********************************************************************
MODULE hamil

!! USE moffload

  USE prec
  USE wave_high

  INTERFACE
     SUBROUTINE PW_NORM_WITH_METRIC(WDES1, C, FNORM, FMETRIC, METRIC)
       USE prec
       USE wave
       REAL(q) FNORM
       TYPE (wavedes1)    WDES1
       COMPLEX(q) :: C
       REAL(q), OPTIONAL :: FMETRIC
       REAL(q), OPTIONAL :: METRIC(WDES1%NPL)
     END SUBROUTINE PW_NORM_WITH_METRIC
  END INTERFACE

  INTERFACE
     SUBROUTINE KINHAMIL( WDES1, GRID, CVR,  LADD, DATAKE, EVALUE, CPTWFP, CH)
       USE mgrid
       USE wave

       TYPE (wavedes1)    WDES1
       TYPE (grid_3d)     GRID
       COMPLEX(q) :: CVR
       LOGICAL    :: LADD
       REAL(q)    :: DATAKE
       REAL(q)    :: EVALUE
       COMPLEX(q) :: CPTWFP
       COMPLEX(q) :: CH
     END SUBROUTINE KINHAMIL
  END INTERFACE

  INTERFACE
     SUBROUTINE VHAMIL(WDES1,GRID,SV,CR,CVR)
       USE prec
       USE mgrid
       USE wave
       TYPE (grid_3d)     GRID
       TYPE (wavedes1)    WDES1
       COMPLEX(q)   SV
       COMPLEX(q) :: CR
       COMPLEX(q) :: CVR
     END SUBROUTINE VHAMIL
  END INTERFACE

  INTERFACE
     SUBROUTINE PW_CHARGE(WDES1,  CHARGE, NDIM, CR1, CR2, WEIGHT)
       USE prec
       USE mgrid
       USE wave
       TYPE (grid_3d)     GRID
       TYPE (wavedes1)    WDES1
       INTEGER NDIM
       COMPLEX(q)   CHARGE
       COMPLEX(q) :: CR1,CR2
       REAL(q) :: WEIGHT
     END SUBROUTINE PW_CHARGE
  END INTERFACE

  INTERFACE
     SUBROUTINE PW_CHARGE_CMPLX(WDES1,  CHARGE, NDIM, CR1, CR2)
       USE prec
       USE mgrid
       USE wave
       TYPE (grid_3d)     GRID
       TYPE (wavedes1)    WDES1
       INTEGER NDIM
       COMPLEX(q)   CHARGE
       COMPLEX(q) :: CR1,CR2
       REAL(q) :: WEIGHT
     END SUBROUTINE PW_CHARGE_CMPLX
  END INTERFACE

  INTERFACE
     SUBROUTINE VHAMIL_TRACE(WDES1, GRID, SV, CR, CVR, WEIGHT)
       USE prec
       USE mgrid
       USE wave

       TYPE (grid_3d)     GRID
       TYPE (wavedes1)    WDES1

       COMPLEX(q)   SV
       COMPLEX(q) :: CR,CVR
       REAL(q)    :: WEIGHT
     END SUBROUTINE VHAMIL_TRACE
  END INTERFACE

  INTERFACE
     SUBROUTINE PW_CHARGE_TRACE(WDES1,  CHARGE, CR1, CR2)
       USE prec
       USE mgrid
       USE wave
       TYPE (grid_3d)     GRID
       TYPE (wavedes1)    WDES1
       COMPLEX(q)   CHARGE
       COMPLEX(q) :: CR1,CR2
     END SUBROUTINE PW_CHARGE_TRACE
  END INTERFACE

  INTERFACE
     SUBROUTINE PW_CHARGE_TRACE_GDEF(WDES1,  CHARGE, CR1, CR2)
       USE prec
       USE mgrid
       USE wave
       TYPE (grid_3d)     GRID
       TYPE (wavedes1)    WDES1
       COMPLEX(q)   CHARGE
       COMPLEX(q) :: CR1,CR2
     END SUBROUTINE PW_CHARGE_TRACE_GDEF
  END INTERFACE

  INTERFACE
     SUBROUTINE PW_CHARGE_TRACE_NO_CONJG(WDES1,  CHARGE, CR1, CR2)
       USE prec
       USE mgrid
       USE wave
       TYPE (grid_3d)     GRID
       TYPE (wavedes1)    WDES1
       COMPLEX(q)   CHARGE
       COMPLEX(q) :: CR1,CR2
     END SUBROUTINE PW_CHARGE_TRACE_NO_CONJG
  END INTERFACE

  INTERFACE
     SUBROUTINE ECCP_NL_(LMDIM,LMMAXC,CDIJ,CQIJ,EVALUE,CPROJ1,CPROJ2,CNL)
       USE prec
       IMPLICIT NONE
       COMPLEX(q)  CNL
       INTEGER LMMAXC, LMDIM
       COMPLEX(q) CDIJ,CQIJ
       REAL(q) EVALUE
       COMPLEX(q) CPROJ1,CPROJ2
     END SUBROUTINE ECCP_NL_
  END INTERFACE

CONTAINS

!************************* SUBROUTINE ECCP   ***************************
! RCS:  $Id: hamil.F,v 1.3 2002/08/14 13:59:39 kresse Exp $
!
!> this subroutine calculates the expectation value of <c|H|cp>
!> where c and cp are two wavefunctions
!> FFT and non-local projections of wavefunctions st be supplied
!
!> @details @ref openmp :
!> loops over grid points and basis vectors are distributed over
!> the available OpenMP threads.
!> @n
!> A nested loop over \"types\" + \"ions-per-type\" is replaced
!> by a loop over \"all ions\", and is distributed over the
!> available OpenMP threads.
!
!***********************************************************************

  SUBROUTINE ECCP(WDES1,W1,W2,LMDIM,CDIJ,GRID,SV,CE)
    USE prec
    USE mpimy
    USE mgrid
    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (wavefun1) :: W1,W2
    TYPE (wavedes1) :: WDES1
    TYPE (grid_3d)  :: GRID

    INTEGER NGVECTOR, ISPINOR, ISPINOR_
    INTEGER NPRO,NPRO_

    COMPLEX(q)    CNL
    COMPLEX(q)   SV(GRID%MPLWV,WDES1%NRSPINORS*WDES1%NRSPINORS) !< local potential
    COMPLEX(q) CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
!$ACC ROUTINE (ECCP_NL) VECTOR

    

!=======================================================================
! calculate the local contribution
!=======================================================================
    NGVECTOR=WDES1%NGVECTOR

!$ACC ENTER DATA CREATE(CE) 
!$ACC KERNELS PRESENT(CE) 
    CE=0
!$ACC END KERNELS

!$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(GRID,W1,W2,WDES1,SV,CE) &
!$ACC PRIVATE(MM,MM_) REDUCTION(+:CE) 

!$OMP PARALLEL &
!$OMP SHARED(WDES1,GRID,SV,W1,W2,NGVECTOR,LMDIM,CDIJ) &
!$OMP PRIVATE(ISPINOR,ISPINOR_,M,MM,MM_,NPRO,NPRO_,NI,NT,LMMAXC,CNL) &
!$OMP REDUCTION(+:CE)

    DO ISPINOR =0,WDES1%NRSPINORS-1
       DO ISPINOR_=0,WDES1%NRSPINORS-1
 !$OMP DO SCHEDULE(STATIC)
          DO M=1,GRID%RL%NP
             MM =M+ISPINOR *GRID%MPLWV
             MM_=M+ISPINOR_*GRID%MPLWV
             CE=CE+SV(M,1+ISPINOR_+2*ISPINOR) *W1%CR(MM_)*CONJG(W2%CR(MM))
          ENDDO
 !$OMP END DO
       ENDDO
    ENDDO

!$ACC KERNELS PRESENT(CE,GRID) 
    CE=CE/GRID%NPLWV
!$ACC END KERNELS
!=======================================================================
! kinetic energy contribution
!=======================================================================
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(W1,W2,WDES1,CE) &
!$ACC PRIVATE(MM) REDUCTION(+:CE) 
    DO ISPINOR=0,WDES1%NRSPINORS-1
 !$OMP DO SCHEDULE(STATIC)
       DO M=1,NGVECTOR
          MM=M+ISPINOR*NGVECTOR
          CE=CE+W1%CPTWFP(MM)*CONJG(W2%CPTWFP(MM))*WDES1%DATAKE(M,ISPINOR+1)
       ENDDO
 !$OMP END DO
    ENDDO
!=======================================================================
! non local contribution
!=======================================================================
!$ACC PARALLEL LOOP COLLAPSE(3) GANG PRESENT(W1,W2,CDIJ,WDES1,CE) &
!$ACC PRIVATE(NPRO,NPRO_,NT,LMMAXC,CNL) REDUCTION(+:CE) 
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
       DO ISPINOR_=0,WDES1%NRSPINORS-1
 !$OMP DO SCHEDULE(STATIC)
          DO NI=1,WDES1%NIONS
             NT=WDES1%ITYP(NI)
             LMMAXC=WDES1%LMMAX(NT)
             NPRO =ISPINOR *(WDES1%NPRO/2)+WDES1%LMBASE(NI)
             NPRO_=ISPINOR_*(WDES1%NPRO/2)+WDES1%LMBASE(NI)
             IF (LMMAXC==0) CYCLE
             CNL=0; CALL ECCP_NL(LMDIM,LMMAXC,CDIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),W1%CPROJ(NPRO_+1),W2%CPROJ(NPRO+1),CNL)
             CE=CE+CNL
          ENDDO
 !$OMP END DO
       ENDDO
    ENDDO spinor
 !$OMP END PARALLEL

!$ACC KERNELS PRESENT(CE) 
    CE=(CE)
!$ACC END KERNELS
    CALL M_sum_z(WDES1%COMM_INB, CE, 1)
!$ACC UPDATE SELF(CE) WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
!$ACC EXIT DATA DELETE(CE) 

    

  END SUBROUTINE ECCP

!************************* SUBROUTINE ECCP_VEC *************************
! RCS:  $Id: hamil.F,v 1.3 2002/08/14 13:59:39 kresse Exp $
!
!> this subroutine calculates the expectation value of <c|H|cp>
!> where c and cp are two wavefunctions
!> FFT and non-local projections of wavefunctions must be supplied
!
!> @details @ref openmp :
!> loops over grid points and basis vectors are distributed over
!> the available OpenMP threads.
!> @n
!> A nested loop over \"types\" + \"ions-per-type\" is replaced
!> by a loop over \"all ions\", and is distributed over the
!> available OpenMP threads.
!> @n
!> Furthermore the workarrays CWORK1 and CWORK2 acquire a second
!> dimension (size=3).  This additional workspace improves the
!> performance of ::kinhamil_vec under OpenMP parallelization.
!
!***********************************************************************

  SUBROUTINE ECCP_VEC(WDES1,W1,W2,LMDIM,CDIJ,GRID,SV,AVEC, CE)
    USE prec
    USE mpimy
    USE mgrid
    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)


    TYPE (wavefun1) :: W1,W2
    TYPE (wavedes1) :: WDES1
    TYPE (grid_3d)  :: GRID

    INTEGER NGVECTOR, ISPINOR
    INTEGER NPRO,NPRO_
!$  INTEGER NPRO2,NPRO2_

    COMPLEX(q)      CNL
    COMPLEX(q)   SV(GRID%MPLWV,WDES1%NRSPINORS*WDES1%NRSPINORS)
    COMPLEX(q)   AVEC(:,:)

    COMPLEX(q) :: CVR(GRID%MPLWV*WDES1%NRSPINORS)
    COMPLEX(q) :: CWORK1(WDES1%NRPLWV     )
    COMPLEX(q) :: CWORK2(WDES1%GRID%MPLWV )
# 313


    COMPLEX(q) CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
!$ACC ROUTINE (ECCP_NL) VECTOR

    

!=======================================================================
! calculate the local contribution
!=======================================================================
    NGVECTOR=WDES1%NGVECTOR

!$ACC ENTER DATA CREATE(CLOCAL) 

!$ACC KERNELS PRESENT(CLOCAL) 
    CLOCAL=0
!$ACC END KERNELS

!$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(GRID,WDES1,SV,W1,W2,CLOCAL) &
!$ACC PRIVATE(MM,MM_) REDUCTION(+:CLOCAL) 

!$OMP PARALLEL SHARED(WDES1,GRID,SV,W1,W2) &
!$OMP PRIVATE(ISPINOR,ISPINOR_,M,MM,MM_) &
!$OMP REDUCTION(+:CLOCAL)

    DO ISPINOR =0,WDES1%NRSPINORS-1
       DO ISPINOR_=0,WDES1%NRSPINORS-1
 !$OMP DO SCHEDULE(STATIC)
          DO M=1,GRID%RL%NP
             MM =M+ISPINOR *GRID%MPLWV
             MM_=M+ISPINOR_*GRID%MPLWV
             CLOCAL=CLOCAL+SV(M,1+ISPINOR_+2*ISPINOR) *W1%CR(MM_)*CONJG(W2%CR(MM))
          ENDDO
 !$OMP END DO
       ENDDO
    ENDDO
 !$OMP END PARALLEL

!$ACC KERNELS PRESENT(CLOCAL,GRID) 
    CLOCAL=CLOCAL/GRID%NPLWV
!$ACC END KERNELS
!$ACC EXIT DATA COPYOUT(CLOCAL) 

!=======================================================================
! kinetic energy contribution
!=======================================================================
!$ACC ENTER DATA CREATE(CVR,CWORK1,CWORK2) 

!$ACC KERNELS PRESENT(CKIN,CVR) 
    CKIN=0; CVR=0
!$ACC END KERNELS
    CALL KINHAMIL_VEC( WDES1, GRID, CVR,  .FALSE., & 
         WDES1%DATAKE(1,1), WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), AVEC(1,1), &
         CWORK1(1 ), CWORK2(1 ), 0.0_q, W1%CPTWFP(1), CWORK1(1 ))

!$ACC EXIT DATA DELETE(CVR,CWORK2) 

!$ACC ENTER DATA CREATE(CTMP) 
!$ACC KERNELS PRESENT(CTMP) 
    CTMP=0
!$ACC END KERNELS

!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CKIN,CWORK1,W2) &
!$ACC PRIVATE(MM) REDUCTION(+:CKIN) 

!$OMP PARALLEL &
!$OMP SHARED(WDES1,NGVECTOR,CWORK1,W1,W2,LMDIM,CDIJ) &
!$OMP PRIVATE(ISPINOR,ISPINOR_,M,MM,NRPO,NPRO_,NI,NT,LMMAXC,CNL) &
!$OMP REDUCTION(+:CKIN,CTMP)

    DO ISPINOR=0,WDES1%NRSPINORS-1
 !$OMP DO SCHEDULE(STATIC)
       DO M=1,NGVECTOR
          MM=M+ISPINOR*NGVECTOR
          CKIN=CKIN+CWORK1(MM )*CONJG(W2%CPTWFP(MM))
       ENDDO
 !$OMP END DO
    ENDDO
!$ACC EXIT DATA COPYOUT(CKIN) 
!$ACC EXIT DATA DELETE(CWORK1) 

!=======================================================================
! non local contribution
!=======================================================================
!$ACC PARALLEL LOOP COLLAPSE(3) GANG PRESENT(W1,W2,CDIJ,WDES1,CTMP) &
!$ACC PRIVATE(NPRO,NPRO_,NT,LMMAXC,CNL) REDUCTION(+:CTMP) 
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
       DO ISPINOR_=0,WDES1%NRSPINORS-1
 !$OMP DO SCHEDULE(STATIC)
          DO NI=1,WDES1%NIONS
             NT=WDES1%ITYP(NI)
             LMMAXC=WDES1%LMMAX(NT)
             NPRO =ISPINOR *(WDES1%NPRO/2)+WDES1%LMBASE(NI)
             NPRO_=ISPINOR_*(WDES1%NPRO/2)+WDES1%LMBASE(NI)
             IF (LMMAXC==0) CYCLE
             CNL=0; CALL ECCP_NL(LMDIM,LMMAXC,CDIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),W1%CPROJ(NPRO_+1),W2%CPROJ(NPRO+1),CNL)
             CTMP=CTMP+CNL
          ENDDO
 !$OMP END DO
       ENDDO
    ENDDO spinor
 !$OMP END PARALLEL

!$ACC EXIT DATA COPYOUT(CTMP) 
!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)

    CE=(CLOCAL+CKIN+CTMP)
    CALL M_sum_z(WDES1%COMM_INB, CE, 1)

    

  END SUBROUTINE ECCP_VEC


!************************* SUBROUTINE ECCP_TAU *************************
! RCS:  $Id: hamil.F,v 1.3 2002/08/14 13:59:39 kresse Exp $
!
!> this subroutine calculates the expectation value of <c|H|cp>
!> where c and cp are two wavefunctions
!> FFT and non-local projections of wavefunctions must be supplied
!
!> @details @ref openmp :
!> loops over grid points and basis vectors are distributed over
!> the available OpenMP threads.
!> @n
!> A nested loop over \"types\" + \"ions-per-type\" is replaced
!> by a loop over \"all ions\", and is distributed over the
!> available OpenMP threads.
!> @n
!> Furthermore the workarrays CWORK1 and CWORK2 acquire a second
!> dimension (size=3 and 6, respectively). This additional workspace
!> improves the performance of ::kinhamil_tau under OpenMP
!> parallelization.
!
!***********************************************************************

  SUBROUTINE ECCP_TAU(WDES1,W1,W2,LMDIM,CDIJ,GRID,SV,LATT_CUR,MU,CE)
    USE prec
    USE mpimy
    USE mgrid
    USE lattice
    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (wavefun1) :: W1,W2
    TYPE (wavedes1) :: WDES1
    TYPE (grid_3d)  :: GRID
    TYPE (latt)     :: LATT_CUR

    INTEGER NGVECTOR, ISPINOR
    INTEGER NPRO,NPRO_
!$  INTEGER NPRO2,NPRO2_

    COMPLEX(q)    CNL
    COMPLEX(q)   SV(GRID%MPLWV,WDES1%NRSPINORS*WDES1%NRSPINORS)
    COMPLEX(q)   MU(GRID%MPLWV,WDES1%NRSPINORS*WDES1%NRSPINORS)

    COMPLEX(q) :: CVR(GRID%MPLWV*WDES1%NRSPINORS)
    COMPLEX(q) :: CH(WDES1%NRPLWV)
    COMPLEX(q) :: CWORK1(WDES1%NRPLWV     )
    COMPLEX(q) :: CWORK2(WDES1%GRID%MPLWV )
# 476


    COMPLEX(q) CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
!$ACC ROUTINE (ECCP_NL) VECTOR

    

!=======================================================================
! calculate the local contribution
!=======================================================================
    NGVECTOR=WDES1%NGVECTOR

!$ACC ENTER DATA CREATE(CE,CLOCAL,CKIN,CTMP) 

!$ACC KERNELS PRESENT(CLOCAL) 
    CLOCAL=0
!$ACC END KERNELS

!$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(GRID,WDES1,SV,W1,W2,CLOCAL) &
!$ACC PRIVATE(MM,MM_) REDUCTION(+:CLOCAL) 

!$OMP PARALLEL SHARED(WDES1,GRID,SV,W1,W2) &
!$OMP PRIVATE(ISPINOR,ISPINOR_,M,MM,MM_) &
!$OMP REDUCTION(+:CLOCAL)

    DO ISPINOR =0,WDES1%NRSPINORS-1
       DO ISPINOR_=0,WDES1%NRSPINORS-1
 !$OMP DO SCHEDULE(STATIC)
          DO M=1,GRID%RL%NP
             MM =M+ISPINOR *GRID%MPLWV
             MM_=M+ISPINOR_*GRID%MPLWV
             CLOCAL=CLOCAL+SV(M,1+ISPINOR_+2*ISPINOR) *W1%CR(MM_)*CONJG(W2%CR(MM))
          ENDDO
 !$OMP END DO
       ENDDO
    ENDDO
 !$OMP END PARALLEL

!$ACC KERNELS PRESENT(CLOCAL,GRID) 
    CLOCAL=CLOCAL/GRID%NPLWV
!$ACC END KERNELS

!=======================================================================
! kinetic energy contribution
!=======================================================================
!$ACC ENTER DATA COPYIN(LATT_CUR,LATT_CUR%B) 
!$ACC ENTER DATA CREATE(CVR,CH,CWORK1,CWORK2) 

!$ACC KERNELS PRESENT(CVR,CH,CKIN) 
    CKIN=0; CVR=0; CH=0
!$ACC END KERNELS

    CALL KINHAMIL_TAU( WDES1, GRID, CVR,  .FALSE., .TRUE., & 
         WDES1%DATAKE(1,1), WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), LATT_CUR, MU(1,1), &
         CWORK1(1 ), CWORK2(1 ), 0.0_q, W1%CPTWFP(1), CH(1))

!$ACC EXIT DATA DELETE(LATT_CUR%B,LATT_CUR,CVR,CWORK1,CWORK2) 

!$ACC KERNELS PRESENT(CTMP) 
    CTMP=0
!$ACC END KERNELS

!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(WDES1,W2,CH,CKIN) &
!$ACC PRIVATE(MM) REDUCTION(+:CKIN) 

!$OMP PARALLEL &
!$OMP SHARED(WDES1,NGVECTOR,CH,W1,W2,LMDIM,CDIJ) &
!$OMP PRIVATE(ISPINOR,ISPINOR_,M,MM,NRPO,NPRO_,NI,NT,LMMAXC,CNL) &
!$OMP REDUCTION(+:CKIN,CTMP)

    DO ISPINOR=0,WDES1%NRSPINORS-1
 !$OMP DO SCHEDULE(STATIC)
       DO M=1,NGVECTOR
          MM=M+ISPINOR*NGVECTOR
          CKIN=CKIN+CH(MM)*CONJG(W2%CPTWFP(MM))
       ENDDO
 !$OMP END DO
    ENDDO
!$ACC EXIT DATA DELETE(CH) 

!=======================================================================
! non local contribution
!=======================================================================
!$ACC PARALLEL LOOP COLLAPSE(3) GANG PRESENT(W1,W2,CDIJ,WDES1,CTMP) &
!$ACC PRIVATE(NPRO,NPRO_,NT,LMMAXC,CNL) REDUCTION(+:CTMP) 
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
       DO ISPINOR_=0,WDES1%NRSPINORS-1
 !$OMP DO SCHEDULE(STATIC)
          DO NI=1,WDES1%NIONS
             NT=WDES1%ITYP(NI)
             LMMAXC=WDES1%LMMAX(NT)
             NPRO =ISPINOR *(WDES1%NPRO/2)+WDES1%LMBASE(NI)
             NPRO_=ISPINOR_*(WDES1%NPRO/2)+WDES1%LMBASE(NI)
             IF (LMMAXC==0) CYCLE
             CNL=0; CALL ECCP_NL(LMDIM,LMMAXC,CDIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),W1%CPROJ(NPRO_+1),W2%CPROJ(NPRO+1),CNL)
             CTMP=CTMP+CNL
          ENDDO
 !$OMP END DO
       ENDDO
    ENDDO spinor
 !$OMP END PARALLEL


!$ACC KERNELS PRESENT(CE,CLOCAL,CKIN,CTMP) 
    CE=(CLOCAL+CKIN+CTMP)
!$ACC END KERNELS
!$ACC EXIT DATA DELETE(CLOCAL,CKIN,CTMP) 
    CALL M_sum_z(WDES1%COMM_INB, CE, 1)
!$ACC UPDATE SELF(CE) WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
!$ACC EXIT DATA DELETE(CE) 

    

  END SUBROUTINE ECCP_TAU


!*********************************************************************
!
!> calculate the band structure energy
!
!*********************************************************************

  FUNCTION BANDSTRUCTURE_ENERGY(WDES, W) RESULT(E)
    IMPLICIT NONE
    TYPE (wavedes)     WDES
    TYPE (wavespin)    W      ! wavefunction
    REAL(q)            E
    INTEGER            ISP, NK, NB

    E=0
    DO ISP=1,WDES%ISPIN
       DO NK=1,WDES%NKPTS

          IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE

          DO NB=1,WDES%NB_TOT
             E=E+WDES%RSPIN* REAL( W%CELTOT(NB,NK,ISP) ,KIND=q) *WDES%WTKPT(NK)*W%FERTOT(NB,NK,ISP)
          ENDDO
       ENDDO
    ENDDO
    CALL M_sum_d(WDES%COMM_KINTER, E, 1)
  END FUNCTION BANDSTRUCTURE_ENERGY


!**********************************************************************
!
!> calculate the kinetic energy of each wavefunction
!
!**********************************************************************

  SUBROUTINE KINETIC_ENERGY(W )
    USE wave_high
    IMPLICIT NONE
    TYPE (wavespin), TARGET :: W
!    local
    TYPE (wavedes1)    WDES1
    TYPE (wavedes), POINTER :: WDES
    
    INTEGER ISP, NK, NB, ISPINOR, MM, M
    REAL(q) :: CKIN
    
    
    WDES=>W%WDES
    
    spin:  DO ISP=1,WDES%ISPIN
       kpoint: DO NK=1,WDES%NKPTS


          IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) THEN
!PK MRG_AUX already zeros contributions from non-kpoint nodes, but be cautious and (0._q,0._q) workspace anyway
             DO NB=1,WDES%NBANDS
                W%AUX(NB,NK,ISP)=0.
             END DO
          ELSE

          CALL SETWDES(WDES,WDES1,NK)
          
          DO NB=1,WDES%NBANDS
             CKIN=0
             DO ISPINOR=0,WDES1%NRSPINORS-1
                DO M=1,WDES1%NGVECTOR
                   MM=M+ISPINOR*WDES1%NGVECTOR
                   CKIN=CKIN+W%CPTWFP(MM,NB,NK,ISP)*CONJG(W%CPTWFP(MM,NB,NK,ISP))*WDES1%DATAKE(M,ISPINOR+1)
                ENDDO
             ENDDO
             W%AUX(NB,NK,ISP)=CKIN
          ENDDO

          ENDIF

       END DO kpoint
    END DO spin
    
    CALL MRG_AUX(WDES,W)
  END SUBROUTINE KINETIC_ENERGY
  
!************************* SUBROUTINE SETUP_PRECOND ********************
!
!> subroutine to set up a preconditioning matrix
!> IALGO determines the type of preconditioning
!>
!>-   0,6,8 TAP preconditioning
!>-    9    Jacobi like preconditioning
!>-   else  no preconditioning
!>
!>- the Jacobi like preconditer required the eigenvalue and
!>- a complex shift
!
!***********************************************************************


  SUBROUTINE SETUP_PRECOND( W1, IALGO, IDUMP, PRECON, EVALUE, DE_ATT )
    IMPLICIT NONE
    TYPE (wavefun1)    W1
    INTEGER IALGO         !< chosen algorithm
    INTEGER IDUMP         !< dump flag
    REAL(q) :: EVALUE     !< eigenvalue minus average local potential
    REAL(q) :: DE_ATT     !< complex shift
    REAL(q) :: PRECON(*)  !< return: preconditioning matrix
! local
    INTEGER M, MM, NK, ISPINOR
    REAL(q) :: EKIN, FAKT, X, X2
    COMPLEX(q) :: CPT

    
!$ACC ENTER DATA COPYIN(W1,W1%CPTWFP,W1%WDES1) 

    IF (IALGO==0 .OR. IALGO==8 .OR. IALGO==6) THEN
       EKIN=0
!$ACC PARALLEL LOOP COLLAPSE(2) GANG VECTOR REDUCTION(+:EKIN) PRIVATE(MM,CPT) &
!$ACC PRESENT(W1) 
       DO ISPINOR=0,W1%WDES1%NRSPINORS-1
!DIR$ IVDEP
!OCL NOVREL
 !$OMP PARALLEL DO PRIVATE(M,MM,CPT) SHARED(W1,ISPINOR) REDUCTION(+:EKIN)
          DO M=1,W1%WDES1%NGVECTOR
             MM=M+ISPINOR*W1%WDES1%NGVECTOR
             CPT=W1%CPTWFP(MM)
             EKIN =EKIN+ REAL( CPT*CONJG(CPT) ,KIND=q) * W1%WDES1%DATAKE(M,ISPINOR+1)
          ENDDO
 !$OMP END PARALLEL DO
       ENDDO
!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
       CALL M_sum_d(W1%WDES1%COMM_INB, EKIN, 1)

       IF (EKIN<2.0_q) EKIN=2.0_q
       EKIN=EKIN*1.5_q
       IF (IDUMP==2)  WRITE(*,'(E9.2,"K")',ADVANCE='NO') EKIN

       FAKT=2._q/EKIN
!$ACC PARALLEL LOOP COLLAPSE(2) GANG VECTOR PRIVATE(MM,X,X2) &
!$ACC PRESENT(PRECON,W1) 
       DO ISPINOR=0,W1%WDES1%NRSPINORS-1
!DIR$ IVDEP
!OCL NOVREL
 !$OMP PARALLEL DO PRIVATE(M,MM,X,X2) SHARED(W1,ISPINOR,EKIN,FAKT,PRECON)
          DO M=1,W1%WDES1%NGVECTOR
             MM=M+ISPINOR*W1%WDES1%NGVECTOR
             X=W1%WDES1%DATAKE(M,ISPINOR+1)/EKIN
             X2= 27+X*(18+X*(12+8*X))
             PRECON(MM)=X2/(X2+16*X*X*X*X)*FAKT
          ENDDO
 !$OMP END PARALLEL DO
       ENDDO
    ELSE IF (IALGO==9) THEN
!$ACC PARALLEL LOOP COLLAPSE(2) GANG VECTOR PRIVATE(MM,X) &
!$ACC PRESENT(PRECON,W1) 
       DO ISPINOR=0,W1%WDES1%NRSPINORS-1
!DIR$ IVDEP
!OCL NOVREL
 !$OMP PARALLEL DO PRIVATE(M,MM,X) SHARED(W1,EVALUE,DE_ATT,PRECON)
          DO M=1,W1%WDES1%NGVECTOR
             MM=M+ISPINOR*W1%WDES1%NGVECTOR
             X=MAX(W1%WDES1%DATAKE(M,ISPINOR+1)-EVALUE,0._q)
             PRECON(MM)= REAL( 1._q/(X+ CMPLX( 0._q , DE_ATT ,KIND=q) ) ,KIND=q) !new
          ENDDO
 !$OMP END PARALLEL DO
       ENDDO
    ELSE
!$ACC PARALLEL LOOP COLLAPSE(2) GANG VECTOR PRIVATE(MM) &
!$ACC PRESENT(PRECON,W1) 
       DO ISPINOR=0,W1%WDES1%NRSPINORS-1
!DIR$ IVDEP
!OCL NOVREL
 !$OMP PARALLEL DO PRIVATE(M,MM) SHARED(W1,PRECON)
          DO M=1,W1%WDES1%NGVECTOR
             MM=M+ISPINOR*W1%WDES1%NGVECTOR
             PRECON(MM)=1.0_q
          ENDDO
 !$OMP END PARALLEL DO
       ENDDO
    ENDIF

!$ACC EXIT DATA DELETE(W1,W1%CPTWFP,W1%WDES1) 
    

  END SUBROUTINE SETUP_PRECOND

!************************* SUBROUTINE APPLY_PRECOND ********************
!
!> apply W2%CPTWFP <- W1%CPTWFP * PRECON * MUL
!
!***********************************************************************

  SUBROUTINE APPLY_PRECOND( W1, W2, PRECON, MUL)
    USE gridq
    IMPLICIT NONE
    TYPE (wavefun1)    W1, W2
    REAL(q) :: PRECON(*)
    REAL(q), OPTIONAL :: MUL
    REAL(q) :: MUL_

    

    IF (PRESENT(MUL)) THEN
       MUL_=MUL
    ELSE
       MUL_=1._q
    ENDIF

    CALL REAL_CMPLX_CMPLX_MUL( W1%WDES1%NPL,  PRECON(1), W1%CPTWFP(1) , MUL_, W2%CPTWFP(1), 0.0_q, W2%CPTWFP(1))

    
    
  END SUBROUTINE APPLY_PRECOND

!************************* SUBROUTINE ADD_PRECOND **********************
!
!> apply W2%CPTWFP < W2%CPTWFP + W1%CPTWFP * PRECON * MUL
!
!***********************************************************************

  SUBROUTINE ADD_PRECOND( W1, W2, PRECON, MUL)
    USE gridq
    IMPLICIT NONE
    TYPE (wavefun1)    W1, W2
    REAL(q) :: PRECON(*)
    REAL(q), OPTIONAL :: MUL
    REAL(q) :: MUL_

    

    IF (PRESENT(MUL)) THEN
       MUL_=MUL
    ELSE
       MUL_=1.0
    ENDIF

    CALL REAL_CMPLX_CMPLX_MUL( W1%WDES1%NPL,  PRECON(1), W1%CPTWFP(1) , MUL_, W2%CPTWFP(1), 1.0_q, W2%CPTWFP(1))

    
    
  END SUBROUTINE ADD_PRECOND


  SUBROUTINE SIMPLE_PRECOND( W1 )
    USE gridq
    IMPLICIT NONE
    TYPE (wavefun1)    W1

    REAL(q) :: EKIN, FAKT, X, X2
    INTEGER :: ISPINOR, M, MM

! old version
    EKIN=10
    FAKT=2._q/EKIN
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(W1) PRIVATE(MM,X,X2) 
    DO ISPINOR=0,W1%WDES1%NRSPINORS-1
 !$OMP PARALLEL DO &
 !$OMP SHARED(W1,ISPINOR,EKIN,FAKT) &
 !$OMP PRIVATE(M,MM,X,X2)
       DO M=1,W1%WDES1%NGVECTOR
          MM=M+ISPINOR*W1%WDES1%NGVECTOR
          X=W1%WDES1%DATAKE(M,ISPINOR+1)/EKIN
          X2=27+X*(18+X*(12+8*X))
          W1%CPTWFP(MM)= W1%CPTWFP(MM)*X2/(X2+16*X*X*X*X)*FAKT
       ENDDO
 !$OMP END PARALLEL DO
    ENDDO

  END SUBROUTINE SIMPLE_PRECOND


!***********************************************************************
!
!> truncate high frequency components when LDELAY is set
!> and/or in case of spin-spiral calculations (LSPIRAL=.TRUE.)
!
!***********************************************************************

  SUBROUTINE TRUNCATE_HIGH_FREQUENCY_W1( W1, LDELAY, ENINI)
    USE prec
    IMPLICIT NONE
    LOGICAL LDELAY     !< usually truncation only during delay phase (LDELAY set)
    REAL(q) ENINI      !< cutoff at which the wavefunction is truncated
    TYPE (wavefun1)    W1
! local
    INTEGER ISPINOR

    

    IF (LDELAY.OR.W1%WDES1%LSPIRAL) THEN
       DO ISPINOR=1,W1%WDES1%NRSPINORS
          CALL TRUNCATE_HIGH_FREQUENCY_ONE(W1%WDES1%NGVECTOR, &
               W1%CPTWFP(1+(ISPINOR-1)*W1%WDES1%NGVECTOR), W1%WDES1%DATAKE(1,ISPINOR), ENINI)
       ENDDO
    END IF

    

  END SUBROUTINE TRUNCATE_HIGH_FREQUENCY_W1


!***********************************************************************
!
!> determine norm of a wavefunction and norm with a metric
!
!***********************************************************************

  SUBROUTINE PW_NORM_WITH_METRIC_W1(W1, FNORM, FMETRIC, METRIC)
    IMPLICIT NONE
    REAL(q) FNORM
    TYPE (wavefun1)    W1
    REAL(q), OPTIONAL :: FMETRIC
    REAL(q), OPTIONAL :: METRIC(W1%WDES1%NPL)

    IF (PRESENT(FMETRIC) .AND. PRESENT(METRIC)) THEN
       CALL PW_NORM_WITH_METRIC( W1%WDES1, W1%CPTWFP(1), FNORM, FMETRIC, METRIC)
    ELSE
       CALL PW_NORM_WITH_METRIC( W1%WDES1, W1%CPTWFP(1), FNORM)
    ENDIF
  END SUBROUTINE PW_NORM_WITH_METRIC_W1



!************************* SUBROUTINE HAMILT ***************************
!
!> this subroutine calculates the H acting onto a wavefuntion
!> the  wavefunction must be given in reciprocal space C and real
!> space CR
!> CH contains the result
!>   H- EVALUE Q |phi>
!> where Q is the nondiagonal part of the overlap matrix S
!***********************************************************************

  SUBROUTINE HAMILT( &
       &    W1, NONLR_S, NONL_S, EVALUE, CDIJ, CQIJ, SV, ISP, CH)
    USE mpimy
    USE mgrid
    
    USE nonl_high
    IMPLICIT NONE

    TYPE (wavefun1)    W1
    TYPE (nonlr_struct)NONLR_S
    TYPE (nonl_struct) NONL_S

    COMPLEX(q) CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    COMPLEX(q)      SV(:,:)
    INTEGER    ISP
    COMPLEX(q) CH(:)
    REAL(q)    EVALUE
! local variables
    COMPLEX(q) :: CWORK1(W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS)
# 942


    IF (NONLR_S%LREAL) THEN
! calculate the local contribution (result in CWORK1)
       CALL VHAMIL(W1%WDES1,W1%WDES1%GRID,SV(1,ISP),W1%CR(1),CWORK1(1)) 
! non-local contribution in real-space

       CALL RACC(NONLR_S, W1, CDIJ, CQIJ, ISP, EVALUE, CWORK1)

       CALL KINHAMIL( W1%WDES1, W1%WDES1%GRID, CWORK1(1), .FALSE., &
            W1%WDES1%DATAKE(1,1), 0.0_q, W1%CPTWFP(1), CH(1))
    ELSE
! calculate the local contribution (result in CWORK1)
       CALL VHAMIL(W1%WDES1,W1%WDES1%GRID,SV(1,ISP),W1%CR(1),CWORK1(1)) 

! calculate the non local contribution in reciprocal space
       CALL VNLACC(NONL_S, W1, CDIJ, CQIJ, ISP, EVALUE, CH)
       CALL KINHAMIL( W1%WDES1, W1%WDES1%GRID, CWORK1(1), .TRUE., &
            W1%WDES1%DATAKE(1,1), 0.0_q, W1%CPTWFP(1), CH(1))

    ENDIF
    RETURN
  END SUBROUTINE HAMILT

!************************* SUBROUTINE HAMILTMU *************************
!
!> this subroutine calculates the H acting onto a set of wavefuntions
!>
!> the  wavefunction must be given in reciprocal space C and real
!> space CR
!> CH contains the result
!>   H- EVALUE S |phi>
!> where S is the overlap matrix S
!
!***********************************************************************

  SUBROUTINE HAMILTMU( &
       &    WDES1, W1, NONLR_S, NONL_S, EVALUE, &
       &    CDIJ, CQIJ, SV, ISP, WRESULT)
    USE mpimy
    USE mgrid
    USE nonl_high
    IMPLICIT NONE

    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1(:)
    TYPE (nonlr_struct)NONLR_S
    TYPE (nonl_struct) NONL_S    
    REAL(q)    EVALUE(:)               ! eigenvalues
    COMPLEX(q) CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    COMPLEX(q)      SV(:,:)
    INTEGER    ISP
    TYPE(wavefuna)     WRESULT
    
! local variables
    COMPLEX(q) :: CWORK1(WDES1%GRID%MPLWV*WDES1%NRSPINORS,SIZE(W1))
    INTEGER NP,NSTRIP,NSTRIP_ACT
# 1001


    

!$ACC ENTER DATA CREATE(CWORK1) IF(OFFLOAD_ON)

! calculate the local contribution (result in CWORK1)
    DO NP=1,SIZE(W1)
       IF ( W1(NP)%LDO ) THEN

          

          CALL VHAMIL(WDES1, WDES1%GRID, SV(1,ISP), W1(NP)%CR(1), CWORK1(1,NP))
       ENDIF
    ENDDO

    IF (NONLR_S%LREAL) THEN
! non-local contribution in real-space
       CALL RACCMU_(NONLR_S, WDES1, W1, CDIJ, CQIJ, ISP, EVALUE, CWORK1)

       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN

             

             CALL KINHAMIL( WDES1, WDES1%GRID, CWORK1(1,NP), .FALSE., &
                  WDES1%DATAKE(1,1), EVALUE(NP), W1(NP)%CPTWFP(1), WRESULT%CPTWFP(1,NP))
          ENDIF
       ENDDO
    ELSE
! calculate the non local contribution in reciprocal space
       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN

             

             CALL VNLACC(NONL_S, W1(NP), CDIJ, CQIJ, ISP, EVALUE(NP),  WRESULT%CPTWFP(:,NP))
             CALL KINHAMIL( WDES1, WDES1%GRID, CWORK1(1,NP), .TRUE., &
                  WDES1%DATAKE(1,1), EVALUE(NP), W1(NP)%CPTWFP(1), WRESULT%CPTWFP(1,NP))
          ENDIF
       ENDDO
    ENDIF

# 1048


    

    RETURN
  END SUBROUTINE HAMILTMU


!************************* SUBROUTINE HAMILTMU_VEC *********************
!
!> this subroutine calculates the H acting onto a set of wavefuntions
!>
!> the  wavefunction must be given in reciprocal space C and real
!> space CR
!> CH contains the result
!>   H- EVALUE S |phi>
!> where S is the overlap matrix S
!> this version includes a vector potential A
!
!> @details @ref openmp :
!> under OpenMP the workarrays CWORK2 and CWORK3 acquire a second
!> dimension (size=3). This additional workspace improves the
!> performance of ::kinhamil_vec under OpenMP parallelization.
!
!***********************************************************************

  SUBROUTINE HAMILTMU_VEC( &
       &    WDES1, W1, NONLR_S, NONL_S, EVALUE, &
       &    CDIJ, CQIJ, SV, AVEC, ISP, WRESULT)
    USE mpimy
    USE mgrid
    USE nonl_high
    USE mopenmp_struct_def, ONLY : omp_nthreads_acc
    IMPLICIT NONE

    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1(:)
    TYPE (nonlr_struct)NONLR_S
    TYPE (nonl_struct) NONL_S    
    REAL(q)    EVALUE(:)               ! eigenvalues
    COMPLEX(q) CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    COMPLEX(q)      SV(:,:)
    COMPLEX(q)      AVEC(:,:)
    INTEGER    ISP
    TYPE(wavefuna)     WRESULT
    
! local variables
    COMPLEX(q) :: CWORK1(WDES1%GRID%MPLWV*WDES1%NRSPINORS,SIZE(W1))
# 1099

    COMPLEX(q) :: CWORK2(WDES1%NRPLWV     )
    COMPLEX(q) :: CWORK3(WDES1%GRID%MPLWV )

# 1105

    INTEGER NP

    

!$ACC ENTER DATA CREATE(CWORK1,CWORK2,CWORK3) IF(OFFLOAD_ON)

! calculate the local contribution (result in CWORK1)
!!!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(NP) NUM_THREADS(omp_nthreads_acc)
    DO NP=1,SIZE(W1)
       IF ( W1(NP)%LDO ) THEN

          

          CALL VHAMIL(WDES1, WDES1%GRID, SV(1,ISP), W1(NP)%CR(1), CWORK1(1,NP))
       ENDIF
    ENDDO
!!!$OMP END PARALLEL DO

    IF (NONLR_S%LREAL) THEN
! non-local contribution in real-space
       CALL RACCMU_(NONLR_S, WDES1, W1, CDIJ, CQIJ, ISP, EVALUE,CWORK1)
!!!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(NP) NUM_THREADS(omp_nthreads_acc)
       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN

             

             CALL KINHAMIL_VEC( WDES1, WDES1%GRID, CWORK1(1,NP), .FALSE., &
                  WDES1%DATAKE(1,1), WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), AVEC(1,1), &
# 1137

                  CWORK2(1 ),  CWORK3(1 ), &

                  EVALUE(NP), W1(NP)%CPTWFP(1), WRESULT%CPTWFP(1,NP))
          ENDIF
       ENDDO
!!!$OMP END PARALLEL DO
    ELSE
!!!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(NP) NUM_THREADS(omp_nthreads_acc)
       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN

             

! calculate the non local contribution in reciprocal space
             CALL VNLACC(NONL_S, W1(NP), CDIJ, CQIJ, ISP, EVALUE(NP),  WRESULT%CPTWFP(:,NP))
             CALL KINHAMIL_VEC( WDES1, WDES1%GRID, CWORK1(1,NP), .TRUE., &
                  WDES1%DATAKE(1,1), WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), AVEC(1,1), &
# 1157

                  CWORK2(1 ),  CWORK3(1 ), &

                  EVALUE(NP), W1(NP)%CPTWFP(1), WRESULT%CPTWFP(1,NP))
          ENDIF
       ENDDO
!!!$OMP END PARALLEL DO
    ENDIF

# 1170


    

    RETURN
  END SUBROUTINE HAMILTMU_VEC


!************************* SUBROUTINE HAMILTMU_TAU *********************
!
!> this subroutine calculates the H acting onto a set of wavefuntions
!> for a kinetic energy functional
!>
!> @details @ref openmp :
!> under OpenMP the workarrays CWORK2 and CWORK3 acquire a second
!> dimension (size=3 and 6, respectively). This additional workspace
!> improves the performance of ::kinhamil_tau under OpenMP
!> parallelization.
!
!***********************************************************************

  SUBROUTINE HAMILTMU_TAU( &
       &    WDES1, W1, NONLR_S, NONL_S, EVALUE, &
       &    CDIJ, CQIJ, SV, LATT_CUR, MU, ISP, WRESULT)
    USE mpimy
    USE mgrid
    USE lattice
    USE nonl_high
    USE mopenmp_struct_def, ONLY : omp_nthreads_acc
    IMPLICIT NONE

    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1(:)
    TYPE (nonlr_struct)NONLR_S
    TYPE (nonl_struct) NONL_S
    TYPE (latt)        LATT_CUR    
    REAL(q)    EVALUE(:)               ! eigenvalues
    COMPLEX(q) CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    COMPLEX(q)      SV(:,:)
    COMPLEX(q)      MU(:,:)
    INTEGER    ISP
    TYPE(wavefuna)     WRESULT
    
! local variables
    COMPLEX(q) :: CWORK1(WDES1%GRID%MPLWV*WDES1%NRSPINORS,SIZE(W1))
# 1218

    COMPLEX(q) :: CWORK2(WDES1%NRPLWV     )
    COMPLEX(q) :: CWORK3(WDES1%GRID%MPLWV )

# 1224


    INTEGER NP

    

!$ACC ENTER DATA CREATE(CWORK1,CWORK2,CWORK3) IF(OFFLOAD_ON)
!$ACC ENTER DATA COPYIN(LATT_CUR,LATT_CUR%B) IF(OFFLOAD_ON)

! calculate the local contribution (result in CWORK1)
!!!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(NP) NUM_THREADS(omp_nthreads_acc)
    DO NP=1,SIZE(W1)
       IF ( W1(NP)%LDO ) THEN

          

          CALL VHAMIL(WDES1, WDES1%GRID, SV(1,ISP), W1(NP)%CR(1), CWORK1(1,NP))
       ENDIF
    ENDDO
!!!$OMP END PARALLEL DO


    IF (NONLR_S%LREAL) THEN
! non-local contribution in real-space
       CALL RACCMU_(NONLR_S, WDES1, W1, CDIJ, CQIJ, ISP, EVALUE, CWORK1)
!!!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(NP) NUM_THREADS(omp_nthreads_acc)
       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN

             

             CALL KINHAMIL_TAU( WDES1, WDES1%GRID, CWORK1(1,NP), .FALSE., .TRUE., &
                  WDES1%DATAKE(1,1), WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), LATT_CUR, MU(1,ISP), &
# 1259

                  CWORK2(1 ),  CWORK3(1 ), &

                  EVALUE(NP), W1(NP)%CPTWFP(1), WRESULT%CPTWFP(1,NP))
          ENDIF
       ENDDO
!!!$OMP END PARALLEL DO
    ELSE
!!!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(NP) NUM_THREADS(omp_nthreads_acc)
       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN

             

! calculate the non local contribution in reciprocal space
             CALL VNLACC(NONL_S, W1(NP), CDIJ, CQIJ, ISP, EVALUE(NP),  WRESULT%CPTWFP(:,NP))
             CALL KINHAMIL_TAU( WDES1, WDES1%GRID, CWORK1(1,NP), .TRUE., .TRUE., &
                  WDES1%DATAKE(1,1), WDES1%IGX(1), WDES1%IGY(1), WDES1%IGZ(1), WDES1%VKPT(1), LATT_CUR, MU(1,ISP), &
# 1279

                  CWORK2(1 ),  CWORK3(1 ), &

                  EVALUE(NP), W1(NP)%CPTWFP(1), WRESULT%CPTWFP(1,NP))
          ENDIF
       ENDDO
!!!$OMP END PARALLEL DO
    ENDIF

# 1292


    

    RETURN
  END SUBROUTINE HAMILTMU_TAU


!************************* SUBROUTINE HAMILTMU_C ***********************
!
!> identical to previous subroutine but with a complex eigenvalue EVALUE
!
!***********************************************************************

  SUBROUTINE HAMILTMU_C( &
       &    WDES1, W1, NONLR_S, NONL_S, EVALUE, &
       &    CDIJ, CQIJ, SV, ISP, WRESULT)
    USE mpimy
    USE mgrid
    USE nonl_high
    IMPLICIT NONE

    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1(:)
    TYPE (nonlr_struct)NONLR_S
    TYPE (nonl_struct) NONL_S    
    COMPLEX(q)    EVALUE(:)               ! eigenvalues
    COMPLEX(q) CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    COMPLEX(q)      SV(:,:)
    INTEGER    ISP
    TYPE(wavefuna)     WRESULT
    
! local variables
    COMPLEX(q) :: CWORK1(WDES1%GRID%MPLWV*WDES1%NRSPINORS,SIZE(W1))
    INTEGER NP
# 1329


    

!$ACC ENTER DATA CREATE(CWORK1) IF(OFFLOAD_ON)

! calculate the local contribution (result in CWORK1)
    DO NP=1,SIZE(W1)
       IF ( W1(NP)%LDO ) THEN

          

          CALL VHAMIL(WDES1, WDES1%GRID, SV(1,ISP), W1(NP)%CR(1), CWORK1(1,NP))
       ENDIF
    ENDDO

    IF (NONLR_S%LREAL) THEN
! non-local contribution in real-space
       CALL RACCMU_C_(NONLR_S, WDES1, W1, CDIJ, CQIJ, ISP, EVALUE,CWORK1)

       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN

             

             CALL KINHAMIL_C( WDES1, WDES1%GRID, CWORK1(1,NP), .FALSE., &
                  WDES1%DATAKE(1,1), EVALUE(NP), W1(NP)%CPTWFP(1), WRESULT%CPTWFP(1,NP))
          ENDIF
       ENDDO
    ELSE
! calculate the non local contribution in reciprocal space
       DO NP=1,SIZE(W1)
          IF ( W1(NP)%LDO ) THEN

             

             CALL VNLACC_C(NONL_S, W1(NP), CDIJ, CQIJ, ISP, EVALUE(NP),  WRESULT%CPTWFP(:,NP))
             CALL KINHAMIL_C( WDES1, WDES1%GRID, CWORK1(1,NP), .TRUE., &
                  WDES1%DATAKE(1,1), EVALUE(NP), W1(NP)%CPTWFP(1), WRESULT%CPTWFP(1,NP))
          ENDIF
       ENDDO
    ENDIF

# 1376


    

    RETURN
  END SUBROUTINE HAMILTMU_C

!************************* SUBROUTINE HAMILT_LOCAL *********************
!
!> this subroutine calculates the local and kinetic energy part of H
!> acting onto a wavefuntion
!>
!> the  wavefunction must be given in reciprocal space C and real
!> space CR
!> LADD=.TRUE.
!> ~~~
!>  CH = CH + SV* W1 + T * W1
!> ~~~
!> LADD=.FALSE.
!>  CH = SV* W1 + T * W1
!> LKIN determines whether the kinetic energy contribution is added or not
!> the default is add the kinetic energy contribution
!
!***********************************************************************

  SUBROUTINE HAMILT_LOCAL(W1, SV, ISP, CH, LADD, LKIN)
    USE mgrid
    IMPLICIT NONE

    TYPE (wavefun1)    W1
    COMPLEX(q)   SV(:,:)
    INTEGER :: ISP
    COMPLEX(q) CH(:)
    LOGICAL LADD
    LOGICAL, OPTIONAL :: LKIN
! local variables
    COMPLEX(q) :: CWORK1(W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS)
# 1415

     
!$ACC DATA CREATE(CWORK1) 
    CALL VHAMIL(W1%WDES1,W1%WDES1%GRID,SV(1,ISP),W1%CR(1),CWORK1(1)) 
    IF (.NOT. PRESENT(LKIN)) THEN
       CALL KINHAMIL( W1%WDES1, W1%WDES1%GRID, CWORK1(1), LADD, &
            W1%WDES1%DATAKE(1,1), 0.0_q, W1%CPTWFP(1), CH(1))
    ELSE
       IF (LKIN) THEN
          CALL KINHAMIL( W1%WDES1, W1%WDES1%GRID, CWORK1(1), LADD, &
               W1%WDES1%DATAKE(1,1), 0.0_q, W1%CPTWFP(1), CH(1))
       ELSE
          CALL FFTHAMIL( W1%WDES1, W1%WDES1%GRID, CWORK1(1), LADD, &
                                     0.0_q, W1%CPTWFP(1), CH(1))
       ENDIF
    ENDIF
!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
!$ACC END DATA
    

  END SUBROUTINE HAMILT_LOCAL


  SUBROUTINE HAMILT_LOCAL_TAU(W1, SV, LATT_CUR, MU, ISP, CH, LADD, LKIN)
    USE mgrid
    USE lattice
    IMPLICIT NONE

    TYPE (wavefun1)    W1
    TYPE (latt)        LATT_CUR
    COMPLEX(q)   SV(:,:)
    COMPLEX(q)   MU(:,:)
    INTEGER ISP
    COMPLEX(q) CH(:)
    LOGICAL LADD
    LOGICAL, OPTIONAL :: LKIN
! local variables
    COMPLEX(q) :: CVR(W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS)
    COMPLEX(q) :: CWORK1(W1%WDES1%NRPLWV     )
    COMPLEX(q) :: CWORK2(W1%WDES1%GRID%MPLWV )
    LOGICAL :: LKINETIC=.TRUE.
# 1458

    

    IF (PRESENT(LKIN)) LKINETIC=LKIN
!$ACC ENTER DATA CREATE(CVR,CWORK1,CWORK2) 
!$ACC ENTER DATA COPYIN(LATT_CUR,LATT_CUR%B) 
    CALL VHAMIL(W1%WDES1,W1%WDES1%GRID,SV(1,ISP),W1%CR(1),CVR(1))

    CALL KINHAMIL_TAU( W1%WDES1, W1%WDES1%GRID, CVR(1), LADD, LKINETIC, W1%WDES1%DATAKE(1,1), &
   &   W1%WDES1%IGX(1), W1%WDES1%IGY(1), W1%WDES1%IGZ(1), W1%WDES1%VKPT(1), LATT_CUR, MU(1,ISP), CWORK1, CWORK2, &
   &   0.0_q, W1%CPTWFP(1), CH(1))
!$ACC EXIT DATA DELETE(CVR,CWORK1,CWORK2,LATT_CUR%B,LATT_CUR) 
    

  END SUBROUTINE HAMILT_LOCAL_TAU


!************************* SUBROUTINE ECCP_NL_ALL   ********************
!
!> this subroutine calculates the non local contribution to the
!> expectation value of sum_ij <c2|p_j> D_ij - e Q_ij <p_i|c1>
!> where c and cp are two wavefunctions
!
!***********************************************************************

  SUBROUTINE ECCP_NL_ALL(WDES1,W1,W2,CDIJ,CQIJ,EVALUE,CE)

!! USE moffload_struct_def

    USE mpimy
    USE mgrid
    IMPLICIT NONE

    TYPE (wavefun1) :: W1,W2
    TYPE (wavedes1) :: WDES1
    TYPE (grid_3d)  :: GRID

    COMPLEX(q) CDIJ(:,:,:,:)
    COMPLEX(q) CQIJ(:,:,:,:)
    REAL(q) EVALUE
! return
    COMPLEX(q) CE
! local
    COMPLEX(q) CTMP
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NI, NIS, NT, LMMAXC

!$ACC ROUTINE(ECCP_NL_) VECTOR

!$ACC KERNELS PRESENT(CE) 
    CE=0
!$ACC END KERNELS

!$ACC PARALLEL LOOP REDUCTION(+:CE) PRESENT(CE,CDIJ,CQIJ,W1%CPROJ,W2%CPROJ) &
!$ACC COLLAPSE(3) PRIVATE(CTMP,NPRO,NPRO_,NT,LMMAXC) 
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
       DO ISPINOR_=0,WDES1%NRSPINORS-1

          DO NI=1,WDES1%NIONS
             NT=WDES1%ITYP(NI)
             LMMAXC=WDES1%LMMAX(NT)
             NPRO =ISPINOR *(WDES1%NPRO/2)+WDES1%LMBASE(NI)
             NPRO_=ISPINOR_*(WDES1%NPRO/2)+WDES1%LMBASE(NI)
             IF (LMMAXC==0) CYCLE
             CTMP=0
             CALL ECCP_NL_(SIZE(CDIJ,1),LMMAXC, &
                  CDIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),CQIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),EVALUE, &
                  W1%CPROJ(NPRO_+1),W2%CPROJ(NPRO+1),CTMP)
             CE=CE+CTMP
          ENDDO
       ENDDO
    ENDDO spinor

    CALL M_sum_z(WDES1%COMM_INB, CE, 1)

    RETURN
  END SUBROUTINE ECCP_NL_ALL


!************************* SUBROUTINE ECCP_NL_ALL_ION ******************
!
!> this subroutine calculates the non local contribution to the
!> expectation value of sum_ij <c2|p_j> D_ij - e Q_ij <p_i|c1>
!> where c and cp are two wavefunctions
!
!***********************************************************************

  SUBROUTINE ECCP_NL_ALL_ION(WDES1,W1,W2,CDIJ,CQIJ,EVALUE,CE)
    USE mpimy
    USE mgrid
    IMPLICIT NONE

    TYPE (wavefun1) :: W1,W2
    TYPE (wavedes1) :: WDES1
    TYPE (grid_3d)  :: GRID

    COMPLEX(q) CDIJ(:,:,:,:)
    COMPLEX(q) CQIJ(:,:,:,:)
    REAL(q) EVALUE
! return
    COMPLEX(q) CE(*)
! local
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NI, NIS, NT, LMMAXC

    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
       DO ISPINOR_=0,WDES1%NRSPINORS-1

          NPRO =ISPINOR *(WDES1%NPRO/2)
          NPRO_=ISPINOR_*(WDES1%NPRO/2)

          NIS =1
          DO NT=1,WDES1%NTYP
             LMMAXC=WDES1%LMMAX(NT)
             IF (LMMAXC==0) GOTO 310
             DO NI=NIS,WDES1%NITYP(NT)+NIS-1
                CALL ECCP_NL_(SIZE(CDIJ,1),LMMAXC, &
                     CDIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),CQIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),EVALUE, &
                     W1%CPROJ(NPRO_+1),W2%CPROJ(NPRO+1),CE(NI))
                NPRO = LMMAXC+NPRO
                NPRO_= LMMAXC+NPRO_
             ENDDO
310          NIS = NIS+WDES1%NITYP(NT)
          ENDDO
       ENDDO
    ENDDO spinor

    RETURN
  END SUBROUTINE ECCP_NL_ALL_ION

END MODULE hamil

!************************* SUBROUTINE VHAMIL  **************************
!
!> low level F77 to calculate the product of the local potential
!> times a wavefunction stored in CR and returns the result in CVR
!> ~~~
!>   CVR= SV * CR
!> ~~~
!> this low level routine should be called only from within
!> hamil.F
!
!***********************************************************************

   SUBROUTINE VHAMIL(WDES1,GRID,SV,CR,CVR)

!! USE moffload

     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (grid_3d)     GRID
     TYPE (wavedes1)    WDES1

     COMPLEX(q)   SV(GRID%MPLWV,WDES1%NRSPINORS*WDES1%NRSPINORS) ! local potential
     COMPLEX(q) :: CR(GRID%MPLWV*WDES1%NRSPINORS),CVR(GRID%MPLWV*WDES1%NRSPINORS)
! local variables
     REAL(q) RINPLW
     INTEGER ISPINOR, ISPINOR_, M, MM, MM_
# 1619


     

     RINPLW=1._q/GRID%NPLWV

     IF (WDES1%NRSPINORS==1) THEN
!$ACC PARALLEL LOOP GANG VECTOR PRESENT(GRID,CVR,SV,CR) &
!$ACC  DEFAULT(PRESENT)
!DIR$ IVDEP
!OCL NOVREL
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
 !$OMP PRIVATE(M) SHARED(GRID, CVR, SV, CR, RINPLW)
        DO M=1,GRID%RL%NP
           CVR(M)= SV(M,1) *CR(M)*RINPLW
        ENDDO
 !$OMP END PARALLEL DO
     ELSE
!!!        CVR(1:GRID%MPLWV*2)=0
!!!!$OMP PARALLEL &
!!!!$OMP SHARED(GRID,CVR,SV,CR,RINPLW) &
!!!!$OMP PRIVATE(ISPINOR,ISPINOR_,M,MM,MM_)
!!!        DO ISPINOR =0,1
!!!           DO ISPINOR_=0,1
!!!!DIR$ IVDEP
!!!!OCL NOVREL
!!!!$OMP DO SCHEDULE(STATIC)
!!!              DO M=1,GRID%RL%NP
!!!                 MM =M+ISPINOR *GRID%MPLWV
!!!                 MM_=M+ISPINOR_*GRID%MPLWV
!!!                 CVR(MM)= CVR(MM)+ SV(M,1+ISPINOR_+2*ISPINOR) *CR(MM_)*RINPLW
!!!              ENDDO
!!!!$OMP END DO
!!!           ENDDO
!!!        ENDDO
!!!!$OMP END PARALLEL

!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(GRID,CVR,CR) PRIVATE(MM) 
 !$OMP PARALLEL SHARED(GRID,CVR,SV,CR,RINPLW) PRIVATE(ISPINOR,M,MM)
        DO ISPINOR =0,1
!DIR$ IVDEP
!OCL NOVREL
 !$OMP DO SCHEDULE(STATIC)
           DO M=1,GRID%RL%NP
              MM =M+ISPINOR *GRID%MPLWV
              CVR(MM)=SV(M,1+2*ISPINOR) *CR(M)*RINPLW
           ENDDO
 !$OMP END DO
!!   ENDDO

!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(GRID,CVR,SV,CR) PRIVATE(MM) 
!!   DO ISPINOR =0,1
!DIR$ IVDEP
!OCL NOVREL
 !$OMP DO SCHEDULE(STATIC)
           DO M=1,GRID%RL%NP
              MM =M+ISPINOR *GRID%MPLWV
              CVR(MM)= CVR(MM)+ SV(M,2+2*ISPINOR) *CR(M+GRID%MPLWV)*RINPLW
           ENDDO
 !$OMP END DO
        ENDDO
 !$OMP END PARALLEL
     ENDIF

     

   END SUBROUTINE VHAMIL

   SUBROUTINE VHAMIL_TRACE(WDES1, GRID, SV, CR, CVR, WEIGHT)

!! USE moffload

     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (grid_3d)  :: GRID
     TYPE (wavedes1) :: WDES1

     COMPLEX(q) :: CR( GRID%MPLWV*WDES1%NRSPINORS), &
                   CVR(GRID%MPLWV*WDES1%NRSPINORS)
     COMPLEX(q)       :: SV(GRID%MPLWV)                  ! local potential

     REAL(q) :: WEIGHT

! local variables
     INTEGER :: M, MM

     

 !$OMP PARALLEL SHARED(GRID,CVR,SV,CR,WEIGHT) PRIVATE(M,MM)
!DIR$ IVDEP
!OCL NOVREL
 !$OMP DO SCHEDULE(STATIC)
!$ACC PARALLEL LOOP PRESENT(GRID,CVR,SV,CR) 
     DO M=1,GRID%RL%NP
        CVR(M)= CVR(M)+SV(M) *CR(M)*WEIGHT
     ENDDO
 !$OMP END DO
     IF (WDES1%NRSPINORS==2) THEN
!DIR$ IVDEP
!OCL NOVREL
 !$OMP DO SCHEDULE(STATIC)
!$ACC PARALLEL LOOP PRESENT(GRID,CVR,SV,CR) PRIVATE(MM) 
        DO M=1,GRID%RL%NP
           MM =M+GRID%MPLWV
           CVR(MM)= CVR(MM)+ SV(M) *CR(MM)*WEIGHT
        ENDDO
 !$OMP END DO
     ENDIF
 !$OMP END PARALLEL

     

   END SUBROUTINE VHAMIL_TRACE

   SUBROUTINE VHAMIL_TRACE_MU(WDES1, GRID, SV, LD, CR, WXI, WEIGHT, NSIM)

!! USE moffload

     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (grid_3d)  :: GRID
     TYPE (wavedes1) :: WDES1
     TYPE (wavefun1) :: WXI(NSIM)

     COMPLEX(q) :: CR(GRID%MPLWV*WDES1%NRSPINORS)
     COMPLEX(q)       :: SV(LD,NSIM)                     ! local potential

     REAL(q) :: WEIGHT
     INTEGER :: LD,NSIM

! local variables
     INTEGER ISPINOR, M, MM, N, IBLOCK
     INTEGER, PARAMETER :: BLOCKSIZE=128

     

!DIR$ IVDEP
!OCL NOVREL
!$ACC PARALLEL LOOP GANG COLLAPSE(2) PRESENT(WXI,SV,CR,GRID) 
     DO N=1,NSIM
        DO IBLOCK=1,WDES1%GRID%RL%NP,BLOCKSIZE

!$ACC CACHE(CR(IBLOCK:MIN(IBLOCK+BLOCKSIZE-1,WDES1%GRID%RL%NP)))

!$ACC LOOP VECTOR
           DO M=IBLOCK,MIN(IBLOCK+BLOCKSIZE-1,WDES1%GRID%RL%NP)
              WXI(N)%CR(M)= WXI(N)%CR(M)+SV(M,N) *CR(M)*WEIGHT
           ENDDO
        ENDDO
     ENDDO
     IF (WDES1%NRSPINORS==2) THEN
!DIR$ IVDEP
!OCL NOVREL
!$ACC PARALLEL LOOP GANG COLLAPSE(2) PRESENT(WXI,SV,CR,GRID) 
        DO N=1,NSIM
           DO IBLOCK=1,WDES1%GRID%RL%NP,BLOCKSIZE

!$ACC CACHE(CR(GRID%MPLWV+IBLOCK:GRID%MPLWV+MIN(IBLOCK+BLOCKSIZE-1,WDES1%GRID%RL%NP)))

!$ACC LOOP VECTOR PRIVATE(MM)
              DO M=IBLOCK,MIN(IBLOCK+BLOCKSIZE-1,WDES1%GRID%RL%NP)
                 MM =M+GRID%MPLWV
                 WXI(N)%CR(MM)= WXI(N)%CR(MM)+SV(M,N) *CR(MM)*WEIGHT
              ENDDO
           ENDDO
        ENDDO
     ENDIF

     

   END SUBROUTINE VHAMIL_TRACE_MU
# 1869

!************************* SUBROUTINE PW_CHARGE ************************
!
!> low level F77 to calculate charge density from two plane wave states
!> and add the weighted spinor like density to an array
!> this low level routine should be called only from within
!> hamil.F
!> ~~~
!>  CHARGE = CR1*CONJG(CR2) * WEIGHT + CHARGE
!> ~~~
!
!***********************************************************************

!
! this version yields an COMPLEX(q) array spinor density
!

   SUBROUTINE PW_CHARGE(WDES1,  CHARGE, NDIM, CR1, CR2, WEIGHT)

!! USE moffload

     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (grid_3d)     GRID
     TYPE (wavedes1)    WDES1
     INTEGER NDIM
     COMPLEX(q)   CHARGE(NDIM,WDES1%NRSPINORS*WDES1%NRSPINORS)
     COMPLEX(q) :: CR1(WDES1%GRID%MPLWV*WDES1%NRSPINORS),CR2(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
     REAL(q) :: WEIGHT
! local variables
     INTEGER ISPINOR,ISPINOR_,M,MM,MM_

     

     IF (WDES1%NRSPINORS==1) THEN
!$ACC PARALLEL LOOP PRESENT(CHARGE,CR1,CR2) 
!DIR$ IVDEP
!OCL NOVREL
 !$OMP PARALLEL DO PRIVATE(M) SHARED(WDES1,CHARGE,CR1,CR2,WEIGHT)
        DO M=1,WDES1%GRID%RL%NP
           CHARGE(M,1)=CHARGE(M,1)+ & 
                CONJG(CR2(M))*CR1(M)*WEIGHT
        ENDDO
 !$OMP END PARALLEL DO
     ELSE
!$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(CHARGE,CR1,CR2) 
 !$OMP PARALLEL &
 !$OMP SHARED(WDES1,CHARGE,CR1,CR2,WEIGHT) &
 !$OMP PRIVATE(ISPINOR,ISPINOR_,M,MM,MM_)
        spinor: DO ISPINOR=0,1
           DO ISPINOR_=0,1
!DIR$ IVDEP
!OCL NOVREL
 !$OMP DO
              DO M=1,WDES1%GRID%RL%NP
                 MM =M+ISPINOR *WDES1%GRID%MPLWV
                 MM_=M+ISPINOR_*WDES1%GRID%MPLWV
                 CHARGE(M,1+ISPINOR_+2*ISPINOR)=CHARGE(M,1+ISPINOR_+2*ISPINOR)+ & 
                      CONJG(CR2(MM_))*CR1(MM)*WEIGHT
              ENDDO
 !$OMP END DO NOWAIT
           ENDDO
        ENDDO spinor
 !$OMP END PARALLEL
     ENDIF

     

   END SUBROUTINE PW_CHARGE

!
!> this version yields a COMPLEX(q) charge array, spinor density
!> (used by exchange routines)
!
   SUBROUTINE PW_CHARGE_CMPLX(WDES1,  CHARGE, NDIM, CR1, CR2 )
     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (grid_3d)     GRID
     TYPE (wavedes1)    WDES1
     INTEGER NDIM
     COMPLEX(q)   CHARGE(NDIM,WDES1%NRSPINORS*WDES1%NRSPINORS)
     COMPLEX(q) :: CR1(WDES1%GRID%MPLWV*WDES1%NRSPINORS),CR2(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
     REAL(q) :: WEIGHT
! local variables
     INTEGER ISPINOR,ISPINOR_,M,MM,MM_

     IF (WDES1%NRSPINORS==1) THEN
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%GRID%RL%NP
           CHARGE(M,1)=CONJG(CR2(M))*CR1(M)
        ENDDO
     ELSE
        CHARGE=0
        spinor: DO ISPINOR=0,1
           DO ISPINOR_=0,1
!DIR$ IVDEP
!OCL NOVREL
              DO M=1,WDES1%GRID%RL%NP
                 MM =M+ISPINOR *WDES1%GRID%MPLWV
                 MM_=M+ISPINOR_*WDES1%GRID%MPLWV
                 CHARGE(M,1+ISPINOR_+2*ISPINOR)=CHARGE(M,1+ISPINOR_+2*ISPINOR)+ & 
                      CONJG(CR2(MM_))*CR1(MM)
              ENDDO
           ENDDO
        ENDDO spinor
     ENDIF

   END SUBROUTINE PW_CHARGE_CMPLX

!
!> this version traces of spinors to yield only the density
!> results array is COMPLEX(q) (used by exchange routines)
!
   SUBROUTINE PW_CHARGE_TRACE(WDES1, CHARGE, CR1, CR2)

!! USE moffload_struct_def

     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (grid_3d)  :: GRID
     TYPE (wavedes1) :: WDES1

     COMPLEX(q) :: CR1(WDES1%GRID%MPLWV*WDES1%NRSPINORS),CR2(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
     COMPLEX(q)       :: CHARGE(WDES1%GRID%MPLWV)

! local variables
     REAL(q) RINPLW
     INTEGER M,MM

     

 !$OMP PARALLEL SHARED(WDES1,CHARGE,CR1,CR2) PRIVATE(M,MM)
!DIR$ IVDEP
!OCL NOVREL
 !$OMP DO SCHEDULE(STATIC)
!$ACC PARALLEL LOOP PRESENT(WDES1,CHARGE,CR2,CR1) 
     DO M=1,WDES1%GRID%RL%NP
        CHARGE(M)=CONJG(CR2(M))*CR1(M)
     ENDDO
 !$OMP END DO
     IF (WDES1%NRSPINORS==2) THEN
!DIR$ IVDEP
!OCL NOVREL
 !$OMP DO SCHEDULE(STATIC)
!$ACC PARALLEL LOOP PRESENT(WDES1,CHARGE,CR2,CR1) PRIVATE(MM) 
        DO M=1,WDES1%GRID%RL%NP
           MM =M+WDES1%GRID%MPLWV
           CHARGE(M)=CHARGE(M)+CONJG(CR2(MM))*CR1(MM)
        ENDDO
 !$OMP END DO
     ENDIF
 !$OMP END PARALLEL

     

   END SUBROUTINE PW_CHARGE_TRACE

   SUBROUTINE PW_CHARGE_TRACE_MU(WDES1, CHARGE, LD, W1, CR2, NSIM)

!! USE moffload_struct_def

     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (wavedes1) :: WDES1
     TYPE (wavefun1) :: W1(NSIM)

     COMPLEX(q) :: CR2(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
     COMPLEX(q)       :: CHARGE(LD,NSIM)

     INTEGER :: LD,NSIM

     

# 2061

     CALL PW_CHARGE_TRACE_MU_(WDES1, CHARGE(1,1), LD, W1, CR2(1), NSIM)

     

     RETURN
   END SUBROUTINE PW_CHARGE_TRACE_MU

   SUBROUTINE PW_CHARGE_TRACE_MU_(WDES1, CHARGE, LD, W1, CR2, NSIM)
     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (wavedes1) :: WDES1
     TYPE (wavefun1) :: W1(NSIM)

     COMPLEX(q) :: CR2(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
     COMPLEX(q)       :: CHARGE(LD,NSIM)

     INTEGER :: LD,NSIM

! local variables
     INTEGER :: N,M,MM

     

!DIR$ IVDEP
!OCL NOVREL
     DO N=1,NSIM
        DO M=1,WDES1%GRID%RL%NP
           CHARGE(M,N)=CONJG(CR2(M))*W1(N)%CR(M)
        ENDDO
     ENDDO
     IF (WDES1%NRSPINORS==2) THEN
!DIR$ IVDEP
!OCL NOVREL
        DO N=1,NSIM
           DO M=1,WDES1%GRID%RL%NP
              MM =M+WDES1%GRID%MPLWV
              CHARGE(M,N)=CHARGE(M,N)+CONJG(CR2(MM))*W1(N)%CR(MM)
           ENDDO
        ENDDO
     ENDIF

     

   END SUBROUTINE PW_CHARGE_TRACE_MU_

   SUBROUTINE PW_CHARGE_TRACE_NO_CONJG_MU(WDES1, CHARGE, LD, W1, CR2, NSIM)

!! USE moffload

     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

    TYPE (wavedes1) :: WDES1
    TYPE (wavefun1) :: W1(NSIM)

    COMPLEX(q) :: CR2(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
    COMPLEX(q)       :: CHARGE(LD,NSIM)

    INTEGER :: LD,NSIM

! local variables
    INTEGER :: N,M,MM,IBLOCK
    INTEGER, PARAMETER :: BLOCKSIZE=128

    

!DIR$ IVDEP
!OCL NOVREL
!$ACC PARALLEL LOOP GANG COLLAPSE(2) PRESENT(CHARGE,CR2,W1) 
    DO N=1,NSIM
       DO IBLOCK=1,WDES1%GRID%RL%NP,BLOCKSIZE

!$ACC CACHE(CR2(IBLOCK:MIN(IBLOCK+BLOCKSIZE-1,WDES1%GRID%RL%NP)))

!$ACC LOOP VECTOR
          DO M=IBLOCK,MIN(IBLOCK+BLOCKSIZE-1,WDES1%GRID%RL%NP)
             CHARGE(M,N)=CR2(M)*W1(N)%CR(M)
          ENDDO
       ENDDO
    ENDDO
    IF (WDES1%NRSPINORS==2) THEN
!DIR$ IVDEP
!OCL NOVREL
!$ACC PARALLEL LOOP GANG COLLAPSE(2) PRESENT(CHARGE,CR2,W1) 
       DO N=1,NSIM
          DO IBLOCK=1,WDES1%GRID%RL%NP,BLOCKSIZE

!$ACC CACHE(CR2(WDES1%GRID%MPLWV+IBLOCK:WDES1%GRID%MPLWV+MIN(IBLOCK+BLOCKSIZE-1,WDES1%GRID%RL%NP)))

!$ACC LOOP VECTOR PRIVATE(MM)
             DO M=IBLOCK,MIN(IBLOCK+BLOCKSIZE-1,WDES1%GRID%RL%NP)
                MM =M+WDES1%GRID%MPLWV
                CHARGE(M,N)=CHARGE(M,N)+CR2(MM)*W1(N)%CR(MM)
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    

  END SUBROUTINE PW_CHARGE_TRACE_NO_CONJG_MU

# 2317


   SUBROUTINE PW_CHARGE_TRACE_MU2(WDES1, CHARGE, LDC, CR1, LD1, CR2, LD2, NSIM)

!! USE moffload_struct_def

     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (wavedes1) :: WDES1

     COMPLEX(q)       :: CHARGE(LDC,NSIM)
     COMPLEX(q) :: CR1(LD1,NSIM),CR2(LD1,NSIM)

     INTEGER :: LDC,LD1,LD2,NSIM

! local variables
     INTEGER :: N

     

# 2346

     DO N=1,NSIM
        CALL PW_CHARGE_TRACE(WDES1, CHARGE(1,N), CR1(1,N), CR2(1,N))
     ENDDO

     

     RETURN
   END SUBROUTINE PW_CHARGE_TRACE_MU2

# 2397

!
!> this version traces of spinors to yield only the density
!> result array as well as orbital input is COMPLEX(q)

   SUBROUTINE PW_CHARGE_TRACE_GDEF(WDES1,  CHARGE, CR1, CR2)
     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (grid_3d)     GRID
     TYPE (wavedes1)    WDES1
     COMPLEX(q) :: CHARGE(WDES1%GRID%MPLWV)
     COMPLEX(q) :: CR1(WDES1%GRID%MPLWV*WDES1%NRSPINORS),CR2(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
! local variables
     REAL(q) RINPLW
     INTEGER M,MM

!DIR$ IVDEP
!OCL NOVREL
      DO M=1,WDES1%GRID%RL%NP
         CHARGE(M)=CONJG(CR2(M))*CR1(M)
      ENDDO
      IF (WDES1%NRSPINORS==2) THEN
!DIR$ IVDEP
!OCL NOVREL
         DO M=1,WDES1%GRID%RL%NP
            MM =M+WDES1%GRID%MPLWV
            CHARGE(M)=CHARGE(M)+CONJG(CR2(MM))*CR1(MM)
         ENDDO
      ENDIF

   END SUBROUTINE PW_CHARGE_TRACE_GDEF

!
!> this version traces of spinors to yield only the density
!> result array with a REAL(q) layout (used by chi_super.F)
!
   SUBROUTINE PW_CHARGE_TRACE_REAL(WDES1,  CHARGE, CR1, CR2)

!! USE moffload_struct_def

     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (grid_3d)  :: GRID
     TYPE (wavedes1) :: WDES1
     COMPLEX(q) :: CR1(WDES1%GRID%MPLWV*WDES1%NRSPINORS),CR2(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
     REAL(q)    :: CHARGE(WDES1%GRID%MPLWV)

! local variables
     REAL(q) RINPLW
     INTEGER M,MM

     

 !$OMP PARALLEL SHARED(WDES1,CHARGE,CR1,CR2) !$OMP PRIVATE(M,MM)
!DIR$ IVDEP
!OCL NOVREL
 !$OMP DO SCHEDULE(STATIC)
!$ACC PARALLEL LOOP PRESENT(WDES1,CHARGE,CR2,CR1) 
     DO M=1,WDES1%GRID%RL%NP
        CHARGE(M)=CONJG(CR2(M))*CR1(M)
     ENDDO
 !$OMP END DO
     IF (WDES1%NRSPINORS==2) THEN
!DIR$ IVDEP
!OCL NOVREL
 !$OMP DO SCHEDULE(STATIC)
!$ACC PARALLEL LOOP PRESENT(WDES1,CHARGE,CR2,CR1) PRIVATE(MM) 
        DO M=1,WDES1%GRID%RL%NP
           MM =M+WDES1%GRID%MPLWV
           CHARGE(M)=CHARGE(M)+CONJG(CR2(MM))*CR1(MM)
        ENDDO
 !$OMP END DO
     ENDIF
 !$OMP END PARALLEL

     

   END SUBROUTINE PW_CHARGE_TRACE_REAL


   SUBROUTINE PW_CHARGE_TRACE_NO_CONJG(WDES1,  CHARGE, CR1, CR2)
     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (grid_3d)     GRID
     TYPE (wavedes1)    WDES1
     COMPLEX(q)   CHARGE(WDES1%GRID%MPLWV)
     COMPLEX(q) :: CR1(WDES1%GRID%MPLWV*WDES1%NRSPINORS),CR2(WDES1%GRID%MPLWV*WDES1%NRSPINORS)
! local variables
     REAL(q) RINPLW
     INTEGER M,MM

!DIR$ IVDEP
!OCL NOVREL
      DO M=1,WDES1%GRID%RL%NP
         CHARGE(M)=CR2(M)*CR1(M)
      ENDDO
      IF (WDES1%NRSPINORS==2) THEN
!DIR$ IVDEP
!OCL NOVREL
         DO M=1,WDES1%GRID%RL%NP
            MM =M+WDES1%GRID%MPLWV
            CHARGE(M)=CHARGE(M)+CR2(MM)*CR1(MM)
         ENDDO
      ENDIF

   END SUBROUTINE PW_CHARGE_TRACE_NO_CONJG


!************************* SUBROUTINE KINHAMIL  ************************
!
!> low level F77 to calculate
!> the kinetic energy part of a Hamiltonian
!> possibly substract the wavefunction times an eigenvalue
!> and adds a component after FFT to real space
!> this low level routine should be called only from within
!> hamil.F
!> ~~~
!> LADD=.TRUE.
!> ~~~
!> ~~~
!>  CH= CH + FFT(CVR) + DATAKE*CPTWFP - EVALUE*CPTWFP
!> ~~~
!> ~~~
!> LADD=.FALSE.
!> ~~~
!> ~~~
!>  CH= FFT(CVR) + DATAKE*CPTWFP - EVALUE*CPTWFP
!> ~~~
!>
!> FFTHAMIL exclude the kinetic energy contribution
!
!***********************************************************************

   SUBROUTINE KINHAMIL( WDES1, GRID, CVR,  LADD, DATAKE, EVALUE, CPTWFP, CH)

!! USE moffload

     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (wavedes1)    WDES1
     TYPE (grid_3d)     GRID
     COMPLEX(q) :: CVR(GRID%MPLWV*WDES1%NRSPINORS) !< usually potential times wavefunction
     LOGICAL    :: LADD                            !< if .TRUE. add results to CH
     REAL(q)    :: DATAKE(WDES1%NGDIM,WDES1%NRSPINORS)
     REAL(q)    :: EVALUE                          !< subtract EVALUE*wavefunction
     COMPLEX(q) :: CPTWFP(WDES1%NRPLWV)                !< wavefunction
     COMPLEX(q) :: CH(WDES1%NRPLWV)                !< result
! local
     INTEGER ISPINOR, M, MM
# 2560


     

     DO ISPINOR=0,WDES1%NRSPINORS-1
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CH(1+ISPINOR*WDES1%NGVECTOR),GRID,LADD)
!$ACC PARALLEL LOOP PRESENT(CH,CPTWFP,DATAKE,WDES1) PRIVATE(MM) 
!DIR$ IVDEP
!OCL NOVREL
 !$OMP PARALLEL DO SHARED(WDES1,CH,ISPINOR,DATAKE,EVALUE) PRIVATE (M,MM)
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)+CPTWFP(MM)*(DATAKE(M,ISPINOR+1)-EVALUE)
        ENDDO
 !$OMP END PARALLEL DO
     ENDDO

     

   END SUBROUTINE KINHAMIL


   SUBROUTINE FFTHAMIL( WDES1, GRID, CVR, LADD, EVALUE, CPTWFP, CH)

!! USE moffload

     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (wavedes1)    WDES1
     TYPE (grid_3d)     GRID
     COMPLEX(q) :: CVR(GRID%MPLWV*WDES1%NRSPINORS) !< usually potential times wavefunction
     LOGICAL    :: LADD                            !< if .TRUE. add results to CH
     REAL(q)    :: EVALUE                          !< subtract EVALUE*wavefunction
     COMPLEX(q) :: CPTWFP(WDES1%NRPLWV)                !< wavefunction
     COMPLEX(q) :: CH(WDES1%NRPLWV)                !< result
! local
     INTEGER ISPINOR, M, MM

     

     DO ISPINOR=0,WDES1%NRSPINORS-1
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CH(1+ISPINOR*WDES1%NGVECTOR),GRID,LADD)
!$ACC PARALLEL LOOP PRESENT(CH,CPTWFP) PRIVATE(MM) 
!DIR$ IVDEP
!OCL NOVREL
 !$OMP PARALLEL DO SHARED(WDES1,CH,ISPINOR,EVALUE) PRIVATE(M,MM)
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CPTWFP(MM)*EVALUE
        ENDDO
 !$OMP END PARALLEL DO
     ENDDO

     

   END SUBROUTINE FFTHAMIL

!************************* SUBROUTINE KINHAMIL_VEC *********************
!
!> low level F77 to calculate
!> similar to #KINHAMIL
!> but adds the action of the vector potential on the wavefunction
!>
!> ~~~
!> -(e A p) e(iGr)
!>  |e|A p  e(iGr) = (|e| A hbar) -i nabla e(iGr)
!> (|e| A hbar) G e(iGr)
!> ~~~
!>
!> the coefficients |e| A hbar are included in the stored AVEC
!>
!> @details @ref openmp :
!> the workarrays CWORK1 and CWORK2 acquire a second dimension (size=3).
!> The additional workspace is used to reorder the routine. First fill
!> all three components of CWORK1, then do all FFTs, then add the
!> relevant contributions of all components of CWORK2 to CVR inside
!> a single loop.
!> This improves vectorization of the loops over the elements of these
!> arrays. Aforementioned loops are distributed over all available
!> OpenMP threads.
!
!***********************************************************************

   SUBROUTINE KINHAMIL_VEC( WDES1, GRID, CVR,  LADD, & 
                  DATAKE, IGX, IGY, IGZ, VKPT, AVEC, CWORK1, CWORK2, EVALUE, CPTWFP, CH)

!! USE moffload

     USE prec
     USE mgrid
     USE wave
     USE constant
     IMPLICIT NONE

     TYPE (wavedes1)    WDES1
     TYPE (grid_3d)     GRID
     COMPLEX(q) :: CVR(GRID%MPLWV*WDES1%NRSPINORS) !< usually potential times wavefunction
     LOGICAL    :: LADD                            !< if .TRUE. add results to CH
     REAL(q)    :: DATAKE(WDES1%NGDIM,WDES1%NRSPINORS)
     REAL(q)    :: VKPT(3)                         !< k-point in reciprocal lattice
     INTEGER    :: IGX(WDES1%NGDIM)                !< component of each G vector in direction of first rec lattice vector
     INTEGER    :: IGY(WDES1%NGDIM)                !< component of each G vector in direction of second rec lattice vector
     INTEGER    :: IGZ(WDES1%NGDIM)                !< component of each G vector in direction of third rec lattice vector
     COMPLEX(q)      :: AVEC(GRID%MPLWV,3)     !< magnetic vector potential
     COMPLEX(q) :: CWORK1(WDES1%NRPLWV )
     COMPLEX(q) :: CWORK2(GRID%MPLWV   )
     REAL(q)    :: EVALUE                          ! subtract EVALUE*wavefunction
     COMPLEX(q) :: CPTWFP(WDES1%NRPLWV)                !< wavefunction
     COMPLEX(q) :: CH(WDES1%NRPLWV)                !< result
     COMPLEX(q), PARAMETER :: CI=(0._q,1._q)
! local
     COMPLEX(q) :: CTMP1,CTMP2,CTMP3,CTMP4
     REAL(q)    :: RINPLW
     INTEGER    :: ISPINOR,M,MM
# 2681


     

     RINPLW=2*PI/GRID%NPLWV
 
     DO ISPINOR=0,WDES1%NRSPINORS-1

! component along first axis
!$ACC PARALLEL LOOP PRESENT(WDES1,CWORK1,CPTWFP,IGX) PRIVATE(MM) 
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CWORK1(M)=CPTWFP(MM)*(IGX(M)+VKPT(1))
        ENDDO
! FFT to real space
        CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)
!$ACC PARALLEL LOOP PRESENT(GRID,WDES1,CVR,AVEC,CWORK2) PRIVATE(MM) 
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,GRID%RL%NP
           MM =M+ISPINOR *WDES1%GRID%MPLWV
           CVR(MM)= CVR(MM)+AVEC(M,1) *CWORK2(M) *RINPLW
        ENDDO

! component along second axis
!$ACC PARALLEL LOOP PRESENT(WDES1,CWORK1,CPTWFP,IGY) PRIVATE(MM) 
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CWORK1(M)=CPTWFP(MM)*(IGY(M)+VKPT(2))
        ENDDO
! FFT to real space
        CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)
!$ACC PARALLEL LOOP PRESENT(GRID,WDES1,CVR,AVEC,CWORK2) PRIVATE(MM) 
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,GRID%RL%NP
           MM =M+ISPINOR *WDES1%GRID%MPLWV
           CVR(MM)= CVR(MM)+AVEC(M,2) *CWORK2(M)*RINPLW
        ENDDO

! component along third axis
!$ACC PARALLEL LOOP PRESENT(WDES1,CWORK1,CPTWFP,IGZ) PRIVATE(MM) 
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CWORK1(M)=CPTWFP(MM)*(IGZ(M)+VKPT(3))
        ENDDO
! FFT to real space
        CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)
!$ACC PARALLEL LOOP PRESENT(GRID,WDES1,CVR,AVEC,CWORK2) PRIVATE(MM) 
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,GRID%RL%NP
           MM =M+ISPINOR *WDES1%GRID%MPLWV
           CVR(MM)= CVR(MM)+AVEC(M,3) *CWORK2(M)*RINPLW
        ENDDO

        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CH(1+ISPINOR*WDES1%NGVECTOR),GRID,LADD)
!$ACC PARALLEL LOOP PRESENT(WDES1,CH,CPTWFP,DATAKE) PRIVATE(MM) 
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CPTWFP(MM)*EVALUE+CPTWFP(MM)*DATAKE(M,ISPINOR+1)
        ENDDO
# 2790

     ENDDO

     

   END SUBROUTINE KINHAMIL_VEC


!************************* SUBROUTINE KINHAMIL_TAU *********************
!
! low level F77 to calculate
!> similar to #KINHAMIL
!> but adds the action of   - nabla mu  nabla
!>
!> - nabla mu nabla e(iGr)
!>
!> @details @ref openmp :
!> the workarrays CWORK1 and CWORK2 acquire a second dimension
!> (size=3 and 6, respectively).
!> The additional workspace is used to reorder the routine. First fill
!> all three components of CWORK1, FFT(CWORK1), build the product
!> with the field MU in CWORK2, FFT(CWORK2), and finally add the
!> relevant contributions to CH.
!> This improves vectorization of the loops over the elements of these
!> arrays. Aforementioned loops are distributed over all available
!> OpenMP threads.
!
!***********************************************************************

   SUBROUTINE KINHAMIL_TAU( WDES1, GRID, CVR, LADD, LKINETIC, DATAKE, & 
                  IGX, IGY, IGZ, VKPT, LATT_CUR, MU, CWORK1, CWORK2, EVALUE, CPTWFP, CH)

!! USE moffload

     USE prec
     USE mgrid
     USE wave
     USE lattice
     USE constant
     IMPLICIT NONE

     TYPE (wavedes1)    WDES1
     TYPE (grid_3d)     GRID
     TYPE (latt)        LATT_CUR
     COMPLEX(q) :: CVR(GRID%MPLWV*WDES1%NRSPINORS) !< usually potential times wavefunction
     LOGICAL    :: LADD                            !< if .TRUE. add results to CH
     REAL(q)    :: DATAKE(WDES1%NGDIM,WDES1%NRSPINORS)
     REAL(q)    :: VKPT(3)                         !< k-point in reciprocal lattice
     INTEGER    :: IGX(WDES1%NGDIM)                !< component of each G vector in direction of first rec lattice vector
     INTEGER    :: IGY(WDES1%NGDIM)                !< component of each G vector in direction of second rec lattice vector
     INTEGER    :: IGZ(WDES1%NGDIM)                !< component of each G vector in direction of third rec lattice vector
     COMPLEX(q)      :: MU(GRID%MPLWV*WDES1%NRSPINORS*WDES1%NRSPINORS) ! \mu = d f_xc / d \tau
     REAL(q)    :: EVALUE                          !< subtract EVALUE*wavefunction
     COMPLEX(q) :: CPTWFP(WDES1%NRPLWV)                !< wavefunction
     COMPLEX(q) :: CH(WDES1%NRPLWV)                !< result
     LOGICAL    :: LKINETIC
! local
     INTEGER ISPINOR,ISPINOR_,M,MM,MM_
     REAL(q)    :: RINPLW
     REAL(q)    :: GX,GY,GZ
     COMPLEX(q) :: CWORK1(WDES1%NRPLWV )
     COMPLEX(q) :: CWORK2(GRID%MPLWV   )
!$   COMPLEX(q) :: CTMP1,CTMP2,CTMP3,CTMP4
     COMPLEX(q), PARAMETER :: CI=(0._q,1._q)
# 2858


     

     RINPLW=1.0_q/GRID%NPLWV
     DO ISPINOR=0,WDES1%NRSPINORS-1

        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CH(1+ISPINOR*WDES1%NGVECTOR),GRID,LADD)
        IF (LKINETIC) THEN
!$ACC PARALLEL LOOP PRESENT(WDES1,CH,CPTWFP,DATAKE) PRIVATE(MM) 
!DIR$ IVDEP
!OCL NOVREL
           DO M=1,WDES1%NGVECTOR
              MM=M+ISPINOR*WDES1%NGVECTOR
              CH(MM)=CH(MM)-CPTWFP(MM)*EVALUE+CPTWFP(MM)*DATAKE(M,ISPINOR+1)
           ENDDO
        ELSE
!$ACC PARALLEL LOOP PRESENT(WDES1,CH,CPTWFP) PRIVATE(MM) 
!DIR$ IVDEP
!OCL NOVREL
           DO M=1,WDES1%NGVECTOR
              MM=M+ISPINOR*WDES1%NGVECTOR
              CH(MM)=CH(MM)-CPTWFP(MM)*EVALUE
           ENDDO
        ENDIF

! component along x-axis
!$ACC KERNELS PRESENT(CVR) 
        CVR(ISPINOR*WDES1%GRID%MPLWV+1:(ISPINOR+1)*WDES1%GRID%MPLWV)=0
!$ACC END KERNELS
        DO ISPINOR_=0,WDES1%NRSPINORS-1
!$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CWORK1,CPTWFP) PRIVATE(GX,MM) 
!DIR$ IVDEP
!OCL NOVREL
           DO M=1,WDES1%NGVECTOR
              GX=((IGX(M)+VKPT(1))*LATT_CUR%B(1,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(1,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(1,3))
              MM=M+ISPINOR_*WDES1%NGVECTOR
              CWORK1(M)=CPTWFP(MM)*GX*CITPI
           ENDDO
! FFT to real space
           CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)
!$ACC PARALLEL LOOP PRESENT(GRID,WDES1,CVR,MU,CWORK2) PRIVATE(MM,MM_) 
!DIR$ IVDEP
!OCL NOVREL
           DO M=1,GRID%RL%NP
              MM =M+ISPINOR *WDES1%GRID%MPLWV
              MM_=M+(ISPINOR_+2*ISPINOR)*GRID%MPLWV
              CVR(MM)=CVR(MM)+MU(MM_) *CWORK2(M) *RINPLW
           ENDDO
        ENDDO
! FFT to reciprocal space
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CWORK1(1),GRID,.FALSE.)
!$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CH,CWORK1,CPTWFP)  &
!$ACC PRIVATE(GX,MM)
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           GX=((IGX(M)+VKPT(1))*LATT_CUR%B(1,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(1,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(1,3))
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CWORK1(M)*GX*CITPI
        ENDDO

! component along y-axis
!$ACC KERNELS PRESENT(CVR) 
        CVR(ISPINOR*WDES1%GRID%MPLWV+1:(ISPINOR+1)*WDES1%GRID%MPLWV)=0
!$ACC END KERNELS
        DO ISPINOR_=0,WDES1%NRSPINORS-1
!$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CWORK1,CPTWFP) PRIVATE(GY,MM) 
!DIR$ IVDEP
!OCL NOVREL
           DO M=1,WDES1%NGVECTOR
              GY=((IGX(M)+VKPT(1))*LATT_CUR%B(2,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(2,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(2,3))
              MM=M+ISPINOR_*WDES1%NGVECTOR
              CWORK1(M)=CPTWFP(MM)*GY*CITPI
           ENDDO
! FFT to real space
           CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)
!$ACC PARALLEL LOOP PRESENT(GRID,WDES1,CVR,MU,CWORK2) PRIVATE(MM,MM_) 
!DIR$ IVDEP
!OCL NOVREL
           DO M=1,GRID%RL%NP
              MM =M+ISPINOR *WDES1%GRID%MPLWV
              MM_=M+(ISPINOR_+2*ISPINOR)*GRID%MPLWV
              CVR(MM)=CVR(MM)+MU(MM_) *CWORK2(M)*RINPLW
           ENDDO
        ENDDO
! FFT to reciprocal space
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CWORK1(1),GRID,.FALSE.)
!$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CH,CWORK1,CPTWFP)  &
!$ACC PRIVATE(GY,MM)
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           GY=((IGX(M)+VKPT(1))*LATT_CUR%B(2,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(2,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(2,3))
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CWORK1(M)*GY*CITPI
        ENDDO

! component along z-axis
!$ACC KERNELS PRESENT(CVR) 
        CVR(ISPINOR*WDES1%GRID%MPLWV+1:(ISPINOR+1)*WDES1%GRID%MPLWV)=0
!$ACC END KERNELS
        DO ISPINOR_=0,WDES1%NRSPINORS-1
!$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CWORK1,CPTWFP) PRIVATE(GZ,MM) 
!DIR$ IVDEP
!OCL NOVREL
           DO M=1,WDES1%NGVECTOR
              GZ=((IGX(M)+VKPT(1))*LATT_CUR%B(3,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(3,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(3,3))
              MM=M+ISPINOR_*WDES1%NGVECTOR
              CWORK1(M)=CPTWFP(MM)*GZ*CITPI
           ENDDO
! FFT to real space
           CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)
!$ACC PARALLEL LOOP PRESENT(GRID,WDES1,CVR,MU,CWORK2) PRIVATE(MM,MM_) 
!DIR$ IVDEP
!OCL NOVREL
           DO M=1,GRID%RL%NP
              MM =M+ISPINOR *WDES1%GRID%MPLWV
              MM_=M+(ISPINOR_+2*ISPINOR)*GRID%MPLWV
              CVR(MM)=CVR(MM)+MU(MM_) *CWORK2(M)*RINPLW
           ENDDO
        ENDDO
! FFT to reciprocal space
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CWORK1(1),GRID,.FALSE.)
!$ACC PARALLEL LOOP PRESENT(WDES1,IGX,IGY,IGZ,LATT_CUR,CH,CWORK1,CPTWFP)  &
!$ACC PRIVATE(GZ,MM)
!DIR$ IVDEP
!OCL NOVREL
        DO M=1,WDES1%NGVECTOR
           GZ=((IGX(M)+VKPT(1))*LATT_CUR%B(3,1)+(IGY(M)+VKPT(2))*LATT_CUR%B(3,2)+(IGZ(M)+VKPT(3))*LATT_CUR%B(3,3))
           MM=M+ISPINOR*WDES1%NGVECTOR
           CH(MM)=CH(MM)-CWORK1(M)*GZ*CITPI
        ENDDO
# 3081

     ENDDO

     

   END SUBROUTINE KINHAMIL_TAU


!************************* SUBROUTINE PSCURRENT ************************
!
!> gradient psi contribution to current: Im[psi^*(r) nabla psi(r)]
!> density contribution to current: psi^*(r) * psi(r)
!> for a supplied weighted k-point and psi
!> with psi(r) = e(iGr)
!> we get contributions proportional to (iG) e(iGr)
!
!***********************************************************************

   SUBROUTINE PSCURRENT( WDES1, GRID, CR, & 
                  IGX, IGY, IGZ, VKPT, WEIGHT, AVEC, CWORK1, CWORK2, CPTWFP, J, C)
     USE prec
     USE mgrid
     USE wave
     USE constant
     IMPLICIT NONE

     TYPE (wavedes1)    WDES1
     TYPE (grid_3d)     GRID
     COMPLEX(q) :: CR(GRID%MPLWV*WDES1%NRSPINORS)  !< wavefunction in real space
     REAL(q)    :: VKPT(3)                         !< k-point in reciprocal lattice
     REAL(q)    :: WEIGHT
     INTEGER    :: IGX(WDES1%NGDIM)                !< component of each G vector in direction of first rec lattice vector
     INTEGER    :: IGY(WDES1%NGDIM)                !< component of each G vector in direction of second rec lattice vector
     INTEGER    :: IGZ(WDES1%NGDIM)                !< component of each G vector in direction of third rec lattice vector
     COMPLEX(q)      :: AVEC(GRID%MPLWV,3)     !< magnetic vector potential
     COMPLEX(q) :: CWORK1(WDES1%NRPLWV)
     COMPLEX(q) :: CWORK2(GRID%MPLWV)
     COMPLEX(q) :: CPTWFP(WDES1%NRPLWV)                !< wavefunction
     COMPLEX(q)      :: J(GRID%RL%NP,3)                 !< three component vector storing current density in real space
     COMPLEX(q)      :: C(GRID%RL%NP)                   !< charge density
! local
     INTEGER ISPINOR, M, MM

     DO ISPINOR=0,WDES1%NRSPINORS-1
! component along first axis
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CWORK1(M)=CPTWFP(MM)*(IGX(M)+VKPT(1))*CITPI
        ENDDO
! FFT to real space
        CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

        DO M=1,GRID%RL%NP
           MM =M+ISPINOR *WDES1%GRID%MPLWV
           J(M,1)=J(M,1)+AIMAG(CWORK2(M)*CONJG(CR(MM)))*WEIGHT
        ENDDO

! component along second axis
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CWORK1(M)=CPTWFP(MM)*(IGY(M)+VKPT(2))*CITPI
        ENDDO
! FFT to real space
        CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

        DO M=1,GRID%RL%NP
           MM =M+ISPINOR *WDES1%GRID%MPLWV
           J(M,2)=J(M,2)+AIMAG(CWORK2(M)*CONJG(CR(MM)))*WEIGHT
        ENDDO

! component along third axis
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
           CWORK1(M)=CPTWFP(MM)*(IGZ(M)+VKPT(3))*CITPI
        ENDDO
! FFT to real space
        CALL FFTWAV(WDES1%NGVECTOR,WDES1%NINDPW(1),CWORK2(1),CWORK1(1),GRID)

        DO M=1,GRID%RL%NP
           MM =M+ISPINOR *WDES1%GRID%MPLWV
           J(M,3)=J(M,3)+AIMAG(CWORK2(M)*CONJG(CR(MM)))*WEIGHT
        ENDDO

        DO M=1,GRID%RL%NP
           MM =M+ISPINOR *WDES1%GRID%MPLWV
           C(M)=C(M)+CR(MM)*CONJG(CR(MM))*WEIGHT
        ENDDO
     ENDDO
   END SUBROUTINE PSCURRENT



!************************* SUBROUTINE KINHAMIL_C ***********************
!
!> low level routine for complex eigenvalues
!
!***********************************************************************

   SUBROUTINE KINHAMIL_C( WDES1, GRID, CVR,  LADD, DATAKE, EVALUE, CPTWFP, CH)

!! USE moffload_struct_def

     USE prec
     USE mgrid
     USE wave
     IMPLICIT NONE

     TYPE (wavedes1)    WDES1
     TYPE (grid_3d)     GRID
     COMPLEX(q) :: CVR(GRID%MPLWV*WDES1%NRSPINORS) !< usually potential times wavefunction
     LOGICAL    :: LADD                            !< if .TRUE. add results to CH
     REAL(q)    :: DATAKE(WDES1%NGDIM,WDES1%NRSPINORS)
     COMPLEX(q) :: EVALUE                          !< subtract EVALUE*wavefunction
     COMPLEX(q) :: CPTWFP(WDES1%NRPLWV)                !< wavefunction
     COMPLEX(q) :: CH(WDES1%NRPLWV)                !< result
! local
     INTEGER ISPINOR, M, MM
# 3200


     DO ISPINOR=0,WDES1%NRSPINORS-1
        CALL FFTEXT(WDES1%NGVECTOR,WDES1%NINDPW(1),CVR(1+ISPINOR*WDES1%GRID%MPLWV),CH(1+ISPINOR*WDES1%NGVECTOR),GRID,LADD)
!$ACC PARALLEL LOOP PRESENT(CH,CPTWFP,DATAKE,WDES1) PRIVATE(MM) 
!DIR$ IVDEP
!OCL NOVREL
 !$OMP PARALLEL DO DEFAULT(NONE) SHARED(WDES1,ISPINOR,CH,CPTWFP,EVALUE,DATAKE) PRIVATE(M,MM)
        DO M=1,WDES1%NGVECTOR
           MM=M+ISPINOR*WDES1%NGVECTOR
!          CH(MM)=CH(MM)-CPTWFP(MM)*EVALUE+CPTWFP(MM)*DATAKE(M,ISPINOR+1)
           CH(MM)=CH(MM)+CPTWFP(MM)*(DATAKE(M,ISPINOR+1)-EVALUE)
        ENDDO
 !$OMP END PARALLEL DO
     ENDDO
   END SUBROUTINE KINHAMIL_C

!***********************************************************************
!
!> low level F77 version
!> determine norm of a wavefunction and norm with a metric
!> (if supplied)
!
!***********************************************************************

  SUBROUTINE PW_NORM_WITH_METRIC(WDES1, C, FNORM, FMETRIC, METRIC)

!! USE moffload

    USE wave
    IMPLICIT NONE
    REAL(q) FNORM
    TYPE (wavedes1)    WDES1
    COMPLEX(q) :: C(WDES1%NPL)
    REAL(q), OPTIONAL :: FMETRIC
    REAL(q), OPTIONAL :: METRIC(WDES1%NPL)
! local
    INTEGER ISPINOR, M, MM

    IF (PRESENT (METRIC).AND. PRESENT(FMETRIC)) THEN
!$ACC DATA COPYOUT(FMETRIC,FNORM) 
!$ACC KERNELS 
       FNORM=0
       FMETRIC=0
!$ACC END KERNELS

!$ACC PARALLEL LOOP COLLAPSE(2) REDUCTION(+:FNORM,FMETRIC) PRIVATE(MM) &
!$ACC PRESENT(C,METRIC,FMETRIC,FNORM) 
       DO ISPINOR=0,WDES1%NRSPINORS-1
 !$OMP PARALLEL DO PRIVATE(M,MM) SHARED(WDES1,ISPINOR,C,METRIC) REDUCTION(+:FNORM,FMETRIC)
          DO M=1,WDES1%NGVECTOR
             MM=M+ISPINOR*WDES1%NGVECTOR
             FNORM =FNORM+C(MM)*CONJG(C(MM))
             FMETRIC =FMETRIC+C(MM)*CONJG(C(MM))*METRIC(MM)
          ENDDO
 !$OMP END PARALLEL DO
       ENDDO
!$ACC END DATA
!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON.AND.(.NOT.(ACC_IS_PRESENT(FNORM,1).AND.ACC_IS_PRESENT(FMETRIC,1))))
       CALL M_sum_2(WDES1%COMM_INB, FNORM, FMETRIC)
    ELSE
!$ACC DATA COPYOUT(FNORM) 
!$ACC KERNELS 
       FNORM=0
!$ACC END KERNELS

!$ACC PARALLEL LOOP COLLAPSE(2) REDUCTION(+:FNORM) PRIVATE(MM) &
!$ACC PRESENT(C,FNORM) 
       DO ISPINOR=0,WDES1%NRSPINORS-1
 !$OMP PARALLEL DO PRIVATE(M,MM) SHARED(WDES1,ISPINOR,C) REDUCTION(+:FNORM)
          DO M=1,WDES1%NGVECTOR
             MM=M+ISPINOR*WDES1%NGVECTOR
             FNORM =FNORM+C(MM)*CONJG(C(MM))
          ENDDO
 !$OMP END PARALLEL DO
       ENDDO
!$ACC END DATA
!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON.AND.(.NOT.ACC_IS_PRESENT(FNORM,1)))
       CALL M_sum_d(WDES1%COMM_INB, FNORM, 1)
    ENDIF

  END SUBROUTINE PW_NORM_WITH_METRIC

!***********************************************************************
!
!> low level F77 routine
!> truncate high frequency components if LDELAY is set
!> or if use_eini_without_ldelay is defined
!> first version is scheduled for removal
!
!***********************************************************************


  SUBROUTINE TRUNCATE_HIGH_FREQUENCY( WDES1, C, LDELAY, ENINI)
    USE prec
    USE wave
    IMPLICIT NONE
    TYPE (wavedes1)    WDES1
    COMPLEX(q) C(WDES1%NPL)
    LOGICAL LDELAY     !< usually truncation only during delay phase (LDELAY set)
    REAL(q) ENINI      !< cutoff at which the wavefunction is truncated
! local
    INTEGER ISPINOR

    

    IF (LDELAY.OR.WDES1%LSPIRAL) THEN
       DO ISPINOR=1,WDES1%NRSPINORS
          CALL TRUNCATE_HIGH_FREQUENCY_ONE(WDES1%NGVECTOR, &
               C(1+(ISPINOR-1)*WDES1%NGVECTOR), WDES1%DATAKE(1,ISPINOR), ENINI)
       ENDDO
    END IF

    

  END SUBROUTINE TRUNCATE_HIGH_FREQUENCY


  SUBROUTINE TRUNCATE_HIGH_FREQUENCY_ONE(NP, CPTWFP, DATAKE, ENINI)

!! USE moffload

    USE prec
    IMPLICIT NONE
    INTEGER NP
    COMPLEX(q) :: CPTWFP(NP)
    REAL(q)    :: DATAKE(NP)
    REAL(q)    :: ENINI
! local
    INTEGER N
!$ACC PARALLEL LOOP GANG VECTOR PRESENT(DATAKE,CPTWFP) 
 !$OMP PARALLEL DO PRIVATE(N) SHARED(NP,DATAKE,ENINI,CPTWFP)
    DO N=1,NP
       IF (DATAKE(N)>ENINI) CPTWFP(N)=0
    ENDDO
 !$OMP END PARALLEL DO
  END SUBROUTINE TRUNCATE_HIGH_FREQUENCY_ONE


!************************* SUBROUTINE ECCP_NL_ *************************
!
!> this subroutine calculates the expectation value of
!>    <c| D- evalue Q |cp>
!> where c and cp are two wavefunctions; non local part only
!> for (1._q,0._q) ion only
!
!***********************************************************************

  SUBROUTINE ECCP_NL_(LMDIM,LMMAXC,CDIJ,CQIJ,EVALUE,CPROJ1,CPROJ2,CNL)
!$ACC ROUTINE VECTOR
    USE prec
    IMPLICIT NONE
    COMPLEX(q)  CNL
    INTEGER LMMAXC, LMDIM
    COMPLEX(q) CDIJ(LMDIM,LMDIM),CQIJ(LMDIM,LMDIM)
    REAL(q) EVALUE
    COMPLEX(q) CPROJ1(LMMAXC),CPROJ2(LMMAXC)
! local variables
    INTEGER L, LP

!$ACC LOOP COLLAPSE(2) REDUCTION(+:CNL) VECTOR
    DO L=1,LMMAXC
       DO LP=1,LMMAXC
          CNL=CNL+(CDIJ(LP,L)-EVALUE*CQIJ(LP,L))*CPROJ1(LP)*CONJG(CPROJ2(L))
       ENDDO
    ENDDO
  END SUBROUTINE ECCP_NL_
