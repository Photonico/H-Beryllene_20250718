# 1 "opergrid.F"
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


# 2 "opergrid.F" 2 

!***********************************************************************
! RCS:  $Id: opergrid.F,v 1.3 2003/06/27 13:22:21 kresse Exp kresse $
!
!>  Sums all points on the grid and dump onto screen
!
!***********************************************************************

      SUBROUTINE SUMGRD(CSRC,GRID)
      USE prec
      USE mgrid

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d) GRID
      DIMENSION CSRC(GRID%NGX_rd,GRID%NGY_rd,GRID%NGZ_rd)

      CSUM=0

      DO NZ=1,GRID%NGZ_rd
      DO NY=1,GRID%NGY_rd
      DO NX=1,GRID%NGX_rd
        CSUM=CSUM+CSRC(NX,NY,NZ)
      ENDDO
      ENDDO
      ENDDO
      WRITE(*,*) CSUM
      RETURN
      END

!***********************************************************************
!
!>  Sums all points on the real grid and dump onto screen
!
!***********************************************************************

      SUBROUTINE SUMRGR(CSRC, GRID)
      USE prec
      USE mgrid

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d) GRID
      COMPLEX(q) CSRC(GRID%NGX,GRID%NGY,GRID%NGZ)
      COMPLEX(q) CSUM
      CSUM=0

      DO NZ=1,GRID%NGZ
      DO NY=1,GRID%NGY
      DO NX=1,GRID%NGX
        CSUM=CSUM+CSRC(NX,NY,NZ)
      ENDDO
      ENDDO
      ENDDO
      WRITE(*,*) CSUM
      RETURN
      END

!***********************SUBROUTINE DUMP1******************************
!
!>  This debug-routine dumps a compressed array for (1._q,0._q) K-point
!
!**********************************************************************

      SUBROUTINE DUMP1(NK,NRPLWV,IGX,IGY,IGZ,NPLWKP,C1)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION IGX(NRPLWV,1),IGY(NRPLWV,1),IGZ(NRPLWV,1)
      DIMENSION NPLWKP(1)
      REAL(q)      C1(1)

      DO 50 K=1,NPLWKP(NK)
        WRITE(*,10)K,IGX(K,NK),IGY(K,NK),IGZ(K,NK), &
     & REAL( C1(K) ,KIND=q)
   50 CONTINUE
   10 FORMAT(I3,2X,3I2,6E14.7)
      WRITE(*,*)
      RETURN
      END SUBROUTINE

!***********************SUBROUTINE DUMPA ******************************
!
!>  This debug-routine dumps a uncompressed array for (1._q,0._q) K-Point
!
!**********************************************************************

      SUBROUTINE DUMPA(NK,NGX,NGY,NGZ,CA)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      COMPLEX(q)   CA(NGX,NGY,NGZ)

      NDUMP=MIN(NGX,NGY,NGZ,10)

      DO NZ=1,NDUMP
      WRITE(*,*)'NZ=',NZ
      DO NY=1,NDUMP
        WRITE(*,'(I2,24E10.3)') NY,(ABS(CA(NX,NY,NZ)),NX=1,NDUMP)
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE

!*********************************************************************
!
!> Routine to dump a (NxN) Matrix
!
!*********************************************************************

      SUBROUTINE DDUMP(C,NBANDS)
      USE prec
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      COMPLEX(q) C(NBANDS,NBANDS)

      NBA=MIN(10,NBANDS)

      DO 830 N1=NBANDS-NBA+1,NBANDS
      WRITE(*,1)N1, &
     &         (REAL( C(N1,N2) ,KIND=q) ,N2=NBANDS-NBA+1,NBANDS)
 830  CONTINUE

      WRITE(*,*)
      DO 840 N1=NBANDS-NBA+1,NBANDS
      WRITE(*,2)N1, &
     &         (AIMAG(C(N1,N2)),N2=NBANDS-NBA+1,NBANDS)
 840  CONTINUE

      WRITE(*,*)

    1 FORMAT(1I2,3X,20F9.5)
    2 FORMAT(1I2,3X,20F9.5)
      RETURN
      END SUBROUTINE

