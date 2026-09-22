# 1 "stufak.F"
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


# 2 "stufak.F" 2 
!************************ SUBROUTINE STUFAK ****************************
! RCS:  $Id: stufak.F,v 1.1 2000/11/15 08:13:54 kresse Exp $
!
!> This subroutine calculates the structure factor on the grid of
!> reciprocal lattice vectors
!> cstrf(g) = sum over ions (-exp(ig.r)) where r is the position of the
!> ion.
!***********************************************************************

      SUBROUTINE STUFAK(GRIDC,T_INFO,CSTRF)

!! USE moffload_struct_def

      USE prec

      USE mpimy
      USE mgrid
      USE poscar
      USE constant
      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (type_info)   T_INFO

      COMPLEX(q) CSTRF(GRIDC%MPLWV,T_INFO%NTYP)

      

! loop over all types of atoms
      NIS=1
      typ: DO NT=1,T_INFO%NTYP
      CALL STUFAK_ONE(GRIDC,T_INFO%NITYP(NT),T_INFO%POSION(1,NIS),T_INFO%VCA(NT),CSTRF(1,NT))
      NIS=NIS+T_INFO%NITYP(NT)

      ENDDO typ

!$ACC UPDATE DEVICE(CSTRF) IF_PRESENT 
      

      RETURN
      END


!************************ SUBROUTINE STUFAK_ONE ************************
!
!> This subroutine calculates the structure factor on the grid of
!> for (1._q,0._q) species (i.e. partial structure factor)
!> cstrf(g) = sum over ions (-exp(ig.r)) where r is the position of the
!> ion.
!
!***********************************************************************

      SUBROUTINE STUFAK_ONE(GRIDC,NIONS,POSION,VCA,CSTRF)
      USE prec

      USE mpimy
      USE mgrid
      USE poscar
      USE constant

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      REAL(q) POSION(3,NIONS)
      REAL(q) VCA

      TYPE (grid_3d)     GRIDC

      COMPLEX(q) CSTRF(GRIDC%RC%NP)

      CSTRF=0


      ion: DO NI=1,NIONS
!=======================================================================
! loop over all grid points
!=======================================================================
# 97

!-----------------------------------------------------------------------
! more envolved version which is faster on most (scalar) machines
! and includes support for parallel machines
!-----------------------------------------------------------------------
         CX =EXP(-CITPI*POSION(1,NI))
         G1 =POSION(1,NI)*(-(GRIDC%NGX/2-1))

         col: DO NC=1,GRIDC%RC%NCOL
            N=(NC-1)*GRIDC%RC%NROW+1

            N2= GRIDC%RC%I2(NC)
            N3= GRIDC%RC%I3(NC)
            G2=POSION(2,NI)*GRIDC%LPCTY(N2)
            G3=POSION(3,NI)*GRIDC%LPCTZ(N3)
            CE=EXP(-CITPI*(G3+G2+G1))*VCA
!DIR$ IVDEP
!$DIR FORCE_VECTOR
!OCL NOVREC
            DO N1P=0,GRIDC%RC%NROW-1
               N1=MOD(N1P+(-(GRIDC%NGX/2-1))+GRIDC%NGX,GRIDC%NGX)
               CSTRF(N+N1)=CSTRF(N+N1)+CE
               CE=CE*CX
            ENDDO
         ENDDO col

!-----------------------------------------------------------------------
!  next ion
!-----------------------------------------------------------------------
      ENDDO ion
# 163

      RETURN
      END

