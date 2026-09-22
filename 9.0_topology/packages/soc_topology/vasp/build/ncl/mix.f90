# 1 "mix.F"
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


# 2 "mix.F" 2 
!***********************SUBROUTINE CHCJG ******************************
! RCS:  $Id: mix.F,v 1.1 2000/11/15 08:13:54 kresse Exp $
!
!> Performes a steeptes descent or conjugate gradient
!> step to optimize the chargedensity
!>
!> the method introduced by G. Kresse
!>
!> The routine uses (1._q,0._q) flags to exchange information with the main
!> program
!> ~~~
!>  IMIX       0  initialisation main program supplies the old
!>                chargedensity and old energy\n
!>             1  conventional Kerker-mixing\n
!>  B          recip. lattice  vectors\n
!>  RMST       norm of residual vector
!> ~~~
!> Mind: RMST must be set to 0 before call
!
!**********************************************************************

      SUBROUTINE MIX_SIMPLE(GRIDC,MIX,ISPIN, CHTOT,CHTOTL, &
                       N_MIX_PAW, RHOLM, RHOLM_LAST, B, OMEGA, DE2H, RMST)
      USE prec

      USE mpimy
      USE mgrid
      USE base
      USE constant
      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d) GRIDC
      TYPE (mixing)  MIX

      COMPLEX(q) CHTOT(GRIDC%MPLWV,ISPIN),CHTOTL(GRIDC%MPLWV,ISPIN)
      REAL(q)    B(3,3)
      REAL(q)    RHOLM(N_MIX_PAW,ISPIN),RHOLM_LAST(N_MIX_PAW,ISPIN)
      REAL(q)    OMEGA
      COMPLEX(q), SAVE, ALLOCATABLE :: CHVEL(:,:)
      REAL(q), SAVE, ALLOCATABLE    :: RHOVEL(:,:)
      REAL(q) :: DE2H
      REAL(q) :: MU

      LOGICAL,SAVE :: INI=.TRUE.
      IF (INI .AND. MIX%IMIX==2) THEN
        ALLOCATE(CHVEL(GRIDC%RC%NP,ISPIN)) ;  INI=.FALSE.
        CHVEL=0
        IF (N_MIX_PAW>0) THEN
           ALLOCATE(RHOVEL(N_MIX_PAW,ISPIN))
           RHOVEL=0
        ENDIF
      ENDIF

      DE2H=0

      FAKT=1._q/OMEGA
!=======================================================================
! IMIX==0
! on init just copy CHTOT to CHTOTL
!=======================================================================
      IF (MIX%IMIX==0) THEN
          CHTOTL(1:GRIDC%RC%NP,:)=CHTOT(1:GRIDC%RC%NP,:)
          RHOLM_LAST=RHOLM
!======================================================================
! IMIX==1
! conventional Kerker-mixing
! (IMIX=4 is usually reserved for Broyden, so handle that case as
!  well, in case the routine is called instead of BRMIX)
!======================================================================
      ELSE IF (MIX%IMIX==1 .OR. MIX%IMIX==4) THEN
      AMIX=MIX%AMIX
      BMIX=MIX%BMIX

      DO ISP=1,ISPIN
      FLAM=BMIX**2
      DO K=1,GRIDC%RC%NP
        N1= MOD((K-1),GRIDC%RC%NROW) +1
        NC= (K-1)/GRIDC%RC%NROW+1
        N2= GRIDC%RC%I2(NC)
        N3= GRIDC%RC%I3(NC)

        FACTM=1
        

        GX= (GRIDC%LPCTX(N1)*B(1,1)+GRIDC%LPCTY(N2)*B(1,2)+GRIDC%LPCTZ(N3)*B(1,3))
        GY= (GRIDC%LPCTX(N1)*B(2,1)+GRIDC%LPCTY(N2)*B(2,2)+GRIDC%LPCTZ(N3)*B(2,3))
        GZ= (GRIDC%LPCTX(N1)*B(3,1)+GRIDC%LPCTY(N2)*B(3,2)+GRIDC%LPCTZ(N3)*B(3,3))

        GSQU =(GX**2+GY**2+GZ**2)*TPI*TPI
        IF (N1==1 .AND. N2==1 .AND. N3==1) GSQU=1E30_q

        CHNEW= CHTOTL(K,ISP)+AMIX*GSQU/(GSQU+FLAM)*(CHTOT(K,ISP)-CHTOTL(K,ISP))

        RTMP = FAKT* CONJG((CHTOT(K,ISP)-CHTOTL(K,ISP)))*(CHTOT(K,ISP)-CHTOTL(K,ISP))

        RMST = RMST+ RTMP
        DE2H = DE2H+ RTMP/GSQU

        CHTOTL(K,ISP)=CHTOT(K,ISP)
        CHTOT (K,ISP)=CHNEW
      ENDDO

      AMIX_PAW=AMIX
      DO K=1,N_MIX_PAW
         DEL =(RHOLM(K,ISP)-RHOLM_LAST(K,ISP))
         RMST=RMST+DEL*DEL
         RHO_NEW=RHOLM_LAST(K,ISP)+DEL*AMIX_PAW
         RHOLM_LAST(K,ISP)=RHOLM(K,ISP)
         RHOLM(K,ISP)=RHO_NEW
      ENDDO

      AMIX=MIX%AMIX_MAG
      BMIX=MIX%BMIX_MAG

      ENDDO

      CALL M_sum_2(GRIDC%COMM, RMST, DE2H)
      RMST=SQRT(RMST)
      DE2H=DE2H*EDEPS

!======================================================================
! IMIX==2
! Tshebyshef Kerker-mixing (velocity damping algorithm)
!======================================================================
      ELSE IF (MIX%IMIX==2) THEN
      AMIX=MIX%AMIX
      BMIX=MIX%BMIX
      MU =MIX%AMIN

      DO ISP=1,ISPIN
      FLAM=MIX%BMIX**2
      DO K=1,GRIDC%RC%NP
        N1= MOD((K-1),GRIDC%RC%NROW) +1
        NC= (K-1)/GRIDC%RC%NROW+1
        N2= GRIDC%RC%I2(NC)
        N3= GRIDC%RC%I3(NC)

        FACTM=1
        

        GX= (GRIDC%LPCTX(N1)*B(1,1)+GRIDC%LPCTY(N2)*B(1,2)+GRIDC%LPCTZ(N3)*B(1,3))
        GY= (GRIDC%LPCTX(N1)*B(2,1)+GRIDC%LPCTY(N2)*B(2,2)+GRIDC%LPCTZ(N3)*B(2,3))
        GZ= (GRIDC%LPCTX(N1)*B(3,1)+GRIDC%LPCTY(N2)*B(3,2)+GRIDC%LPCTZ(N3)*B(3,3))

        RTMP = FAKT*  CONJG((CHTOT(K,ISP)-CHTOTL(K,ISP)))*(CHTOT(K,ISP)-CHTOTL(K,ISP))

        RMST = RMST+ RTMP

        GSQU =(GX**2+GY**2+GZ**2)*TPI*TPI
        IF (N1==1 .AND. N2==1 .AND. N3==1) GSQU=1E30_q

        CACC = AMIX*GSQU/(GSQU+FLAM)*(CHTOT(K,ISP)-CHTOTL(K,ISP))
        CHVEL(K,ISP)= ((1-MU/2) * CHVEL(K,ISP) + 2* CACC)/(1+MU/2)
        CHNEW=CHTOTL(K,ISP)+CHVEL(K,ISP)

        DE2H = DE2H+ RTMP/GSQU

        CHTOTL(K,ISP)=CHTOT(K,ISP)
        CHTOT (K,ISP)=CHNEW
      ENDDO
      AMIX_PAW=AMIX
      DO K=1,N_MIX_PAW
         DEL =(RHOLM(K,ISP)-RHOLM_LAST(K,ISP))
         RMST=RMST+DEL*DEL
         RHOVEL(K,ISP)= ((1-MU/2) * RHOVEL(K,ISP) + 2* DEL*AMIX_PAW)/(1+MU/2)

         RHO_NEW=RHOLM_LAST(K,ISP)+RHOVEL(K,ISP)
         RHOLM_LAST(K,ISP)=RHOLM(K,ISP)
         RHOLM(K,ISP)=RHO_NEW
      ENDDO

      AMIX=MIX%AMIX_MAG
      BMIX=MIX%BMIX_MAG
      ENDDO
      CALL M_sum_2(GRIDC%COMM, RMST, DE2H)

      RMST=SQRT(RMST)
      DE2H=DE2H*EDEPS
      ELSE
!======================================================================
! IMIX not sensibly set
!======================================================================
         CALL vtutor%bug("internal error in vasp: MIX_SIMPLE MIX%IMIX is not 0-2; sorry no mixing " &
            // "selected " // str(MIX%IMIX), "mix.F", 187)
      ENDIF

      RETURN
      END
