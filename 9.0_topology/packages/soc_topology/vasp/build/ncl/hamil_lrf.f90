# 1 "hamil_lrf.F"
!#define debug
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


# 3 "hamil_lrf.F" 2 
MODULE hamil_lrf
  USE prec
  USE hamil_lr

  CONTAINS

!************************ SUBROUTINE LRF_COMMUTATOR *********************
!
! this subroutine evaluates
!        |xi_nk>=  ([H(0), r] - e(0) [S(0),r]- e(1) S(0)) | phi(0)_nk>
!               =-i (d/dk [H_k - e(0) S_k- d e_k/d k S(0)]) | phi(0)_nk >
! where H_k (S_k, e_k) is the Hamiltonian (overlap, eigenenergy) in
! dependency of k
! more specifically:
!    |xi>  = - hbar^2/ m_e nabla
!                 + sum_ij | p_i> D(1)_ij <p_j | phi >
!                 + sum_ij | p_i> D(0)_ij - e(0) Q(0)_ij < p_j r| phi >
!                 - sum_ij |r p_i>D(0)_ij - e(0) Q(0)_ij < p_j  | phi >
!                 - sum_ij | p_i> e(1) Q(0)_ij <p_j | phi >
!                 - e(1) |phi>
!
! the last two terms involving e(1) are added to make the xi orthogonal to phi_0
!  (<xi | phi_0> = 0)
! D(1) is usually 0, but could be set to
!  - hbar^2/ m_e < psi_i| nabla |psi_j > - < tilde psi_i| nabla |tilde psi_j >
! for LNABLA ==.TRUE., terms involving <p_i r| are then ommited
! IDIR supplies the cartesian index of r (r_i)
!
! for LI=.TRUE.,  (d/dk [H_k - e(0) S_k] - d e_k/d k S(0)]) | phi(0)_nk >
! is returned
! in this case the change of the eigenvalues is real, whereas in the
! previous case it is complex
!
! set upon return:
! <G | xi>        is returned in WXI%CPTWFP
! <p_j r| phi >   is returned in WXI%CPROJ,
! e(1)            is returned in WXI%CELEN
!
!***********************************************************************

    SUBROUTINE LRF_COMMUTATOR(GRID,INFO,LATT_CUR, &
         NONLR_S, NONLR_D, NONL_S, NONL_D, W0, WXI, WDES, &
          LMDIM, CDIJ, CDIJ0, CQIJ, LNABLA, LI, IDIR, CSHIFT, RMS, ICOUEV, LXI_SET)
      USE prec
      
      USE wave
      USE wave_high
      USE base
      USE lattice
      USE mpimy
      USE mgrid
    
      USE nonl_high
      USE hamil
      USE constant
      USE wave_mpi
      USE tutor, ONLY: vtutor

      IMPLICIT NONE

      TYPE (grid_3d)     GRID
      TYPE (info_struct) INFO
      TYPE (latt)        LATT_CUR
      TYPE (nonlr_struct) NONLR_S, NONLR_D
      TYPE (nonl_struct) NONL_S, NONL_D
      TYPE (wavespin)    W0, WXI, W1
      TYPE (wavedes)     WDES

      INTEGER LMDIM
      COMPLEX(q) CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      COMPLEX(q) CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      COMPLEX(q) CDIJ0(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      REAL(q) RMS        ! magnitude of the residual vector
      INTEGER ICOUEV     ! number of H | phi> evaluations
      LOGICAL LNABLA
      LOGICAL LI
      INTEGER IDIR
      REAL(q) CSHIFT
      LOGICAL, OPTIONAL :: LXI_SET
! local work arrays and structures
      TYPE (wavedes1)  WDES1                     ! descriptor for (1._q,0._q) k-point
      TYPE (wavefun1), TARGET :: W0_1(WDES%NSIM) ! current wavefunction
      TYPE (wavefun1)  WTMP(WDES%NSIM)           ! see below

      COMPLEX(q),ALLOCATABLE:: CF(:,:)
      COMPLEX(q), TARGET, ALLOCATABLE ::  CPROJ(:,:)

      INTEGER :: NSIM                  ! number of bands treated simultaneously
      INTEGER :: NODE_ME, IONODE
      INTEGER :: NP, ISP, NK, NPL, NGVECTOR, NB_DONE, N, IDUMP, ISPINOR, LD, M, MM
      REAL(q) :: FNORM
      COMPLEX(q) :: ORTH
      INTEGER :: NB(WDES%NSIM)         ! contains a list of bands currently optimized
      REAL(q) :: EVALUE(WDES%NSIM)     ! eigenvalue during optimization
      LOGICAL :: LDO(WDES%NSIM)        ! band finished
      LOGICAL :: LSTOP,LNEWB
      COMPLEX(q) ::C
      COMPLEX(q) ::CI                ! imaginary

      

      EVALUE=0.0_q
      NODE_ME=0
      IONODE =0

      NODE_ME=WDES%COMM%NODE_ME
      IONODE =WDES%COMM%IONODE


      IF (LI) THEN
         CI=(0.0_q,1.0_q)
# 117

      ELSE
         CI=(1.0_q,0.0_q)
      ENDIF
    
!=======================================================================
!  INITIALISATION:
! maximum  number of bands treated simultaneously
!=======================================================================
      RMS   =0
      ICOUEV=0
      NSIM=WDES%NSIM

      ALLOCATE(CF(WDES%NRPLWV,NSIM),CPROJ(WDES%NPROD,NSIM))
! CF=(0.0_q, 0.0_q); CPROJ=0 ! Breaks bulk_InP_SOC_G0W0_sym (gnu), may be copied uninitialized below.
      LD=WDES%NRPLWV

      DO NP=1,NSIM
         ALLOCATE(W0_1(NP)%CR(GRID%MPLWV*WDES%NRSPINORS))
      ENDDO

      WXI%CPROJ=0
!=======================================================================
      spin:    DO ISP=1,WDES%ISPIN
      kpoints: DO NK=1,WDES%NKPTS

      IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE

      CALL SETWDES(WDES,WDES1,NK); CALL SETWGRID_OLD(WDES1,GRID) 

      NPL=WDES1%NPL
      NGVECTOR=WDES1%NGVECTOR

      IF (INFO%LREAL) THEN
        CALL PHASER (GRID,LATT_CUR,NONLR_S,NK,WDES)
        CALL PHASERR(GRID,LATT_CUR,NONLR_D,NK,WDES,IDIR)
      ELSE
        CALL PHASE(WDES,NONL_S,NK)
        NONL_D%NK=NONL_S%NK        ! uses the same phasefactor array
      ENDIF

!=======================================================================
      NB_DONE=0                   ! index of the bands allready optimised
      bands: DO
        NB=0                      ! empty the list of bands, which are optimized currently
!
!  check the NB list, whether there is any empty slot
!  fill in a not yet optimized wavefunction into the slot
!
        IDUMP=0

        newband: DO NP=1,NSIM
        LNEWB=.FALSE.

        IF (NB_DONE<WXI%WDES%NBANDS) THEN
           NB_DONE=NB_DONE+1; LNEWB=.TRUE.
        ENDIF

        IF (LNEWB) THEN
           N     =NB_DONE
           NB(NP)=NB_DONE
           ICOUEV=ICOUEV+1

           CALL SETWAV(W0,W0_1(NP),WDES1,N,ISP)  ! fill band N into W0_1(NP)

# 184


           IF (NODE_ME /= IONODE) IDUMP=0


! fft to real space
           DO ISPINOR=0,WDES%NRSPINORS-1
              CALL FFTWAV(NGVECTOR, WDES%NINDPW(1,NK),W0_1(NP)%CR(1+ISPINOR*GRID%MPLWV),W0_1(NP)%CPTWFP(1+ISPINOR*NGVECTOR),GRID)
           ENDDO
! WTMP is identical to W0_1, except for CPROJ entry
! which will contain the derivative of the projectors
! w.r.t. k after the call HAMILMU_LRF
           WTMP(NP)=W0_1(NP)
           WTMP(NP)%CPROJ => CPROJ(:,NP)

           EVALUE(NP)=W0_1(NP)%CELEN
        ENDIF
        ENDDO newband

!=======================================================================
! if the NB list is now empty end the bands DO loop
!=======================================================================
        LSTOP=.TRUE.
        LDO  =.FALSE.
        DO NP=1,NSIM
           IF ( NB(NP) /= 0 ) THEN
              LSTOP  =.FALSE.
              LDO(NP)=.TRUE.     ! band not finished yet
           ENDIF
        ENDDO
        IF (LSTOP) EXIT bands
!=======================================================================
! determine gradient and store it
!=======================================================================
        
        DO NP=1,NSIM
           N=NB(NP); IF (.NOT. LDO(NP)) CYCLE
           IF (PRESENT (LXI_SET)) THEN
              CF(:,NP)=WXI%CPTWFP(:,N,NK,ISP)
           ELSE
              CF(:,NP)=0
           ENDIF
        ENDDO

!  store H | psi > temporarily
!  to have uniform stride for result array
        CALL HAMILTMU_COMMUTATOR(WDES1,W0_1,WTMP,NONLR_S,NONLR_D,NONL_S,NONL_D, &
                & GRID,  INFO%LREAL, EVALUE, LATT_CUR, &
                & LMDIM,CDIJ(1,1,1,ISP),CDIJ0(1,1,1,ISP),CQIJ(1,1,1,ISP), &
                & CF(1,1),LD, NSIM, LDO, WDES%VKPT(1:3,NK), LNABLA, IDIR, CSHIFT)

        i2: DO NP=1,NSIM
           N=NB(NP); IF (.NOT. LDO(NP)) CYCLE i2
           ORTH =0
           FNORM=0
           DO ISPINOR=0,WDES%NRSPINORS-1
!$OMP PARALLEL DO SHARED(NGVECTOR,ISPINOR,CF,NP,W0_1,CI,WDES,INFO,WXI,NK) &
!$OMP PRIVATE(M,MM) REDUCTION(+:FNORM,ORTH)
           DO M=1,NGVECTOR
              MM=M+ISPINOR*NGVECTOR
              CF(MM,NP)=(CF(MM,NP)-W0_1(NP)%CELEN*W0_1(NP)%CPTWFP(MM))*CI
              IF (WDES%LSPIRAL.AND.WDES%DATAKE(M,ISPINOR+1,NK)>INFO%ENINI) CF(MM,NP)=0
              FNORM =FNORM+CF(MM,NP)*CONJG(CF(MM,NP))
              ORTH  =ORTH +(CF(MM,NP)*CONJG(W0_1(NP)%CPTWFP(MM)))
              WXI%CPTWFP(MM,N,NK,ISP)=CF(MM,NP)
           ENDDO
!$OMP END PARALLEL DO
           WXI%CELEN(N,NK,ISP)  =W0_1(NP)%CELEN*CI
           WXI%CPROJ(:,N,NK,ISP)=WTMP(NP)%CPROJ*CI

           ENDDO
           CALL M_sum_d(WDES%COMM_INB, FNORM, 1)
           CALL M_sum_z(WDES%COMM_INB, ORTH, 1)

           IF (ABS(ORTH)>1E-4) THEN
              WRITE(0,*)'LRF_COMMUTATOR internal error: the vector H(1)-e(1) S(1) |phi(0)> is not orthogonal to |phi(0)>',ORTH
!             CALL M_stop('VASP aborting ...'); stop
           ENDIF

           IF (IDUMP==2) WRITE(*,'(I3,E11.4,"R ",2E11.4,"O ",2E11.4,"E ",2E14.7)') N,SQRT(ABS(FNORM)),ORTH,WXI%CELEN(N,NK,ISP),WXI%FERWE(N,NK,ISP)
           RMS=RMS+WDES%RSPIN*WDES%WTKPT(NK)*W0%FERWE(N,NK,ISP)*SQRT(ABS(FNORM))/WDES%NB_TOT

!=======================================================================
! move onto the next block of bands
!=======================================================================
        ENDDO i2
!=======================================================================
      ENDDO bands
      ENDDO kpoints
      ENDDO spin
!=======================================================================
      CALL M_sum_d(WDES%COMM_INTER, RMS, 1)
      CALL M_sum_d(WDES%COMM_KINTER, RMS, 1)

      CALL M_sum_i(WDES%COMM_INTER, ICOUEV ,1)
      CALL M_sum_i(WDES%COMM_KINTER, ICOUEV ,1)

      DO NP=1,NSIM
         DEALLOCATE(W0_1(NP)%CR)
      ENDDO
      DEALLOCATE(CF, CPROJ)

      

      RETURN
    END SUBROUTINE LRF_COMMUTATOR


!************************* SUBROUTINE HAMILTMU_COMMUTATOR *************
!
! this subroutine evaluates |xi> at a selected k-point
!    |xi>  = - hbar^2/ m_e nabla
!                 + sum_ij | p_i> D(1)_ij <p_j | phi >
!                 + sum_ij | p_i> D(0)_ij - e(0) Q(0)_ij < p_j r| phi >
!                 - sum_ij |r p_i>D(0)_ij - e(0) Q(0)_ij < p_j  | phi >
!                 - sum_ij | p_i> e(1) Q(0)_ij <p_j | phi >
!                [- e(1) |phi> (see NOTE below)]
!          = -i d/dk [H(k) - e(0) S(k)] | phi >
!
! the  wavefunctions must be supplied in reciprocal space C and real
! space CR
!
! D(1) change of strength parameter with k    (usually (0._q,0._q))
! D(0) is the original strength               (CDIJ0)
! e(0) is the (0._q,0._q) order eigenvalue           (EVALUE0)
! e(1) is the first order change of the eigen energy
!
! e(1) is evaluated during the calculation of |xi>:
! e(1) = <phi| - hbar^2/ m_e nabla  | phi > + sum_ij <phi | p_i> D(1)_ij <p_j | phi >
!        +  sum_ij <phi | p_i > D(0)_ij - e(0) Q(0)_ij <p_j r| phi > -
!        -  sum_ij <phi |r p_i >D(0)_ij - e(0) Q(0)_ij <p_j | phi >
!
! IDIR supplies the cartesian index of (r_i)
! for LNABLA ==.TRUE., terms involving <p_i r| are ommited
!
! NOTE: the calling routine has to subtract e(1) |phi> to obtain the
! correct vector xi
!
!***********************************************************************

    SUBROUTINE HAMILTMU_COMMUTATOR( &
         WDES1,W0_1,WTMP,NONLR_S,NONLR_D,NONL_S,NONL_D, &
         GRID, LREAL, EVALUE0, LATT_CUR, &
         LMDIM,CDIJ,CDIJ0,CQIJ, &
         CH,LD, NSIM, LDO, VKPT, LNABLA, IDIR, CSHIFT)
      USE prec
      USE constant
      USE mpimy
      USE mgrid
      USE wave
      
      USE nonl_high
      USE lattice
      USE hamil

      IMPLICIT NONE


      INTEGER NSIM,NP,LD
      INTEGER LMDIM, NGVECTOR, ISPINOR, ISPINOR_, MM, MM_ 
      TYPE (grid_3d)     GRID
      TYPE (nonlr_struct) NONLR_S, NONLR_D
      TYPE (nonl_struct)  NONL_S, NONL_D
      TYPE (wavefun1)    W0_1(NSIM)
      TYPE (wavefun1)    WTMP(NSIM)
      TYPE (wavedes1)    WDES1
      TYPE (latt)        LATT_CUR

      COMPLEX(q)    CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS), &
                 CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS), &
                 CDIJ0(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
      COMPLEX(q) CH(LD,NSIM)
      REAL(q)    EVALUE0(NSIM)
      COMPLEX(q) EVALUE0_(NSIM)
      COMPLEX(q) EVALUE1(NSIM)
      LOGICAL LREAL
      LOGICAL LDO(NSIM)
      LOGICAL LNABLA
      INTEGER IDIR
      REAL(q) VKPT(3)
      REAL(q) CSHIFT
! local variables
      REAL(q) RINPLW; INTEGER M
      COMPLEX(q), ALLOCATABLE :: CWORK1(:,:)
      REAL(q)      ::    G1,G2,G3,GC(LD)
      COMPLEX(q)   ::    CE
      INTEGER      ::    LMBASE, LMBASE_, NIS, LMMAXC, NI, L, LP, NT
      REAL(q)      ::    WEIGHT

      

      EVALUE1=0

      ALLOCATE(CWORK1(GRID%MPLWV*WDES1%NRSPINORS,NSIM)) 
      RINPLW=1._q/GRID%NPLWV
      NGVECTOR=WDES1%NGVECTOR
      EVALUE0_=EVALUE0+CMPLX(0.0_q,CSHIFT,q)
!=======================================================================
! calculate the local contribution (result in CWORK1)
!=======================================================================
      CE=0
!$OMP PARALLEL DO PRIVATE(NP)
      DO NP=1,NSIM
         IF ( LDO(NP) ) THEN
            CWORK1(:,NP) =0
         ENDIF
      ENDDO
!$OMP END PARALLEL DO
!=======================================================================
! add term   [- hbar^2/ 2 m_e laplace , r] = - hbar^2/ m_e nabla
!
! using  phi(r) = sum_G C_G  e^ i(G+k)r
!
! the term becomes  -i hbar^2/ m_e (G+k) C_G which is the change of the
! kinetic energy operator with respect to k times -i
!   -i nabla_k ( hbar^2/ 2 m_e (G+k)^2)
!
!=======================================================================
!$OMP PARALLEL DO &
!$OMP SHARED(NGVECTOR,WDES1,VKPT,LATT_CUR,GC) &
!$OMP PRIVATE(M,G1,G2,G3)
      DO M=1,NGVECTOR
         G1=WDES1%IGX(M)+VKPT(1)
         G2=WDES1%IGY(M)+VKPT(2)
         G3=WDES1%IGZ(M)+VKPT(3)
         GC(M)=(G1*LATT_CUR%B(IDIR,1)+G2*LATT_CUR%B(IDIR,2)+G3*LATT_CUR%B(IDIR,3))*(TPI*(2*HSQDTM))
      ENDDO
!$OMP END PARALLEL DO
      DO NP=1,NSIM
         IF ( LDO(NP) ) THEN
            CE=0
            DO ISPINOR =0,WDES1%NRSPINORS-1
!$OMP PARALLEL DO &
!$OMP SHARED(NGVECTOR,ISPINOR,NP,W0_1,CH,GC) &
!$OMP PRIVATE(M,MM) &
!$OMP REDUCTION(+:CE)
               DO M=1,NGVECTOR
                  MM=M+ISPINOR*NGVECTOR
                  CH(MM,NP)=CH(MM,NP)+(W0_1(NP)%CPTWFP(MM)*GC(M))*(0._q,-1._q)
                  CE=CE+CONJG(W0_1(NP)%CPTWFP(MM))*CH(MM,NP)
               ENDDO
!$OMP END PARALLEL DO
            ENDDO
            CALL M_sum_z(WDES1%COMM_INB, CE, 1)
            W0_1(NP)%CELEN=(CE)
         ENDIF
      ENDDO
!=======================================================================
! non-local contribution in real-space
!=======================================================================
      IF (LREAL) THEN
         IF (.NOT. LNABLA) THEN
! contribution - r_idir | p_i > D_ij - e(0) Q_ij < p_j |
            IF (CSHIFT==0) THEN
               CALL RACCMU(NONLR_D,WDES1,W0_1, LMDIM,CDIJ0,CQIJ,EVALUE0,CWORK1,GRID%MPLWV*WDES1%NRSPINORS, NSIM, LDO)
            ELSE
               CALL RACCMU_C(NONLR_D,WDES1,W0_1, LMDIM,CDIJ0,CQIJ,CONJG(EVALUE0_),CWORK1,GRID%MPLWV*WDES1%NRSPINORS, NSIM, LDO)
            ENDIF
            DO NP=1,NSIM
               IF ( LDO(NP) ) CWORK1(:,NP)=-CWORK1(:,NP)
            ENDDO

! contribution | p_i > D_ij - e(0) Q_ij <  p_j | r_idir
            CALL RPROMU(NONLR_D,WDES1,WTMP, NSIM, LDO)
            IF (CSHIFT==0) THEN
               CALL RACCMU(NONLR_S,WDES1,WTMP, LMDIM,CDIJ0,CQIJ,EVALUE0,CWORK1,GRID%MPLWV*WDES1%NRSPINORS, NSIM, LDO)
            ELSE
               CALL RACCMU_C(NONLR_S,WDES1,WTMP, LMDIM,CDIJ0,CQIJ,EVALUE0_,CWORK1,GRID%MPLWV*WDES1%NRSPINORS, NSIM, LDO)
            ENDIF
         ELSE
!$OMP PARALLEL DO PRIVATE(NP)
            DO NP=1,NSIM
               IF ( LDO(NP) ) THEN
                  CWORK1(:,NP)=0
                  WTMP(NP)%CPROJ=0
               ENDIF
            ENDDO
!$OMP END PARALLEL DO
         ENDIF

! non local contributions to e(1)
         DO NP=1,NSIM
            IF ( LDO(NP) ) THEN
               CALL ECCP_NL_ALL(WDES1,W0_1(NP),WTMP(NP),CDIJ0,CQIJ,EVALUE0(NP),CE)
               W0_1(NP)%CELEN=W0_1(NP)%CELEN-CE
               CALL ECCP_NL_ALL(WDES1,WTMP(NP),W0_1(NP),CDIJ0,CQIJ,EVALUE0(NP),CE)
               W0_1(NP)%CELEN=W0_1(NP)%CELEN+CE
               CALL ECCP_NL_ALL(WDES1,W0_1(NP),W0_1(NP),CDIJ,CQIJ,0.0_q,CE)
               W0_1(NP)%CELEN=W0_1(NP)%CELEN+CE
! should be purely imaginary
               W0_1(NP)%CELEN=CMPLX(0.0_q,AIMAG(W0_1(NP)%CELEN),q)
               EVALUE1(NP)=W0_1(NP)%CELEN
            ENDIF
         ENDDO

! contribution | p_i > d D_ij / d k - e(1) Q_ij < p_j |
         CALL RACCMU_C(NONLR_S,WDES1,W0_1, LMDIM,CDIJ,CQIJ,EVALUE1,CWORK1,GRID%MPLWV*WDES1%NRSPINORS, NSIM, LDO)
!=======================================================================
! calculate the non local contribution in reciprocal space
!=======================================================================
      ELSE
         DO NP=1,NSIM
         IF ( LDO(NP) ) THEN
            IF (.NOT. LNABLA) THEN
! contribution - | -i d p_i/ d k> D_ij - epsilon Q_ij < p_j |
               CH(:,NP)=-CH(:,NP)
               IF (CSHIFT==0) THEN
                  CALL VNLACC_ADD(NONL_D,W0_1(NP),CDIJ0,CQIJ,1,EVALUE0(NP),CH(:,NP))
               ELSE
                  CALL VNLACC_ADD_C(NONL_D,W0_1(NP),CDIJ0,CQIJ,1,CONJG(EVALUE0_(NP)),CH(:,NP))
               ENDIF
               CH(:,NP)=-CH(:,NP)
! contribution | p_i > D_ij - epsilon Q_ij < -i d p_j/ d k |
               CALL PROJ1(NONL_D,WDES1, WTMP(NP))
               IF (CSHIFT==0) THEN
                  CALL VNLACC_ADD(NONL_S,WTMP(NP),CDIJ0,CQIJ,1,EVALUE0(NP),CH(:,NP))
               ELSE
                  CALL VNLACC_ADD_C(NONL_S,WTMP(NP),CDIJ0,CQIJ,1,EVALUE0_(NP),CH(:,NP))
               ENDIF
! non local contributions to e(1)
               CALL ECCP_NL_ALL(WDES1,W0_1(NP),WTMP(NP),CDIJ0,CQIJ,EVALUE0(NP),CE)
               W0_1(NP)%CELEN=W0_1(NP)%CELEN-CE
               CALL ECCP_NL_ALL(WDES1,WTMP(NP),W0_1(NP),CDIJ0,CQIJ,EVALUE0(NP),CE)
               W0_1(NP)%CELEN=W0_1(NP)%CELEN+CE
            ELSE
               WTMP(NP)%CPROJ=0
            ENDIF

            CALL ECCP_NL_ALL(WDES1,W0_1(NP),W0_1(NP),CDIJ,CQIJ,0.0_q,CE)
            W0_1(NP)%CELEN=W0_1(NP)%CELEN+CE

! should be purely imaginary
            W0_1(NP)%CELEN=CMPLX(0.0_q,AIMAG(W0_1(NP)%CELEN),q)
            EVALUE1(NP)=W0_1(NP)%CELEN

! contribution | p_i > d D_ij - e(1) Q_ij < p_j |
            CALL VNLACC_ADD_C(NONL_S,W0_1(NP),CDIJ,CQIJ,1,EVALUE1(NP),CH(:,NP))
         ENDIF
         ENDDO
      ENDIF
      DO NP=1,NSIM
         IF ( LDO(NP) ) THEN
            DO ISPINOR=0,WDES1%NRSPINORS-1
               CALL FFTEXT(NGVECTOR,WDES1%NINDPW(1),CWORK1(1+ISPINOR*WDES1%GRID%MPLWV,NP),CH(1+ISPINOR*NGVECTOR,NP),GRID,.TRUE.)
            ENDDO
         ENDIF
      ENDDO

      DEALLOCATE(CWORK1)

      

      RETURN
    END SUBROUTINE HAMILTMU_COMMUTATOR



!************************ SUBROUTINE LRF_HAMIL **************************
!
! this subroutine evaluates
!
!    |xi>  =  (H_0(1) -e(0) S(1)) |phi> + (H_sc(1) - e(1) S(0)) | phi>
!
! for a general perturbation, where H_0(1) -e(0) S(1) |phi> is supplied
! as input an the array in WXI
!
! H_sc(1), e(1) are the first order change of the
! Hamiltonian, eigenenergies and overlap, respectively
! More specifically for the PAW Hamiltonian this becomes:
!
!    |xi>  = V_sc(1) + sum_ij | p_i> D_sc(1)_ij <p_j | phi >
!                    - sum_ij | p_i> e(1) Q(0)_ij <p_j | phi >
!                    - e(1) |phi>
!                    + |xi(0)> + sum_i xi(0)_i | p_i>
!
! xi(0)_i is usually the change of the wave function character with respect
! to i k (see LRF_COMMUTATOR)
! roughly speeking this is essentially
!
!  xi(0)_i= -i d/dk <p_i |phi>= -i(<p_i | d/dk phi>+< d/dk p_i | phi>)
!                               - tau < p_i | phi_nk>
! see comments at the bottom of LRF_RPHI
!
!***********************************************************************

  SUBROUTINE LRF_HAMIL( &
       GRID,INFO,LATT_CUR,NONLR_S,NONL_S,W0,WXI,WDES,LMDIM,CDIJ,CQIJ,SV,RMS,ICOUEV)

!! USE moffload

    USE prec
      
    USE wave
    USE wave_high
    USE base
    USE lattice
    USE mpimy
    USE mgrid
    
    USE nonl_high
    USE hamil
    USE constant
    USE wave_mpi
    USE string, ONLY: str
    USE tutor, ONLY: vtutor

    IMPLICIT NONE

    TYPE (grid_3d)      GRID
    TYPE (info_struct)  INFO
    TYPE (latt)         LATT_CUR
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct)  NONL_S
    TYPE (wavespin)     W0,WXI
    TYPE (wavedes)      WDES

    COMPLEX(q)   SV(GRID%MPLWV,WDES%NCDIJ) ! local potential
    INTEGER LMDIM
    COMPLEX(q) CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    COMPLEX(q) CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    REAL(q) RMS        ! magnitude of the residual vector
    INTEGER ICOUEV     ! number of H | phi> evaluations
! local work arrays and structures
    TYPE (wavedes1)  WDES1                     ! descriptor for (1._q,0._q) k-point
    TYPE (wavefun1), TARGET :: W0_1(WDES%NSIM) ! current wavefunction
    TYPE (wavefun1)  WXI0(WDES%NSIM)           ! see below

    COMPLEX(q) :: ORTH (WDES%NSIM), ORTH_
    REAL(q)    :: FNORM(WDES%NSIM),FNORM_

    COMPLEX(q),ALLOCATABLE:: CF(:,:)
    INTEGER :: NSIM                  ! number of bands treated simultaneously
    INTEGER :: NODE_ME, IONODE
    INTEGER :: NP, ISP, NK, NGVECTOR, NB_DONE, N, IDUMP, ISPINOR, LD, M, MM
    INTEGER :: NB(WDES%NSIM)         ! contains a list of bands currently optimized
    LOGICAL :: LDO(WDES%NSIM)        ! band finished
    LOGICAL :: LSTOP,LNEWB

      

# 637


      NODE_ME=0
      IONODE =0

      NODE_ME=WDES%COMM%NODE_ME
      IONODE =WDES%COMM%IONODE

!=======================================================================
!  INITIALISATION:
! maximum  number of bands treated simultaneously
!=======================================================================
      RMS   =0
      ICOUEV=0
      NSIM=WDES%NSIM


      ALLOCATE(CF(WDES%NRPLWV,NSIM))
!$ACC ENTER DATA CREATE(CF,WXI0(:),W0_1(:),WDES1,ORTH,FNORM)  
      LD=WDES%NRPLWV

      CALL SETWDES(WDES, WDES1, 0)

      DO NP=1,NSIM
         CALL NEWWAV(WXI0(NP),WDES1,.FALSE.)
         CALL NEWWAV_R(W0_1(NP),WDES1)
      ENDDO

!=======================================================================
      spin:    DO ISP=1,WDES%ISPIN
      kpoints: DO NK=1,WDES%NKPTS

      IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE

      CALL SETWDES(WDES,WDES1,NK)

      NGVECTOR=WDES1%NGVECTOR

      IF (INFO%LREAL) THEN
        CALL PHASER(GRID,LATT_CUR,NONLR_S,NK,WDES)
      ELSE
        CALL PHASE(WDES,NONL_S,NK)
      ENDIF

!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)

!=======================================================================
      NB_DONE=0                 ! index of the bands already optimised
      bands: DO

        NB=0                    ! empty the list of bands, which are optimized currently
!
!  check the NB list, whether there is any empty slot
!  fill in a not yet optimized wavefunction into the slot
!
        IDUMP=0

        

        newband: DO NP=1,NSIM
        LNEWB=.FALSE.

        IF (NB_DONE<WXI%WDES%NBANDS ) THEN
           NB_DONE=NB_DONE+1; LNEWB=.TRUE.
        ENDIF

        IF (LNEWB) THEN

           

           N     =NB_DONE
           NB(NP)=NB_DONE
           ICOUEV=ICOUEV+1

           CALL SETWAV(W0,W0_1(NP),WDES1,N,ISP)  ! fill band N into W0_1(NP)

           CALL ZCOPY(SIZE(WXI%CPTWFP,   1), WXI%CPTWFP   (1,N,NK,ISP), 1, WXI0(NP)%CPTWFP   (1), 1)
           CALL ZCOPY(SIZE(WXI%CPROJ,1), WXI%CPROJ(1,N,NK,ISP), 1, WXI0(NP)%CPROJ(1), 1)
           

# 719


           IF (NODE_ME /= IONODE) IDUMP=0


! fft to real space
           DO ISPINOR=0,WDES%NRSPINORS-1
              CALL FFTWAV(NGVECTOR, WDES%NINDPW(1,NK),W0_1(NP)%CR(1+ISPINOR*GRID%MPLWV),W0_1(NP)%CPTWFP(1+ISPINOR*NGVECTOR),GRID)
           ENDDO
        ENDIF
        ENDDO newband

        

!=======================================================================
! if the NB list is now empty end the bands DO loop
!=======================================================================
        LSTOP=.TRUE.
        LDO  =.FALSE.
        DO NP=1,NSIM
           IF ( NB(NP) /= 0 ) THEN
              LSTOP  =.FALSE.
              LDO(NP)=.TRUE.     ! band not finished yet
           ENDIF
        ENDDO
        IF (LSTOP) EXIT bands

!=======================================================================
! determine gradient and store it
!=======================================================================
!  store H | psi > temporarily
!  to have uniform stride for result array
        CALL HAMILTMU_LRF(WDES1,W0_1,WXI0,NONLR_S,NONL_S, &
                & GRID,  INFO%LREAL, LATT_CUR, &
                & LMDIM,CDIJ(1,1,1,ISP),CQIJ(1,1,1,ISP), &
                & SV(1,ISP), CF(1,1),LD, NSIM, LDO)

        

        i2: DO NP=1,NSIM
            N=NB(NP); IF (.NOT. LDO(NP)) CYCLE i2
# 763

            FNORM_=0
            ORTH_ =0

!$OMP PARALLEL DO COLLAPSE(2) SHARED(NGVECTOR,CF,NP,W0_1,WDES,INFO,WXI) PRIVATE(ISPINOR,M,MM) REDUCTION(+:FNORM_,ORTH_)
            DO ISPINOR=0,WDES%NRSPINORS-1
               DO M=1,NGVECTOR
                  MM=M+ISPINOR*NGVECTOR
                  CF(MM,NP)=CF(MM,NP)-W0_1(NP)%CELEN*W0_1(NP)%CPTWFP(MM)
                  IF (WDES%LSPIRAL.AND.WDES%DATAKE(M,ISPINOR+1,NK)>INFO%ENINI) CF(MM,NP)=0
                  FNORM_ =FNORM_+CF(MM,NP)*CONJG(CF(MM,NP))
                  ORTH_  =ORTH_ +(CF(MM,NP)*CONJG(W0_1(NP)%CPTWFP(MM)))
                  WXI%CPTWFP(MM,N,NK,ISP)=CF(MM,NP)
                ENDDO
            ENDDO
!$OMP END PARALLEL DO

            CALL M_sum_d(WDES%COMM_INB, FNORM_, 1)
            CALL M_sum_z(WDES%COMM_INB, ORTH_ , 1)

            FNORM(NP)=FNORM_ ; ORTH(NP)=ORTH_


! the projector functions do not change when an external field is applied
!$ACC KERNELS PRESENT(WXI) 
            WXI%CPROJ(:,N,NK,ISP)=0
!$ACC END KERNELS
            WXI%CELEN(N,NK,ISP)=W0_1(NP)%CELEN

!=======================================================================
! move onto the next block of bands
!=======================================================================
        ENDDO i2

! mM: this could be used to replace the WAIT(ACC_ASYNC_Q) below as soon
!     as ACC_SYNC_ASYNC_Q can use cudaDeviceSynchronize (which for the
!     moment seems to be buggy).
!!#ifdef ACC_OFFLOAD
!!        CALL ACC_SYNC_ASYNC_Q
!!#endif

        DO NP=1,NSIM
           N=NB(NP); IF (.NOT. LDO(NP)) CYCLE

           
!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)

! at some point I need to double check this why are the
! gradients sometimes not orthogonal in OEP methods
           IF (ABS(ORTH(NP))>1E-1) THEN
              CALL vtutor%bug("LRF_HAMIL internal error: the vector H(1)-e(1) S(1) |phi(0)> is not &
                 &orthogonal to |phi(0)> " // str(ORTH(NP)) // " " // str(W0_1(NP)%CELEN), &
                 "hamil_lrf.F", 815)
           ENDIF

           IF (IDUMP==2) WRITE(*,'(I3,E11.4,"R ",2E11.4,"O ",2E14.7)') N,SQRT(ABS(FNORM(NP))),ORTH(NP),WXI%CELEN(N,NK,ISP)

           RMS=RMS+WDES%RSPIN*WDES%WTKPT(NK)*W0%FERWE(N,NK,ISP)*SQRT(ABS(FNORM(NP)))/WDES%NB_TOT
        ENDDO

        

!=======================================================================
      ENDDO bands
      ENDDO kpoints
      ENDDO spin
!=======================================================================
      CALL M_sum_d(WDES%COMM_INTER, RMS, 1)
      CALL M_sum_d(WDES%COMM_KINTER, RMS, 1)

      CALL M_sum_i(WDES%COMM_INTER, ICOUEV ,1)
      CALL M_sum_i(WDES%COMM_KINTER, ICOUEV ,1)

      DO NP=1,NSIM
         CALL DELWAV(WXI0(NP) ,.FALSE.)
         CALL DELWAV_R(W0_1(NP))
      ENDDO

!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
!! ACC_ASYNC_Q=ACC_ASYNC_ASYNC

!$ACC EXIT DATA DELETE(WXI0(:),W0_1(:),CF,ORTH,FNORM) 
      DEALLOCATE(CF)
# 848

# 864


      
# 898

    END SUBROUTINE LRF_HAMIL


!************************* SUBROUTINE HAMILTMU_LRF ********************
!
! this subroutine calculates the first order change of H
! acting onto a set of wavefuntions
! the  wavefunction must be given in reciprocal space C and real
! space CR
! CH contains the result
!    |xi>  = V_sc(1) + sum_ij | p_i> D_sc(1)_ij <p_j | phi >
!                    - sum_ij | p_i> e(1) Q(0)_ij <p_j | phi >
!                    - e(1) |phi>
!                    + |xi(0)> + sum_i xi(0)_i | p_i>
!
! V(1) is the first order change of the local potential  (SV)
! D(1) is the first order change of the PAW strength     (CDIJ)
! e(1) is the first order change of the eigen energy
!
! e(1) is evaluated during the calculation of |xi>:
! e(1) = <phi| xi(0) > +  sum_i <phi |p_i > xi(0)_i +
!                  <phi|V_sc(1)|phi> +  <phi | p_i> D_sc(1)_ij <p_j | phi >
!
! NOTE: the calling routine has to subtract e(1) |phi> to obtain the
! correct vector xi
!
!***********************************************************************

    SUBROUTINE HAMILTMU_LRF( &
         WDES1,W0_1,WXI0,NONLR_S,NONL_S,GRID,LREAL,LATT_CUR,LMDIM,CDIJ,CQIJ,SV,CH,LD,NSIM,LDO)

!! USE moffload

      USE prec
      USE mpimy
      USE mgrid
      USE wave
      USE nonl_high
      USE lattice
      USE hamil
      IMPLICIT NONE

      TYPE (grid_3d)      :: GRID
      TYPE (nonlr_struct) :: NONLR_S
      TYPE (nonl_struct ) :: NONL_S
      TYPE (wavefun1)     :: W0_1(NSIM)
      TYPE (wavefun1)     :: WXI0(NSIM)
      TYPE (wavedes1)     :: WDES1
      TYPE (latt)         :: LATT_CUR

      COMPLEX(q)      :: SV(GRID%MPLWV,WDES1%NRSPINORS*WDES1%NRSPINORS) ! local potential
      COMPLEX(q)    :: CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS), &
                    CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)

      COMPLEX(q) :: CH(LD,NSIM)

      INTEGER :: LMDIM, LD, NSIM
      LOGICAL :: LREAL
      LOGICAL :: LDO(NSIM)
! local variables
      COMPLEX(q), ALLOCATABLE :: CWORK1(:,:)
      COMPLEX(q) :: CTMP(NSIM), EVALUE1(NSIM), CE
      REAL(q) :: RINPLW, WEIGHT
      INTEGER :: NP, NGVECTOR, ISPINOR, ISPINOR_, M, MM, MM_ 

      

      ALLOCATE(CWORK1(GRID%MPLWV*WDES1%NRSPINORS,NSIM)) 

      RINPLW=1._q/GRID%NPLWV
      NGVECTOR=WDES1%NGVECTOR

!$ACC ENTER DATA CREATE(CWORK1,CTMP,EVALUE1) IF(OFFLOAD_ON)

!=======================================================================
! calculate the local contribution (result in CWORK1)
!=======================================================================
      DO NP=1,NSIM
         IF ( LDO(NP) ) THEN
# 981

            W0_1(NP)%CELEN=0
            CWORK1(:,NP)=0

            DO ISPINOR =0,WDES1%NRSPINORS-1
               CE=0
!$OMP PARALLEL DO COLLAPSE(2) SHARED(WDES1,GRID,ISPINOR,CWORK1,NP,SV,W0_1,RINPLW) PRIVATE(ISPINOR_,M,MM,MM_) REDUCTION(+:CE)
               DO ISPINOR_=0,WDES1%NRSPINORS-1
                  DO M=1,GRID%RL%NP
                     MM =M+ISPINOR *GRID%MPLWV
                     MM_=M+ISPINOR_*GRID%MPLWV
                     CWORK1(MM,NP)=  CWORK1(MM,NP)+(SV(M,1+ISPINOR_+2*ISPINOR) *W0_1(NP)%CR(MM_)*RINPLW)
                     CE=CE + CONJG(W0_1(NP)%CR(MM)) *(SV(M,1+ISPINOR_+2*ISPINOR) *W0_1(NP)%CR(MM_)*RINPLW)
                  ENDDO
               ENDDO
!$OMP END PARALLEL DO
               W0_1(NP)%CELEN=W0_1(NP)%CELEN+CE
            ENDDO
            CALL M_sum_z(WDES1%COMM_INB, W0_1(NP)%CELEN, 1)

         ENDIF
      ENDDO

!=======================================================================
! add term xi_0
!=======================================================================
      DO NP=1,NSIM
         IF ( LDO(NP) ) THEN
# 1012

            CE=0
!$OMP PARALLEL DO COLLAPSE(2) SHARED(WDES1,NGVECTOR,NP,CH,WXI0,W0_1) PRIVATE(ISPINOR,M,MM) REDUCTION(+:CE)
            DO ISPINOR =0,WDES1%NRSPINORS-1
               DO M=1,NGVECTOR
                  MM=M+ISPINOR*NGVECTOR
                  CH(MM,NP)=WXI0(NP)%CPTWFP(MM)
                  CE=CE+CONJG(W0_1(NP)%CPTWFP(MM))*CH(MM,NP)
               ENDDO
            ENDDO
!$OMP END PARALLEL DO
            CTMP(NP)=CE

            CE=0
!$OMP PARALLEL DO SHARED(WDES1,WXI0,W0_1) PRIVATE(M) REDUCTION(+:CE)
            DO M=1,WDES1%NPRO
               CE=CE+CONJG(W0_1(NP)%CPROJ(M))*WXI0(NP)%CPROJ(M)
            ENDDO
!$OMP END PARALLEL DO
            CTMP(NP)=CTMP(NP)+CE

            CALL M_sum_z(WDES1%COMM_INB, CTMP(NP), 1)

!$ACC KERNELS PRESENT(CTMP,W0_1) 
            W0_1(NP)%CELEN=W0_1(NP)%CELEN+CTMP(NP)
!$ACC END KERNELS
         ENDIF
      ENDDO

!=======================================================================
! non-local contributions to e(1)
!=======================================================================
         DO NP=1,NSIM
            IF ( LDO(NP) ) THEN

               

               CALL ECCP_NL_ALL(WDES1,W0_1(NP),W0_1(NP),CDIJ,CQIJ,0.0_q,CTMP(NP))
!$ACC KERNELS PRESENT(CTMP,W0_1,EVALUE1) 
               W0_1(NP)%CELEN=W0_1(NP)%CELEN+CTMP(NP)
               EVALUE1(NP)=W0_1(NP)%CELEN
!$ACC END KERNELS
!$ACC UPDATE SELF(W0_1(NP)%CELEN) 
            ENDIF
         ENDDO

!=======================================================================
! calculate non-local contributions to |xi> in real-space
!=======================================================================
      IF (LREAL) THEN
! contribution | p_i > D_sc(1)_ij - e(1) Q_ij < p_j |

         CALL RACCMU_C(NONLR_S,WDES1,W0_1,LMDIM,CDIJ,CQIJ,EVALUE1,CWORK1,GRID%MPLWV*WDES1%NRSPINORS,NSIM,LDO)

! contribution | p_i > cxi_i
         DO NP=1,NSIM
            IF ( LDO(NP) ) THEN
# 1075

                  CALL RACC0(NONLR_S,WDES1,WXI0(NP)%CPROJ(1),CWORK1(1,NP))
# 1079

               DO ISPINOR=0,WDES1%NRSPINORS-1
                  CALL FFTEXT(NGVECTOR,WDES1%NINDPW(1),CWORK1(1+ISPINOR*WDES1%GRID%MPLWV,NP),CH(1+ISPINOR*NGVECTOR,NP),GRID,.TRUE.)
               ENDDO
            ENDIF
         ENDDO

!=======================================================================
! calculate non-local contribution to |xi> in reciprocal space
!=======================================================================
      ELSE
         DO NP=1,NSIM
            IF ( LDO(NP) ) THEN
! contribution | p_i > cxi_i
# 1099

                   CALL VNLAC0(NONL_S,WDES1,WXI0(NP)%CPROJ(1),CH(1,NP))
# 1103

! contribution | p_i > D_(sc)_ij - e(1) Q_ij < p_j |

               CALL VNLACC_ADD_C(NONL_S,W0_1(NP),CDIJ,CQIJ,1,EVALUE1(NP),CH(:,NP))

               DO ISPINOR=0,WDES1%NRSPINORS-1
                  CALL FFTEXT(NGVECTOR,WDES1%NINDPW(1),CWORK1(1+ISPINOR*WDES1%GRID%MPLWV,NP),CH(1+ISPINOR*NGVECTOR,NP),GRID,.TRUE.)
               ENDDO
            ENDIF
         ENDDO
      ENDIF

# 1119


      DEALLOCATE(CWORK1)

      
# 1175

    END SUBROUTINE HAMILTMU_LRF

  END MODULE hamil_lrf
