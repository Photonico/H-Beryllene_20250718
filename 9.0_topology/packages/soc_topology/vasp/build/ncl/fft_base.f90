# 1 "fft_base.F"
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


# 2 "fft_base.F" 2 


!************************* SUBROUTINE FFTCHK ***************************
!
!> Returns the next correct setting for the three dimensional FFT
!
!***********************************************************************

      SUBROUTINE FFTCHK(NFFT)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION NFFT(3)
      LOGICAL FFTCH1

      DO 100 IND=1,3
  200 CONTINUE
        IF (FFTCH1(NFFT(IND))) GOTO 100
        NFFT(IND)=NFFT(IND)+1
        GOTO 200
  100 CONTINUE
      END SUBROUTINE FFTCHK


      LOGICAL FUNCTION FFTCH1(NIN)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      PARAMETER (NFACT=4)
      DIMENSION IFACT(NFACT),NCOUNT(NFACT)
      DATA      IFACT /2,3,5,7/
      N=NIN
      DO 100 I=1,NFACT
         NCOUNT(I)=0
120      NEXT=N/IFACT(I)
         IF (NEXT*IFACT(I)==N) THEN
            N=NEXT
            NCOUNT(I)=NCOUNT(I)+1
            GOTO 120
         ENDIF
100   ENDDO
      IF (N==1 .AND. (NCOUNT(1)/=0)) &
           &  THEN
         FFTCH1=.TRUE.
      ELSE
         FFTCH1=.FALSE.
      ENDIF
      RETURN
      END FUNCTION FFTCH1


!************************* SUBROUTINE FFTGRIDPLAN **********************
!
!***********************************************************************

    SUBROUTINE FFTGRIDPLAN(GRID)
      USE prec
      USE mgrid_struct_def

      IMPLICIT NONE

      TYPE (grid_3d) GRID

! local variables
      COMPLEX(q), ALLOCATABLE :: CWORK(:)
      INTEGER, PARAMETER :: NADD=1024, NTIMES=4
      INTEGER :: I

      IF (GRID%RL%NFAST/=1 .OR. GRID%RL_FFT%NFAST/=1) THEN
         CALL FFTGRIDPLAN_MPI(GRID)
         RETURN
      ENDIF

      ALLOCATE(CWORK(GRID%MPLWV+NADD))

      DO I=1,NTIMES
         CALL INIDAT(GRID%MPLWV+NADD,CWORK)
         CALL FFTMAKEPLAN(CWORK(I),GRID)
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP MASTER
!$       CALL FFTMAKEPLAN(CWORK(I),GRID)
!$OMP END MASTER
!$OMP END PARALLEL
      ENDDO

      DEALLOCATE(CWORK)

      RETURN
    END SUBROUTINE FFTGRIDPLAN


!************************* SUBROUTINE FFTINI ***************************
!
!>  If necessary this routine performes initialization
!>  for #FFTWAV and #FFTEXT
!>
!>  Usually this is only necessary for the Gamma-only
!>  single k-point version
!>
!>   FFTSCA(.,1) is the scaling factor for extracting the wavefunction
!>               from the FFT grid (#FFTEXT)
!>
!>   FFTSCA(.,2) is the scaling factor for puting the wavefunction on
!>               the grid
!
!***********************************************************************

    SUBROUTINE FFTINI(NINDPW,NPLWKP,NKPTS,NRPLW,GRID)
      USE prec
      USE mgrid_struct_def
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE (grid_3d)  GRID
      INTEGER :: NRPLW               !< maximum number of plane wave coefficients
      INTEGER :: NKPTS               !< number of k-points
      INTEGER :: NINDPW(NRPLW,NKPTS) !< index array from coefficients in sphere to 3D grid
      INTEGER :: NPLWKP(NKPTS)       !< actual number of plane wave coefficients for each k-points
! local
      INTEGER :: NK, NX, NY, NZ, NPL, N, IND, N1, N2, N3, NC
      INTEGER :: N2INV, N3INV, NUMBER_OF_CONJG, NCP
      REAL(q) :: FACTM

      IF (GRID%REAL2CPLX) THEN

         IF (GRID%RL%NFAST/=1) THEN
            CALL FFTINI_MPI(NINDPW,NPLWKP,NKPTS,NRPLW,GRID)
            RETURN
         ENDIF

         IF (NKPTS>1) THEN
            CALL vtutor%error("FFT3D: real version works only for 1 k-point")
         ENDIF

         NK=1
         NX=GRID%NGPTAR(1)
         NY=GRID%NGPTAR(2)
         NZ=GRID%NGPTAR(3)
         NPL=NPLWKP(NK)
         NULLIFY(GRID%FFTSCA)
         ALLOCATE(GRID%FFTSCA(NPL,2))

         NUMBER_OF_CONJG=0
         DO N=1,NPL
            IND=NINDPW(N,NKPTS)
            N1= MOD((IND-1),GRID%RC%NROW)+1
            NC= (IND-1)/GRID%RC%NROW+1
            N2= GRID%RC%I2(NC)
            N3= GRID%RC%I3(NC)
            IF (N1==1 .AND. N/=1) THEN
               NUMBER_OF_CONJG=NUMBER_OF_CONJG+1
! invert second and third index
               N2INV=MOD(-GRID%LPCTY(N2)+GRID%NGY,GRID%NGY)+1
               N3INV=MOD(-GRID%LPCTZ(N3)+GRID%NGZ,GRID%NGZ)+1
! now determine the corresponding column index
               NCP=N2INV+(N3INV-1)*GRID%NGY
! check whether correct
               IF ( GRID%LPCTY(N2)+GRID%LPCTY(GRID%RC%I2(NCP)) /=0 .OR.  GRID%LPCTZ(N3)+GRID%LPCTZ(GRID%RC%I3(NCP))/=0) THEN
                  CALL vtutor%bug("internal error in FFTINI: could not determine conjugated &
                     &coefficient " // str(GRID%LPCTY(N2)) // " " // str(GRID%LPCTY(GRID%RC%I2(NCP))) &
                     // " " // str(GRID%LPCTZ(N3) + GRID%LPCTZ(GRID%RC%I3(NCP))), "fft_base.F", 162)
               ENDIF
            ENDIF
         ENDDO
         NULLIFY(GRID%NINDPWCONJG, GRID%IND_IN_SPHERE)
! allocate the required index arrays
         ALLOCATE(GRID%NINDPWCONJG(NUMBER_OF_CONJG), GRID%IND_IN_SPHERE(NUMBER_OF_CONJG))

         NUMBER_OF_CONJG=0
         DO N=1,NPL
            IND=NINDPW(N,NKPTS)
            N1= MOD((IND-1),GRID%RC%NROW)+1
            NC= (IND-1)/GRID%RC%NROW+1
            N2= GRID%RC%I2(NC)
            N3= GRID%RC%I3(NC)
            FACTM=SQRT(2._q)
            IF (N==1) FACTM=1
! in the real version the coefficients stored in the compressed mode (sphere)
! are multiplied by a factor sqrt(2) compared to complex version (except at Gamma)
! this allows to use DGEMM calls to calculate inproducts
            GRID%FFTSCA(N,1)= FACTM
            GRID%FFTSCA(N,2)= 1/FACTM

            IF (N1==1 .AND. N/=1) THEN
               NUMBER_OF_CONJG=NUMBER_OF_CONJG+1
! invert second and third index
               N2INV=MOD(-GRID%LPCTY(N2)+GRID%NGY,GRID%NGY)+1
               N3INV=MOD(-GRID%LPCTZ(N3)+GRID%NGZ,GRID%NGZ)+1
! now determine the corresponding column index
               NCP=N2INV+(N3INV-1)*GRID%NGY
! store original index in the compressed storage mode (sphere)
               GRID%IND_IN_SPHERE(NUMBER_OF_CONJG)=N
! store index to the 3D grid
               GRID%NINDPWCONJG(NUMBER_OF_CONJG)=N1+(NCP-1)*GRID%RC%NROW
            ENDIF
# 204

         ENDDO
      ENDIF

    END SUBROUTINE FFTINI


!************************* SUBROUTINE FFTWAV ***************************
!
!> Transforms a wavefunction C defined within the cutoff-sphere
!> to real space CR
!>
!> @note For the real version (gamma point only) it is assumed
!> that the wavefunctions at NGX != 0 (wNGXhalf)
!> are multiplied by a factor sqrt(2) on the reduced plane wave grid
!> this factor has to be removed before the FFT transformation
!> (scaling with FFTSCA(M,2))!
!
!***********************************************************************

    SUBROUTINE FFTWAV(NPL,NINDPW,CR,C,GRID)

!! USE moffload_struct_def

      USE prec
      USE mgrid_struct_def
      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      COMPLEX(q) :: C(NPL),CR(GRID%NPLWV)
      INTEGER    :: NINDPW(NPL)
      INTEGER    :: NPL
! local variables
      INTEGER :: M

!! LOGICAL :: ACC_ACTIVE

      

      IF (GRID%RL%NFAST/=1) THEN
         CALL FFTWAV_MPI(NPL,NINDPW,CR,C,GRID)
         
         RETURN
      ENDIF

!! ACC_ACTIVE=ACC_IS_PRESENT(NINDPW).AND.ACC_IS_PRESENT(CR).AND. &
!!     ACC_IS_PRESENT(C).AND.ACC_IS_PRESENT(GRID,1).AND.OFFLOAD_ON

! (0._q,0._q) all elements on the grid
!$ACC PARALLEL LOOP PRESENT(CR,GRID) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
!DIR$ IVDEP
!OCL NOVREC
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(CR,GRID) PRIVATE(M)
      DO M=1,GRID%NGX_rd*GRID%NGY_rd*GRID%NGZ_rd
         CR(M)=(0.0_q,0.0_q)
      ENDDO
 !$OMP END PARALLEL DO
! now fill in non (0._q,0._q) elements from
! within the radial cutoff sphere
      IF (GRID%REAL2CPLX) THEN
!$ACC PARALLEL LOOP PRESENT(CR,C,NINDPW,GRID) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
!DIR$ IVDEP
!OCL NOVREC
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,CR,NINDPW,GRID,C) PRIVATE(M)
         DO M=1,NPL
            CR(NINDPW(M))=C(M)*GRID%FFTSCA(M,2)
         ENDDO
 !$OMP END PARALLEL DO

!$ACC PARALLEL LOOP PRESENT(CR,C,GRID) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
!DIR$ IVDEP
!OCL NOVREC
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(CR,GRID,C) PRIVATE(M)
         DO M=1,SIZE(GRID%IND_IN_SPHERE)
            CR(GRID%NINDPWCONJG(M))=CONJG(C(GRID%IND_IN_SPHERE(M)))*GRID%FFTSCA(GRID%IND_IN_SPHERE(M),2)
         ENDDO
 !$OMP END PARALLEL DO

      ELSE
!$ACC PARALLEL LOOP PRESENT(CR,C,NINDPW) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
!DIR$ IVDEP
!OCL NOVREC
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,CR,NINDPW,C) PRIVATE(M)
         DO M=1,NPL
            CR(NINDPW(M))=C(M)
         ENDDO
 !$OMP END PARALLEL DO
      ENDIF

      CALL FFT3D(CR,GRID,1)

      

      RETURN
    END SUBROUTINE FFTWAV


!************************* SUBROUTINE FFTWAV_USEINV ********************
!
!> Transforms a wavefunction C defined within the cutoff-sphere
!> to real space CR
!>
!> @note This routine differs from #FFTWAV in that it allows for the input
!> data to be only defined on half the grid,
!> and uses an inversion through the mid point to determine
!> the PW coefficients at -G-k
!
!***********************************************************************

    SUBROUTINE FFTWAV_USEINV(NPL,NINDPW,NINDPW_INV,FFTSCA,CR,C,GRID)

!! USE moffload_struct_def

      USE prec
      USE mgrid_struct_def
      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      COMPLEX(q) :: C(NPL),CR(GRID%NPLWV)
      REAL(q)    :: FFTSCA(NPL)
      INTEGER    :: NINDPW(NPL),NINDPW_INV(NPL)
      INTEGER    :: NPL
! local variables
      INTEGER :: M

!! LOGICAL :: ACC_ACTIVE

      

!! ACC_ACTIVE=ACC_IS_PRESENT(NINDPW).AND.ACC_IS_PRESENT(NINDPW_INV).AND.ACC_IS_PRESENT(FFTSCA).AND. &
!!     ACC_IS_PRESENT(CR).AND.ACC_IS_PRESENT(C).AND.ACC_IS_PRESENT(GRID,1).AND.OFFLOAD_ON

!$ACC PARALLEL LOOP PRESENT(CR,GRID) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
!DIR$ IVDEP
!OCL NOVREC
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(CR,GRID) PRIVATE(M)
      DO M=1,GRID%NGX_rd*GRID%NGY_rd*GRID%NGZ_rd
         CR(M)=(0.0_q,0.0_q)
      ENDDO
 !$OMP END PARALLEL DO
! now fill in non (0._q,0._q) elements from
! within the radial cutoff sphere
      IF (GRID%REAL2CPLX) THEN
!$ACC PARALLEL LOOP PRESENT(CR,NINDPW,C,GRID) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
!DIR$ IVDEP
!OCL NOVREC
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,CR,C,NINDPW,GRID) PRIVATE(M)
         DO M=1,NPL
            CR(NINDPW(M))=C(M)*GRID%FFTSCA(M,2)
         ENDDO
 !$OMP END PARALLEL DO

!$ACC PARALLEL LOOP PRESENT(CR,C,GRID) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
!DIR$ IVDEP
!OCL NOVREC
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(CR,GRID,C) PRIVATE(M)
         DO M=1,SIZE(GRID%IND_IN_SPHERE)
            CR(GRID%NINDPWCONJG(M))=CONJG(C(GRID%IND_IN_SPHERE(M)))*GRID%FFTSCA(GRID%IND_IN_SPHERE(M),2)
         ENDDO
 !$OMP END PARALLEL DO

      ELSE
!$ACC PARALLEL LOOP PRESENT(CR,NINDPW_INV,C,FFTSCA) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
!DIR$ IVDEP
!OCL NOVREC
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,CR,C,NINDPW_INV,FFTSCA) PRIVATE(M)
         DO M=1,NPL
            CR(NINDPW_INV(M))=CONJG(C(M))*FFTSCA(M)
         ENDDO
 !$OMP END PARALLEL DO
!$ACC PARALLEL LOOP PRESENT(CR,C,NINDPW,FFTSCA) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
!DIR$ IVDEP
!OCL NOVREC
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,CR,NINDPW,C,FFTSCA) PRIVATE(M)
         DO M=1,NPL
            CR(NINDPW(M))=C(M)*FFTSCA(M)
         ENDDO
 !$OMP END PARALLEL DO
      ENDIF

      CALL FFT3D(CR,GRID,1)

      

      RETURN
    END SUBROUTINE FFTWAV_USEINV


!************************* SUBROUTINE FFTEXT ***************************
!
!> Performes a FFT to reciprocal space and extracts data from the FFT-mesh
!>
!> @note For the real version (gamma point only) it is assumed
!> that the wavefunctions at NGX != 0
!> are multiplied by a factor sqrt(2) on the reduced grid
!> this factor has to be applied after the FFT transformation
!> (scaling with FFTSCA(M))!
!
!***********************************************************************

    SUBROUTINE FFTEXT(NPL,NINDPW,CR,C,GRID,LADD)

!! USE moffload_struct_def

      USE prec
      USE mgrid_struct_def
      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      COMPLEX(q) :: C(NPL),CR(GRID%NPLWV)
      INTEGER    :: NINDPW(NPL)
      INTEGER    :: NPL
      LOGICAL    :: LADD
! local variables
      INTEGER :: M

!! LOGICAL :: ACC_ACTIVE

      

      CALL FFT3D(CR,GRID,-1)

!! ACC_ACTIVE=ACC_IS_PRESENT(NINDPW).AND.ACC_IS_PRESENT(CR).AND. &
!!     ACC_IS_PRESENT(C).AND.ACC_IS_PRESENT(GRID,1).AND.OFFLOAD_ON

!gK LUSEINV TODO: REAL2CPLX could possibly go
      IF (LADD .AND. GRID%REAL2CPLX) THEN
!$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW,GRID) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
!DIR$ IVDEP
!OCL NOVREC
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,C,CR,NINDPW,GRID) PRIVATE(M)
         DO M=1,NPL
            C(M)=C(M)+CR(NINDPW(M))*GRID%FFTSCA(M,1)
         ENDDO
 !$OMP END PARALLEL DO
      ELSE IF (LADD .AND. .NOT. GRID%REAL2CPLX) THEN
!$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
!DIR$ IVDEP
!OCL NOVREC
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,C,CR,NINDPW) PRIVATE(M)
         DO M=1,NPL
            C(M)=C(M)+CR(NINDPW(M))
         ENDDO
 !$OMP END PARALLEL DO
      ELSE IF (GRID%REAL2CPLX) THEN
!$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW,GRID) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
!DIR$ IVDEP
!OCL NOVREC
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,C,CR,NINDPW,GRID) PRIVATE(M)
        DO M=1,NPL
          C(M)=CR(NINDPW(M))*GRID%FFTSCA(M,1)
        ENDDO
 !$OMP END PARALLEL DO
     ELSE
!$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
!DIR$ IVDEP
!OCL NOVREC
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,C,CR,NINDPW) PRIVATE(M)
        DO M=1,NPL
          C(M)=CR(NINDPW(M))
        ENDDO
 !$OMP END PARALLEL DO
      ENDIF

      

      RETURN
    END SUBROUTINE FFTEXT


!************************* SUBROUTINE FFTEXT_USEINV ********************
!
!> Performes a FFT to reciprocal space and extracts data from the FFT-mesh
!>
!> @note For the real version (gamma point only) it is assumed
!> that the wavefunctions at NGX != 0
!> are multiplied by a factor sqrt(2) on the reduced grid
!> this factor has to be applied after the FFT transformation
!> (scaling with FFTSCA(M))!
!
!***********************************************************************

    SUBROUTINE FFTEXT_USEINV(NPL,NINDPW,FFTSCA,CR,C,GRID,LADD)

!! USE moffload_struct_def

      USE prec
      USE mgrid_struct_def
      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      COMPLEX(q) :: C(NPL),CR(GRID%NPLWV)
      REAL(q)    :: FFTSCA(NPL)
      INTEGER    :: NINDPW(NPL)
      INTEGER    :: NPL
      LOGICAL    :: LADD
! local variables
      INTEGER :: M
# 506

      

      CALL FFT3D(CR,GRID,-1)

!gK LUSEINV TODO: REAL2CPLX could possibly go
      IF (LADD .AND. GRID%REAL2CPLX) THEN
!$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW,GRID) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
!DIR$ IVDEP
!OCL NOVREC
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,C,CR,NINDPW,GRID) PRIVATE(M)
         DO M=1,NPL
            C(M)=C(M)+CR(NINDPW(M))*GRID%FFTSCA(M,1)
         ENDDO
 !$OMP END PARALLEL DO
      ELSE IF (LADD .AND. .NOT. GRID%REAL2CPLX) THEN
!$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW,FFTSCA) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
!DIR$ IVDEP
!OCL NOVREC
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,C,CR,NINDPW,FFTSCA) PRIVATE(M)
         DO M=1,NPL
            C(M)=C(M)+CR(NINDPW(M))*FFTSCA(M)
         ENDDO
 !$OMP END PARALLEL DO
      ELSE IF (GRID%REAL2CPLX) THEN
!$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW,GRID) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
!DIR$ IVDEP
!OCL NOVREC
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,C,CR,NINDPW,GRID) PRIVATE(M)
        DO M=1,NPL
          C(M)=CR(NINDPW(M))*GRID%FFTSCA(M,1)
        ENDDO
 !$OMP END PARALLEL DO
     ELSE
!$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW,FFTSCA) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
!DIR$ IVDEP
!OCL NOVREC
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,C,CR,NINDPW,FFTSCA) PRIVATE(M)
        DO M=1,NPL
          C(M)=CR(NINDPW(M))*FFTSCA(M)
        ENDDO
 !$OMP END PARALLEL DO
      ENDIF

      

      RETURN
    END SUBROUTINE FFTEXT_USEINV


!************************* SUBROUTINE FFT3D ****************************
!
!> 3-d fast fourier transform (possibly real to complex and vice versa)
!> for charge densities and potentials
!>
!> @param ISN = +1  q->r   vr= sum(q) vq exp(+iqr) (might be C2R)\n
!>            = -1  r->q   vq= sum(r) vr exp(-iqr) (might be R2C)
!
!***********************************************************************

    SUBROUTINE FFT3D(C,GRID,ISN)
      USE prec
      USE mgrid_struct_def
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      REAL(q) :: C(*)
      INTEGER ::ISN
! local variables
      INTEGER :: NX, NY, NZ
      INTEGER :: IL, NDEST, NSRC

      

      IF (GRID%RL%NFAST/=1 .OR. GRID%RL_FFT%NFAST/=1) THEN
         CALL FFT3D_MPI(C, GRID, ISN)
         
         RETURN
      ENDIF

      NX=GRID%NGPTAR(1)
      NY=GRID%NGPTAR(2)
      NZ=GRID%NGPTAR(3)

!-------------------------------------------------------------------------------
!  complex to complex version
!-------------------------------------------------------------------------------
      IF (.NOT. GRID%REAL2CPLX .AND. .NOT. GRID%LREAL ) THEN
         IF (.NOT. (NX==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ==GRID%NGZ_rd) ) THEN
            CALL vtutor%bug("internal error 1 in FFT3D: something not properly set " // &
               str(GRID%LREAL) // " " // str(GRID%REAL2CPLX) // "\n " // str(NX) // " " // str(NY) // &
               " " // str(NZ) // "\n " // str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " // &
               str(GRID%NGZ_rd), "fft_base.F", 600)
         ENDIF

         CALL FFTBAS_PLAN(C,GRID,ISN)
!-------------------------------------------------------------------------------
!  complex to complex version, but with a real array in real space
!-------------------------------------------------------------------------------
      ELSE IF (.NOT. GRID%REAL2CPLX .AND. GRID%LREAL) THEN
         IF (.NOT. (NX==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ==GRID%NGZ_rd) ) THEN
            CALL vtutor%bug("internal error 2 in FFT3D: something not properly set " // &
               str(GRID%LREAL) // " " // str(GRID%REAL2CPLX) // "\n " // str(NX) // " " // str(NY) // &
               " " // str(NZ) // "\n " // str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " // &
               str(GRID%NGZ_rd), "fft_base.F", 612)
         ENDIF

!     q->r FFT
         IF (ISN==1) THEN
            CALL FFTBAS_PLAN(C,GRID,ISN)

!  go from complex stride 2 to 1
!DIR$ IVDEP
!OCL NOVREC
            DO IL=0,NX*NY*NZ-1
               NDEST=IL+1
               NSRC =IL*2+1
               C(NDEST)=C(NSRC)
            ENDDO
         ELSE

!     r->q FFT
!  go from stride 1 to stride 2
!DIR$ IVDEP
!OCL NOVREC
            DO IL=NX*NY*NZ-1,0,-1
               NSRC =IL+1
               NDEST=IL*2+1
               C(NDEST)=C(NSRC)
               C(NDEST+1)=0
            ENDDO
            CALL FFTBAS_PLAN(C,GRID,ISN)
         ENDIF
!-------------------------------------------------------------------------------
!  real to complex FFT  only half grid mode in X direction supported
!  data are stored as real array in real space
!-------------------------------------------------------------------------------
      ELSE IF (GRID%LREAL) THEN
         IF (.NOT. (NX/2+1==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ==GRID%NGZ_rd) ) THEN
            CALL vtutor%bug("internal error 3 in FFT3D: something not properly set " // &
               str(GRID%LREAL) // " " // str(GRID%REAL2CPLX) // "\n " // str(NX) // " " // str(NY) // &
               " " // str(NZ) // "\n " // str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " // &
               str(GRID%NGZ_rd), "fft_base.F", 650)
         ENDIF

!  in real space the first dimension in VASP is NGX (REAL data)
!  but the FFT required NGX+2 (real data)
!  therefore some data movement is required

         IF (ISN==1) THEN
!     q->r FFT
            CALL FFTBRC_PLAN(C,GRID,ISN)
!  concat  x-lines (go from stride NX+2 to NX)
            CALL RESTRIDE_Q2R(NY*NZ,NX,C)
         ELSE
!     x-lines (go from stride NX to NX+2)
            CALL RESTRIDE_R2Q(NY*NZ,NX,C)
!     r->q FFT
            CALL FFTBRC_PLAN(C,GRID,ISN)
         ENDIF
!-------------------------------------------------------------------------------
! same as above (real to complex FFT) but this time the data layout
! is complex in real space
!-------------------------------------------------------------------------------
      ELSE
         IF (.NOT. (NX/2+1==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ==GRID%NGZ_rd) ) THEN
            CALL vtutor%bug("internal error 4 in FFT3D: something not properly set " // &
               str(GRID%LREAL) // " " // str(GRID%REAL2CPLX) // "\n " // str(NX) // " " // str(NY) // &
               " " // str(NZ) // "\n " // str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " // &
               str(GRID%NGZ_rd), "fft_base.F", 677)
         ENDIF

         IF (ISN==1) THEN
!     q->r FFT
            CALL FFTBRC_PLAN(C,GRID,ISN)
! concat  x-lines (go from "real" stride NX+2 to complex stride NX)
            CALL RESTRIDE_Q2R_CMPLX(NY*NZ,NX,C)
         ELSE
!     x-lines (go from complex stride NX to real stride NX+2)
            CALL RESTRIDE_R2Q_CMPLX(NY*NZ,NX,C)
!     r->q FFT
            CALL FFTBRC_PLAN(C,GRID,ISN)
         ENDIF
      ENDIF

      

      RETURN
    END SUBROUTINE FFT3D


!************************* SUBROUTINE EXTBAS ***************************
!
!> Extracts data from the FFT-mesh
!>
!> @note For the real version (gamma point only) it is assumed
!> that the wavefunctions at NGX != 0
!> are multiplied by a factor sqrt(2) on the reduced grid
!> this factor has to be applied after the FFT transformation
!> (scaling with FFTSCA(M))!
!
!***********************************************************************

    SUBROUTINE EXTBAS(NPL,NINDPW,CR,C,GRID,LADD)

!! USE moffload_struct_def

      USE prec
      USE mgrid_struct_def
      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      COMPLEX(q) :: C(NPL),CR(GRID%MPLWV)
      INTEGER    :: NINDPW(NPL)
      INTEGER    :: NPL
      LOGICAL    :: LADD
! local variables
      INTEGER :: M
# 730

      

!gK LUSEINV TODO: REAL2CPLX could possibly go
      IF (LADD .AND. GRID%REAL2CPLX) THEN
!$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW,GRID) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
!DIR$ IVDEP
!OCL NOVREC
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,C,CR,NINDPW,GRID) PRIVATE(M)
         DO M=1,NPL
            C(M)=C(M)+CR(NINDPW(M))*GRID%FFTSCA(M,1)
         ENDDO
 !$OMP END PARALLEL DO
      ELSE IF (LADD .AND. .NOT. GRID%REAL2CPLX) THEN
!$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
!DIR$ IVDEP
!OCL NOVREC
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,C,CR,NINDPW) PRIVATE(M)
         DO M=1,NPL
            C(M)=C(M)+CR(NINDPW(M))
         ENDDO
 !$OMP END PARALLEL DO
      ELSE IF (GRID%REAL2CPLX) THEN
!$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW,GRID) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
!DIR$ IVDEP
!OCL NOVREC
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,C,CR,NINDPW,GRID) PRIVATE(M)
        DO M=1,NPL
          C(M)=CR(NINDPW(M))*GRID%FFTSCA(M,1)
        ENDDO
 !$OMP END PARALLEL DO
     ELSE
!$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
!DIR$ IVDEP
!OCL NOVREC
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NPL,C,CR,NINDPW) PRIVATE(M)
        DO M=1,NPL
          C(M)=CR(NINDPW(M))
        ENDDO
 !$OMP END PARALLEL DO
      ENDIF

      

      RETURN
    END SUBROUTINE EXTBAS


!************************* SUBROUTINE MULZ *****************************
!
!> Multiplies the Z!=0 components by a factor FACT
!
!***********************************************************************

    SUBROUTINE MULZ(C,NGX,NGY,NGZ,FACT)
      USE prec
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      COMPLEX(q) C(0:NGX-1,0:NGY-1,0:NGZ-1)

      DO N3=1,NGZ/2-1
       DO N2=0,NGY-1
        DO N1=0,NGX-1
          C(N1,N2,N3)= C(N1,N2,N3)*FACT
      ENDDO
      ENDDO
      ENDDO

      RETURN
    END SUBROUTINE MULZ


    SUBROUTINE RESTRIDE_Q2R(NCOL,N,C)

!! USE moffload_struct_def

      USE prec
      USE mgrid_struct_def
      INTEGER :: NCOL,N
      REAL(q) :: C(*)
! local
      INTEGER :: ICOL,NDEST,NSRC,I
# 825

!DIR$ IVDEP
!OCL NOVREC
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(C,CTMP) PRIVATE(NSRC,NDEST) 
      DO ICOL=1,NCOL-1
    NDEST=ICOL* N
    NSRC =ICOL*(N+2)
!DIR$ IVDEP
!OCL NOVREC
         DO I=1,N
!!       NDEST=ICOL* N
!!       NSRC =ICOL*(N+2)
       C(NDEST+I)=   C(NSRC+I)
!!       C(NDEST+I)=CTMP(NSRC+I)
         ENDDO
      ENDDO
# 845

    END SUBROUTINE RESTRIDE_Q2R


    SUBROUTINE RESTRIDE_R2Q(NCOL,N,C)

!! USE moffload_struct_def

      USE prec
      USE mgrid_struct_def
      INTEGER :: NCOL,N
      REAL(q) :: C(*)
! local
      INTEGER :: ICOL,NDEST,NSRC,I
# 872

!DIR$ IVDEP
!OCL NOVREC
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(C,CTMP) PRIVATE(NSRC,NDEST) 
      DO ICOL=NCOL-1,1,-1
    NSRC =ICOL* N
    NDEST=ICOL*(N+2)
! ifc10.1 has trouble vectorizing this statement
!!DIR$ IVDEP
!!OCL NOVREC
         DO I=N,1,-1
!!       NSRC =ICOL* N
!!       NDEST=ICOL*(N+2)
       C(NDEST+I)=   C(NSRC+I)
!!       C(NDEST+I)=CTMP(NSRC+I)
         ENDDO
      ENDDO
# 893

    END SUBROUTINE RESTRIDE_R2Q


    SUBROUTINE RESTRIDE_Q2R_CMPLX(NCOL,N,C)

!! USE moffload_struct_def

      USE prec
      USE mgrid_struct_def
      INTEGER :: NCOL,N
      REAL(q) :: C(*)
! local
      INTEGER :: ICOL,NDEST,NSRC,I
# 920

!DIR$ IVDEP
!OCL NOVREC
!$ACC PARALLEL LOOP COLLAPSE(2) GANG PRESENT(C,CTMP) PRIVATE(NDEST,NSRC) 
      DO ICOL=NCOL-1,0,-1
    NDEST=ICOL* N*2
    NSRC =ICOL*(N+2)
!DIR$ IVDEP
!OCL NOVREC
         DO I=N,1,-1
!!       NDEST=ICOL* N*2
!!       NSRC =ICOL*(N+2)
       C(NDEST+I*2-1)=C(NSRC+I)
!!       C(NDEST+I*2-1)=CTMP(NSRC+I)
            C(NDEST+I*2  )=0
         ENDDO
      ENDDO
# 941

    END SUBROUTINE RESTRIDE_Q2R_CMPLX


    SUBROUTINE RESTRIDE_R2Q_CMPLX(NCOL,N,C)

!! USE moffload_struct_def

      USE prec
      USE mgrid_struct_def
      USE tutor, ONLY: vtutor
      INTEGER :: NCOL,N
      REAL(q) :: C(*)
! local
      INTEGER :: ICOL,NDEST,NSRC,I
# 969

      IF (NCOL<1) RETURN
      IF (N<2) CALL vtutor%bug("RESTRIDE_R2Q_CMPLX: can not restride array for N=1", "fft_base.F", 971)

!DIR$ IVDEP
!OCL NOVREC
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(C,CTMP) PRIVATE(NSRC,NDEST) 
      DO ICOL=0,NCOL-1
    NSRC =ICOL* N*2
    NDEST=ICOL*(N+2)
!DIR$ IVDEP
!OCL NOVREC
         DO I=1,N
!!       NSRC =ICOL* N*2
!!       NDEST=ICOL*(N+2)
       C(NDEST+I)=   C(NSRC+I*2-1)
!!       C(NDEST+I)=CTMP(NSRC+I*2-1)
         ENDDO
      ENDDO
# 992

    END SUBROUTINE RESTRIDE_R2Q_CMPLX


!************************* SUBROUTINE FFT1D ****************************
!
!> One dimensional FFT
!> Currently only used in stm.F (which is pretty irrelevant)
!
!***********************************************************************

    SUBROUTINE FFT1D(C,NFFT,ISIGN)
      USE prec

      IMPLICIT NONE

      COMPLEX(q) :: C(*)
      INTEGER    :: NFFT
      INTEGER    :: ISIGN

      CALL FFT1D_C2C_PLAN(C,NFFT,ISIGN)

    END SUBROUTINE FFT1D


!************************* SUBROUTINE FFTWAV_MU ************************
!
!***********************************************************************

    SUBROUTINE FFTWAV_MU(NPL,N,NINDPW,CR,LDR,C,LDC,GRID)

!! USE moffload_struct_def

      USE prec
      USE mgrid_struct_def

      IMPLICIT NONE

      TYPE(grid_3d) :: GRID

      COMPLEX(q) :: CR(LDR,*),C(LDC,*)
      INTEGER    :: NINDPW(NPL)
      INTEGER    :: NPL, N, LDR, LDC

! local variables
      INTEGER :: I,J

# 1043


      

! (0._q,0._q) all elements on the grid
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CR,GRID) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
 !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,GRID,CR) PRIVATE(I,J)
      DO I=1,N
         DO J=1,GRID%NGX_rd*GRID%NGY_rd*GRID%NGZ_rd
            CR(J,I)=(0.0_q,0.0_q)
         ENDDO
      ENDDO
 !$OMP END PARALLEL DO

! now fill in non (0._q,0._q) elements from
! within the radial cutoff sphere
      IF (GRID%REAL2CPLX) THEN

!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CR,C,NINDPW,GRID) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
 !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,CR,NINDPW,C,GRID) PRIVATE(I,J)
         DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
            DO J=1,NPL
               CR(NINDPW(J),I)=C(J,I)*GRID%FFTSCA(J,2)
            ENDDO
         ENDDO
 !$OMP END PARALLEL DO

!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CR,C,GRID) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
 !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,CR,C,GRID) PRIVATE(I,J)
         DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
            DO J=1,SIZE(GRID%IND_IN_SPHERE)
               CR(GRID%NINDPWCONJG(J),I)=CONJG(C(GRID%IND_IN_SPHERE(J),I))*GRID%FFTSCA(GRID%IND_IN_SPHERE(J),2)
            ENDDO
         ENDDO
 !$OMP END PARALLEL DO


      ELSE

!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CR,C,NINDPW) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
 !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,CR,NINDPW,C) PRIVATE(I,J)
         DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
            DO J=1,NPL
               CR(NINDPW(J),I)=C(J,I)
            ENDDO
         ENDDO
 !$OMP END PARALLEL DO

      ENDIF

      CALL FFT3D_MU(N,CR,LDR,GRID,1)

      

    END SUBROUTINE FFTWAV_MU


!************************* SUBROUTINE FFTWAV_USEINV_MU *****************
!
!***********************************************************************

    SUBROUTINE FFTWAV_USEINV_MU(NPL,N,NINDPW,NINDPW_INV,FFTSCA,CR,LDR,C,LDC,GRID)

!! USE moffload_struct_def

      USE prec
      USE mgrid_struct_def

      IMPLICIT NONE

      TYPE(grid_3d) :: GRID

      COMPLEX(q) :: CR(LDR,*),C(LDC,*)
      REAL(q)    :: FFTSCA(NPL)
      INTEGER    :: NINDPW(NPL),NINDPW_INV(NPL)
      INTEGER    :: NPL, N, LDR, LDC

! local variables
      INTEGER :: I,J

# 1133


      

! (0._q,0._q) all elements on the grid
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CR,GRID) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
 !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,GRID,CR) PRIVATE(I,J)
      DO I=1,N
         DO J=1,GRID%NGX_rd*GRID%NGY_rd*GRID%NGZ_rd
            CR(J,I)=(0.0_q,0.0_q)
         ENDDO
      ENDDO
 !$OMP END PARALLEL DO

! now fill in non (0._q,0._q) elements from
! within the radial cutoff sphere
      IF (GRID%REAL2CPLX) THEN

!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CR,C,NINDPW,GRID) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
 !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,CR,NINDPW,C,GRID) PRIVATE(I,J)
         DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
            DO J=1,NPL
               CR(NINDPW(J),I)=C(J,I)*GRID%FFTSCA(J,2)
            ENDDO
         ENDDO
 !$OMP END PARALLEL DO

!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CR,C,GRID) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
 !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,CR,C,GRID) PRIVATE(I,J)
         DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
            DO J=1,SIZE(GRID%IND_IN_SPHERE)
               CR(GRID%NINDPWCONJG(J),I)=CONJG(C(GRID%IND_IN_SPHERE(J),I))*GRID%FFTSCA(GRID%IND_IN_SPHERE(J),2)
            ENDDO
         ENDDO
 !$OMP END PARALLEL DO


      ELSE

!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CR,C,NINDPW_INV,FFTSCA) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
 !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,CR,NINDPW_INV,C,FFTSCA) PRIVATE(I,J)
         DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
            DO J=1,NPL
               CR(NINDPW_INV(J),I)=CONJG(C(J,I))*FFTSCA(J)
            ENDDO
         ENDDO
 !$OMP END PARALLEL DO

!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CR,C,NINDPW,FFTSCA) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
 !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,CR,NINDPW,C,FFTSCA) PRIVATE(I,J)
         DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
            DO J=1,NPL
               CR(NINDPW(J),I)=C(J,I)*FFTSCA(J)
            ENDDO
         ENDDO
 !$OMP END PARALLEL DO

      ENDIF

      CALL FFT3D_MU(N,CR,LDR,GRID,1)

      

    END SUBROUTINE FFTWAV_USEINV_MU


!************************* SUBROUTINE FFTEXT_MU ************************
!
!***********************************************************************

    SUBROUTINE FFTEXT_MU(NPL,N,NINDPW,CR,LDR,C,LDC,GRID,LADD)

!! USE moffload_struct_def

      USE prec
      USE mgrid_struct_def

      IMPLICIT NONE

      TYPE(grid_3d) :: GRID

      COMPLEX(q) :: CR(LDR,*),C(LDC,*)
      INTEGER    :: NINDPW(NPL)
      INTEGER    :: NPL, N, LDR, LDC
      LOGICAL    :: LADD

! local variables
      INTEGER :: I,J

# 1234


      

      CALL FFT3D_MU(N,CR,LDR,GRID,-1)

      IF (GRID%REAL2CPLX) THEN
         IF (LADD) THEN
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(C,CR,NINDPW,GRID) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
 !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,C,CR,NINDPW,GRID) PRIVATE(I,J)
            DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
               DO J=1,NPL
                  C(J,I)=C(J,I)+CR(NINDPW(J),I)*GRID%FFTSCA(J,1)
               ENDDO
            ENDDO
 !$OMP END PARALLEL DO
         ELSE
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(C,CR,NINDPW,GRID) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
 !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,C,CR,NINDPW,GRID) PRIVATE(I,J)
            DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
               DO J=1,NPL
                  C(J,I)=CR(NINDPW(J),I)*GRID%FFTSCA(J,1)
               ENDDO
            ENDDO
 !$OMP END PARALLEL DO
         ENDIF
      ELSE
         IF (LADD) THEN
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(C,CR,NINDPW) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
 !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,C,CR,NINDPW) PRIVATE(I,J)
            DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
               DO J=1,NPL
                  C(J,I)=C(J,I)+CR(NINDPW(J),I)
               ENDDO
            ENDDO
 !$OMP END PARALLEL DO
         ELSE
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(C,CR,NINDPW) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
 !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,C,CR,NINDPW) PRIVATE(I,J)
            DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
               DO J=1,NPL
                  C(J,I)=CR(NINDPW(J),I)
               ENDDO
            ENDDO
 !$OMP END PARALLEL DO
         ENDIF
      ENDIF

      

    END SUBROUTINE FFTEXT_MU


!************************* SUBROUTINE FFTEXT_USEINV_MU *****************
!
!***********************************************************************

    SUBROUTINE FFTEXT_USEINV_MU(NPL,N,NINDPW,FFTSCA,CR,LDR,C,LDC,GRID,LADD)

!! USE moffload_struct_def

      USE prec
      USE mgrid_struct_def

      IMPLICIT NONE

      TYPE(grid_3d) :: GRID

      COMPLEX(q) :: CR(LDR,*),C(LDC,*)
      REAL(q)    :: FFTSCA(NPL)
      INTEGER    :: NINDPW(NPL)
      INTEGER    :: NPL, N, LDR, LDC
      LOGICAL    :: LADD

! local variables
      INTEGER :: I,J

# 1323


      

      CALL FFT3D_MU(N,CR,LDR,GRID,-1)

      IF (GRID%REAL2CPLX) THEN
         IF (LADD) THEN
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(C,CR,NINDPW,GRID) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
 !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,C,CR,NINDPW,GRID) PRIVATE(I,J)
            DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
               DO J=1,NPL
                  C(J,I)=C(J,I)+CR(NINDPW(J),I)*GRID%FFTSCA(J,1)
               ENDDO
            ENDDO
 !$OMP END PARALLEL DO
         ELSE
!$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW,GRID) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
 !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,C,CR,NINDPW,GRID) PRIVATE(I,J)
            DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
               DO J=1,NPL
                  C(J,I)=CR(NINDPW(J),I)*GRID%FFTSCA(J,1)
               ENDDO
            ENDDO
 !$OMP END PARALLEL DO
         ENDIF
      ELSE
         IF (LADD) THEN
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(C,CR,NINDPW,FFTSCA) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
 !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,C,CR,NINDPW,FFTSCA) PRIVATE(I,J)
            DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
               DO J=1,NPL
                  C(J,I)=C(J,I)+CR(NINDPW(J),I)*FFTSCA(J)
               ENDDO
            ENDDO
 !$OMP END PARALLEL DO
         ELSE
!$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW,FFTSCA) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
 !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,C,CR,NINDPW,FFTSCA) PRIVATE(I,J)
            DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
               DO J=1,NPL
                  C(J,I)=CR(NINDPW(J),I)*FFTSCA(J)
               ENDDO
            ENDDO
 !$OMP END PARALLEL DO
         ENDIF
      ENDIF

      

    END SUBROUTINE FFTEXT_USEINV_MU


!************************* SUBROUTINE FFT3D_MU *************************
!
!***********************************************************************

    SUBROUTINE FFT3D_MU(N,C,LDC,GRID,ISN)
      USE prec
      USE mgrid_struct_def
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE(grid_3d) :: GRID

      COMPLEX(q) :: C(LDC,*)
      INTEGER :: N,LDC,ISN

! local variables
      INTEGER :: NX,NY,NZ

      

      NX=GRID%NGPTAR(1)
      NY=GRID%NGPTAR(2)
      NZ=GRID%NGPTAR(3)

      IF (.NOT. GRID%REAL2CPLX) THEN
!-------------------------------------------------------------------------------
!  complex to complex FFTs
!-------------------------------------------------------------------------------
         IF (.NOT. (NX==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ==GRID%NGZ_rd) ) THEN
            CALL vtutor%bug("FFT3D_MU: grid dimensions not properly set (1) " // &
               str(GRID%LREAL) // " " // str(GRID%REAL2CPLX) // "\n " // &
               str(NX) // " " // str(NY) // " " // str(NZ) // "\n " // &
               str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " // str(GRID%NGZ_rd), &
               "fft_base.F", 1418)
         ENDIF

         IF  (.NOT. GRID%LREAL ) THEN

            CALL FFTBAS_PLAN_MU(N,C,LDC,GRID,ISN)

         ELSE

            CALL vtutor%bug("FFT3D_MU: batched ffts for REAL2CPLX = " // str(GRID%REAL2CPLX) // &
               " and LREAL = " // str(GRID%LREAL) //" not supported yet","fft_base.F",1428)

         ENDIF
      ELSE
!-------------------------------------------------------------------------------
!  real to complex FFTs
!-------------------------------------------------------------------------------
         IF (.NOT. (NX/2+1==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ==GRID%NGZ_rd) ) THEN
            CALL vtutor%bug("FFT3D_MU: grid dimensions not properly set (2)" // &
               str(GRID%LREAL) // " " // str(GRID%REAL2CPLX) // "\n " // &
               str(NX) // " " // str(NY) // " " // str(NZ) // "\n " // &
               str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " //  str(GRID%NGZ_rd), &
                "fft_base.F", 1440)
         ENDIF

         IF (.NOT. GRID%LREAL) THEN
!
! the data in real space is stored as an array of complex numbers (with Im(c)=0)
!
            IF (ISN==1) THEN
! q -> r FFT
               CALL FFTBRC_PLAN_MU(N,C,LDC,GRID,ISN)
! contract columns along x (go from stride NX+2 to stride NX)
! and store real valued transform as a array of complex numbers
! (with Im(c) = 0)
               CALL RESTRIDE_Q2R_CMPLX_MU(N,NY*NZ,NX,C,LDC)
            ELSE
! expand columns along x (go from stride NX to stride NX+2)
! and store array of complex numbers with Im(c)=0 as array
! of reals
               CALL RESTRIDE_R2Q_CMPLX_MU(N,NY*NZ,NX,C,LDC)
! r -> q FFT
               CALL FFTBRC_PLAN_MU(N,C,LDC,GRID,ISN)
            ENDIF

         ELSE
!
! the data in real space is stored as an array of real numbers
!
            IF (ISN==1) THEN
! q -> r FFT
               CALL FFTBRC_PLAN_MU(N,C,LDC,GRID,ISN)
! contract columns along x (go from stride NX+2 to stride NX)
               CALL RESTRIDE_Q2R_MU(N,NY*NZ,NX,C,LDC)
            ELSE
! expand columns along x (go from stride NX to stride NX+2)
               CALL RESTRIDE_R2Q_MU(N,NY*NZ,NX,C,LDC)
! q -> r FFT
               CALL FFTBRC_PLAN_MU(N,C,LDC,GRID,ISN)
            ENDIF

         ENDIF
      ENDIF

      

    END SUBROUTINE FFT3D_MU


!************************* SUBROUTINE RESTRIDE_Q2R_CMPLX_MU ************
!
!> The output of an inplace complex-to-real FFT on a grid of NX*NY*NZ
!> points is an ( NX+2,NY,NZ ) array of real numbers (where the first
!> dimension is padded with zeros.
!>
!> RESTRIDE_Q2R_CMPLX_MU does two things:
!>
!> * The padding is removed, i.e., the data is re-stored as an
!>   ( NX,NY,NZ ) array.
!> * The data is stored as an array of complex numbers (with Im(c)=0).
!
!***********************************************************************

    SUBROUTINE RESTRIDE_Q2R_CMPLX_MU(N,NCOL,NX,C,LDC)

!! USE moffload_struct_def

      USE prec
      IMPLICIT NONE
      INTEGER :: N,NCOL,NX,LDC
      REAL(q) :: C(2*LDC,*)
! local variables
      INTEGER :: ICOL,NDEST,NSRC,I,J

# 1523


!$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(C,CTMP) PRIVATE(NDEST,NSRC) 
 !$OMP PARALLEL DO PRIVATE(J,ICOL,NDEST,NSRC,I)
      DO J=1,N
         DO ICOL=NCOL-1,0,-1
       NDEST=ICOL* NX*2
       NSRC =ICOL*(NX+2)
            DO I=NX,1,-1
!!          NDEST=ICOL* NX*2
!!          NSRC =ICOL*(NX+2)
          C(NDEST+I*2-1,J)=C   (NSRC+I,J)
!!          C(NDEST+I*2-1,J)=CTMP(NSRC+I,J)
               C(NDEST+I*2  ,J)=0
            ENDDO
         ENDDO
      ENDDO
 !$OMP END PARALLEL DO

# 1545

    END SUBROUTINE RESTRIDE_Q2R_CMPLX_MU


!************************* SUBROUTINE RESTRIDE_R2Q_CMPLX_MU ************
!
!> RESTRIDE_R2Q_CMPLX_MU adresses two peculiarities:
!>
!> * The input of an inplace real-to-complex FFT on a grid of NX*NY*NZ
!>   points is an ( NX+2,NY,NZ ) array of real numbers (where the first
!>   dimension is padded with zeros.
!>
!> * In many cases the data (even though it is real) in real space is
!>   stored as an array of complex numbers (with Im(c)=0)
!>
!> RESTRIDE_R2Q_CMPLX_MU does two things:
!>
!> * Padding is added, i.e., the data is re-stored as an ( NX+2,NY,NZ )
!>   array.
!> * The data is stored as an array of real numbers.
!
!***********************************************************************

    SUBROUTINE RESTRIDE_R2Q_CMPLX_MU(N,NCOL,NX,C,LDC)

!! USE moffload_struct_def

      USE prec
      IMPLICIT NONE
      INTEGER :: N,NCOL,NX,LDC
      REAL(q) :: C(2*LDC,*)
! local variables
      INTEGER :: ICOL,NDEST,NSRC,I,J

# 1590


!$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(C,CTMP) PRIVATE(NSRC,NDEST) 
 !$OMP PARALLEL DO PRIVATE(J,ICOL,NSRC,NDEST,I)
      DO J=1,N
         DO ICOL=0,NCOL-1
       NSRC =ICOL* NX*2
       NDEST=ICOL*(NX+2)
            DO I=1,NX
!!          NSRC =ICOL* NX*2
!!          NDEST=ICOL*(NX+2)
          C(NDEST+I,J)=C   (NSRC+I*2-1,J)
!!          C(NDEST+I,J)=CTMP(NSRC+I*2-1,J)
            ENDDO
         ENDDO
      ENDDO
 !$OMP END PARALLEL DO

# 1611

    END SUBROUTINE RESTRIDE_R2Q_CMPLX_MU


!************************* SUBROUTINE RESTRIDE_Q2R_MU ******************
!
!> The output of an inplace complex-to-real FFT on a grid of NX*NY*NZ
!> points is an ( NX+2,NY,NZ ) array of real numbers (where the first
!> dimension is padded with zeros.
!>
!> RESTRIDE_Q2R_MU removes the padding, i.e., the data is re-stored as
!> an ( NX,NY,NZ ) array.
!
!***********************************************************************

    SUBROUTINE RESTRIDE_Q2R_MU(N,NCOL,NX,C,LDC)

!! USE moffload_struct_def

      USE prec
      IMPLICIT NONE
      INTEGER :: N,NCOL,NX,LDC
      REAL(q) :: C(2*LDC,*)
! local variables
      INTEGER :: ICOL,NDEST,NSRC,I,J

# 1648


!$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(C,CTMP) PRIVATE(NDEST,NSRC) 
 !$OMP PARALLEL DO PRIVATE(J,ICOL,NDEST,NSRC,I)
      DO J=1,N
         DO ICOL=1,NCOL-1
       NDEST=ICOL* NX
       NSRC =ICOL*(NX+2)
            DO I=1,NX
!!          NDEST=ICOL* NX
!!          NSRC =ICOL*(NX+2)
          C(NDEST+I,J)=C   (NSRC+I,J)
!!          C(NDEST+I,J)=CTMP(NSRC+I,J)
            ENDDO
         ENDDO
      ENDDO
 !$OMP END PARALLEL DO

# 1669

    END SUBROUTINE RESTRIDE_Q2R_MU


!************************* SUBROUTINE RESTRIDE_R2Q_MU ******************
!
!> The input of an inplace real-to-complex FFT on a grid of NX*NY*NZ
!> points is an ( NX+2,NY,NZ ) array of real numbers (where the first
!> dimension is padded with zeros.
!>
!> RESTRIDE_R2Q_MU add this padding, i.e., the data is re-stored as an
!> ( NX+2,NY,NZ ) array.
!
!***********************************************************************

    SUBROUTINE RESTRIDE_R2Q_MU(N,NCOL,NX,C,LDC)

!! USE moffload_struct_def

      USE prec
      IMPLICIT NONE
      INTEGER :: N,NCOL,NX,LDC
      REAL(q) :: C(2*LDC,*)
! local variables
      INTEGER :: ICOL,NDEST,NSRC,I,J

# 1706


!$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(C,CTMP) PRIVATE(NSRC,NDEST) 
 !$OMP PARALLEL DO PRIVATE(J,ICOL,NSRC,NDEST,I)
      DO J=1,N
         DO ICOL=NCOL-1,1,-1
       NSRC =ICOL* NX
       NDEST=ICOL*(NX+2)
            DO I=NX,1,-1
!!          NSRC =ICOL* NX
!!          NDEST=ICOL*(NX+2)
          C(NDEST+I,J)=C   (NSRC+I,J)
!!          C(NDEST+I,J)=CTMP(NSRC+I,J)
            ENDDO
         ENDDO
      ENDDO
 !$OMP END PARALLEL DO

# 1727

    END SUBROUTINE RESTRIDE_R2Q_MU



!************************* SUBROUTINE FFTGRIDPLAN_MPI ******************
!
!***********************************************************************

    SUBROUTINE FFTGRIDPLAN_MPI(GRID)
      USE prec
      USE mgrid_struct_def

      IMPLICIT NONE

      TYPE (grid_3d) GRID

! local variables
      COMPLEX(q), ALLOCATABLE :: CWORK(:)
      INTEGER, PARAMETER :: NADD=1024, NTIMES=4
      INTEGER :: I

      ALLOCATE(CWORK(GRID%MPLWV+NADD))

      DO I=1,NTIMES
         CALL INIDAT(GRID%MPLWV+NADD,CWORK)
         CALL FFTMAKEPLAN_MPI(CWORK(I),GRID)
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP MASTER
!$       CALL FFTMAKEPLAN_MPI(CWORK(I),GRID)
!$OMP END MASTER
!$OMP END PARALLEL
      ENDDO

      DEALLOCATE(CWORK)

      RETURN
    END SUBROUTINE FFTGRIDPLAN_MPI


!************************* SUBROUTINE FFTINI_MPI ***********************
!
!  if necessary this routine performes initialization
!  for FFTWAV and FFTEXT
!  usually this is only necessary for the Gamma point only
!  1-kpoint version
!
!   FFTSCA(.,1) is the scaling factor for extracting the wavefunction
!               from the FFT grid (FFTEXT)
!   FFTSCA(.,2) is the scaling factor for puting the wavefunction on
!               the grid
!***********************************************************************

    SUBROUTINE FFTINI_MPI(NINDPW,NPLWKP,NKPTS,NRPLW,GRID)
      USE prec
      USE mgrid_struct_def
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      INTEGER :: NKPTS, NRPLW
      INTEGER :: NPLWKP(NKPTS)
      INTEGER :: NINDPW(NRPLW,NKPTS)
! local variables
      REAL(q) :: FACTM
      INTEGER :: NK, NPL, N, IND, N1, NC, N2, N3

      IF (GRID%REAL2CPLX) THEN
         IF (NKPTS>1) THEN
            CALL vtutor%error("FFT3D: real version works only for 1 k-point")
         ENDIF

         NK=1
         NPL=NPLWKP(NK)
         NULLIFY(GRID%FFTSCA)
         ALLOCATE(GRID%FFTSCA(NPL,2))

         DO N=1,NPL
            IND=NINDPW(N,NK)
            N1= MOD((IND-1),GRID%RC%NROW)+1
            NC= (IND-1)/GRID%RC%NROW+1
            N2= GRID%RC%I2(NC)
            N3= GRID%RC%I3(NC)

            FACTM=SQRT(2._q)
            IF (N1==1 .AND. N2==1 .AND. N3==1) FACTM=1
            GRID%FFTSCA(N,1)= FACTM
            GRID%FFTSCA(N,2)= 1/FACTM
! this statment is required
! because for z==0 only half of the FFT components are set
! upon calling FFTWAV
            IF (N3==1) GRID%FFTSCA(N,2)=FACTM
         ENDDO
      END IF
      RETURN
    END SUBROUTINE FFTINI_MPI


!************************* SUBROUTINE FFTWAV_MPI ***********************
!
!  this subroutine transforms a wavefunction C defined  within  the
!  cutoff-sphere to real space CR
! MIND:
! for the real version (gamma point only) it is assumed
! that the wavefunctions at NGZ != 0
! are multiplied by a factor sqrt(2) on the linear grid
! this factor has to be removed before the FFT transformation !
! (scaling with   FFTSCA(M,2))
!
!> @details @ref openmp :
!> the loops that fill CR with the reciprocal space components of
!> an orbital are distributed over all available OpenMP threads.
!
!***********************************************************************

    SUBROUTINE FFTWAV_MPI(NPL,NINDPW,CR,C,GRID)

!! USE moffload_struct_def

      USE prec
      USE mgrid_struct_def
      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      COMPLEX(q) :: C(NPL), CR(GRID%MPLWV)
      INTEGER :: NPL, NINDPW(NPL)
! local variables
      INTEGER :: M

      

      IF (GRID%LREAL) THEN
!DIR$ IVDEP
!OCL NOVREC
!$ACC PARALLEL LOOP PRESENT(CR) IF (OFFLOAD_ON) ASYNC(ACC_ASYNC_Q)
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
 !$OMP PRIVATE(M) SHARED(CR,GRID)
         DO M=1,GRID%RL%NCOL*GRID%NGZ/2
            CR(M)=0.0_q
         ENDDO
 !$OMP END PARALLEL DO
      ELSE
!DIR$ IVDEP
!OCL NOVREC
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
 !$OMP PRIVATE(M) SHARED(CR,GRID)
!$ACC PARALLEL LOOP PRESENT(CR) IF (OFFLOAD_ON) ASYNC(ACC_ASYNC_Q)
         DO M=1,GRID%RL%NCOL*GRID%NGZ
            CR(M)=0.0_q
         ENDDO
 !$OMP END PARALLEL DO
      ENDIF

      IF (GRID%REAL2CPLX) THEN
!DIR$ IVDEP
!OCL NOVREC
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
 !$OMP PRIVATE(M) SHARED(NPL,C,CR,NINDPW,GRID)
!$ACC PARALLEL LOOP PRESENT(CR,C,NINDPW,GRID%FFTSCA) IF (OFFLOAD_ON) ASYNC(ACC_ASYNC_Q)
         DO M=1,NPL
            CR(NINDPW(M))=C(M)*GRID%FFTSCA(M,2)
         ENDDO
 !$OMP END PARALLEL DO
      ELSE
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
 !$OMP PRIVATE(M) SHARED(NPL,C,CR,NINDPW)
!$ACC PARALLEL LOOP PRESENT(CR,C,NINDPW) IF (OFFLOAD_ON) ASYNC(ACC_ASYNC_Q)
         DO M=1,NPL
            CR(NINDPW(M))=C(M)
         ENDDO
 !$OMP END PARALLEL DO
      ENDIF

      CALL FFT3D_MPI(CR,GRID,1)

      

      RETURN
    END SUBROUTINE FFTWAV_MPI


!************************* SUBROUTINE FFTEXT_MPI ***********************
!
! this subroutine performes a FFT to reciprocal space and extracts data
! from the FFT-mesh
! MIND:
! for the real version (gamma point only) it is assumed
! that the wavefunctions at NGX != 0
! are multiplied by a factor sqrt(2) on the linear grid
! this factor has to be applied after the FFT transformation !
!  (scaling with   FFTSCA(M))
!
!> @details @ref openmp :
!> the loops that extract the reciprocal space components of an orbital
!> from the full FFT grid are distributed over all available OpenMP
!> threads.
!
!***********************************************************************

    SUBROUTINE FFTEXT_MPI(NPL,NINDPW,CR,C,GRID,LADD)

!! USE moffload_struct_def

      USE prec
      USE mgrid_struct_def
      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      COMPLEX(q) :: C(NPL),CR(GRID%MPLWV)
      INTEGER :: NPL, NINDPW(NPL)
      LOGICAL :: LADD
! local variables
      INTEGER :: M

      

      CALL FFT3D_MPI(CR,GRID,-1)

      IF (LADD .AND. GRID%REAL2CPLX) THEN
!$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW,GRID%FFTSCA) IF (OFFLOAD_ON) ASYNC(ACC_ASYNC_Q)
!DIR$ IVDEP
!OCL NOVREC
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
 !$OMP PRIVATE(M) SHARED(NPL,C,CR,NINDPW, GRID)
         DO M=1,NPL
            C(M)=C(M)+CR(NINDPW(M))*GRID%FFTSCA(M,1)
         ENDDO
 !$OMP END PARALLEL DO
      ELSE IF (LADD .AND. .NOT. GRID%REAL2CPLX) THEN
!$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW) IF (OFFLOAD_ON) ASYNC(ACC_ASYNC_Q)
!DIR$ IVDEP
!OCL NOVREC
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
 !$OMP PRIVATE(M) SHARED(NPL,CR,NINDPW) REDUCTION(+:C)
         DO M=1,NPL
            C(M)=C(M)+CR(NINDPW(M))
         ENDDO
 !$OMP END PARALLEL DO
      ELSE IF (GRID%REAL2CPLX) THEN
!$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW,GRID%FFTSCA) IF (OFFLOAD_ON) ASYNC(ACC_ASYNC_Q)
!DIR$ IVDEP
!OCL NOVREC
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
 !$OMP PRIVATE(M) SHARED(NPL,C,CR,NINDPW,GRID)
        DO M=1,NPL
          C(M)=CR(NINDPW(M))*GRID%FFTSCA(M,1)
        ENDDO
 !$OMP END PARALLEL DO
      ELSE
!$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW) IF (OFFLOAD_ON) ASYNC(ACC_ASYNC_Q)
!DIR$ IVDEP
!OCL NOVREC
 !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
 !$OMP PRIVATE(M) SHARED(NPL,C,CR,NINDPW)
        DO M=1,NPL
          C(M)=CR(NINDPW(M))
        ENDDO
 !$OMP END PARALLEL DO
      ENDIF

      

      RETURN
    END SUBROUTINE FFTEXT_MPI


!************************* SUBROUTINE FFT3D_MPI ************************
!
!    3-d fast fourier transform (possibly real to complex and vice versa)
!    for chardensities and potentials
!     +1  q->r   vr= sum(q) vq exp(+iqr)    (might be complex to real)
!     -1  r->q   vq= sum(r) vr exp(-iqr)    (might be real to complex)
!
!***********************************************************************

    SUBROUTINE FFT3D_MPI(C,GRID,ISN)
      USE prec
      USE mgrid_struct_def
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      REAL(q) :: C(*)
      INTEGER :: ISN
! local variables
      INTEGER :: NX, NY, NZ

      

      NX=GRID%NGPTAR(1)
      NY=GRID%NGPTAR(2)
      NZ=GRID%NGPTAR(3)

      IF (GRID%RL%NFAST==1) THEN
!-------------------------------------------------------------------------------
! parallel FFT with serial data layout  (GRID%RL%NFAST==1)
!-------------------------------------------------------------------------------
! complex to complex case
         IF (.NOT. GRID%REAL2CPLX) THEN
            IF (.NOT. (NX==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ==GRID%NGZ_rd) ) THEN
               CALL vtutor%bug("internal error 1 in FFT3D_MPI: something not properly set " // &
                  str(GRID%LREAL) // " " // str(GRID%REAL2CPLX) // "\n " // str(NX) // " " // str(NY) &
                  // " " // str(NZ) // "\n " // str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " // &
                  str(GRID%NGZ_rd), "fft_base.F", 2031)
            ENDIF

            IF (ISN==1) THEN
! q->r FFT
               CALL FFTBAS_PLAN_MPI(C,GRID,ISN)
! bring real space data from parallel to serial data layout
               IF (GRID%COMM%NODE_ME==1) THEN
                  CALL FFTPAR_TO_SER(GRID%NGX, GRID%NGY, GRID%NGZ, C)
               ENDIF
            ELSE
! bring real space data from serial to parallel data layout
               IF (GRID%COMM%NODE_ME==1) THEN
                  CALL FFTSER_TO_PAR(GRID%NGX, GRID%NGY, GRID%NGZ, C)
               ENDIF
! r->q FFT
               CALL FFTBAS_PLAN_MPI(C,GRID,ISN)
            ENDIF

! complex to real case
         ELSE IF (GRID%LREAL) THEN
            IF (.NOT. (NX==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ/2+1==GRID%NGZ_rd) ) THEN
               CALL vtutor%bug("internal error 2 in FFT3D_MPI: something not properly set " // &
                  str(GRID%LREAL) // " " // str(GRID%REAL2CPLX) // "\n " // str(NX) // " " // str(NY) &
                  // " " // str(NZ) // "\n " // str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " // &
                  str(GRID%NGZ_rd), "fft_base.F", 2056)
            ENDIF

!  in real space the first dimension in VASP is NGZ (REAL data)
!  but the FFT requires NGZ+2 (real data)
!  therefore some data movement is required
            IF (ISN==1) THEN
! q->r FFT
               CALL FFTBRC_PLAN_MPI(C,GRID,ISN)
! x-lines (go from stride NZ+2 to NZ)
               CALL RESTRIDE_Q2R(GRID%RL_FFT%NCOL,NZ,C)
! bring real space data from parallel to serial data layout
               IF (GRID%COMM%NODE_ME==1) THEN
                  CALL FFTPAR_TO_SER_REAL(GRID%NGX, GRID%NGY, GRID%NGZ, C)
               ENDIF
            ELSE
! bring real space data from serial to parallel data layout
               IF (GRID%COMM%NODE_ME==1) THEN
                  CALL FFTSER_TO_PAR_REAL(GRID%NGX, GRID%NGY, GRID%NGZ, C)
               ENDIF
! x-lines (go from stride NZ to NZ+2)
               CALL RESTRIDE_R2Q(GRID%RL_FFT%NCOL,NZ,C)
! r->q FFT
               CALL FFTBRC_PLAN_MPI(C,GRID,ISN)
            ENDIF
         ELSE
            CALL vtutor%error("ERROR in FFT3D_MPI: this version does not support the required half &
               &grid mode \n " // str(NX) // " " // str(NY) // " " // str(NZ) // "\n " // &
               str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " // str(GRID%NGZ_rd))
         ENDIF
!-------------------------------------------------------------------------------
!  complex parallel FFT
!-------------------------------------------------------------------------------
      ELSE IF (.NOT. GRID%REAL2CPLX) THEN
         IF (.NOT. (NX==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ==GRID%NGZ_rd) ) THEN
            CALL vtutor%bug("internal error 3 in FFT3D_MPI: something not properly set " // &
               str(GRID%LREAL) // " " // str(GRID%REAL2CPLX) // "\n " // str(NX) // " " // str(NY) // &
               " " // str(NZ) // "\n " // str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " // &
               str(GRID%NGZ_rd), "fft_base.F", 2094)
         ENDIF
         CALL FFTBAS_PLAN_MPI(C,GRID,ISN)
!-------------------------------------------------------------------------------
!  real to complex parallel FFT
!-------------------------------------------------------------------------------
      ELSE IF (GRID%LREAL) THEN
         IF (.NOT.(NX==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ/2+1==GRID%NGZ_rd) ) THEN
            CALL vtutor%bug("internal error 4 in FFT3D_MPI: something not properly set " // &
               str(GRID%LREAL) // " " // str(GRID%REAL2CPLX) // "\n " // str(NX) // " " // str(NY) // &
               " " // str(NZ) // "\n " // str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " // &
               str(GRID%NGZ_rd), "fft_base.F", 2105)
         ENDIF

!  in real space the first dimension in VASP is NGZ (REAL data)
!  but the FFT requires NGZ+2 (real data)
!  therefore some data movement is required
         IF (ISN==1) THEN
! q->r FFT
            CALL FFTBRC_PLAN_MPI(C,GRID,ISN)
! concat z-lines (go from stride NZ+2 to NZ)
            CALL RESTRIDE_Q2R(GRID%RL%NCOL,NZ,C)
         ELSE
! concat z-lines (go from stride NZ to NZ+2)
            CALL RESTRIDE_R2Q(GRID%RL%NCOL,NZ,C)
! r->q FFT
            CALL FFTBRC_PLAN_MPI(C,GRID,ISN)
         ENDIF
!-------------------------------------------------------------------------------
!  real to complex parallel FFT with complex storage layout in real space
!-------------------------------------------------------------------------------
      ELSE
         IF (.NOT.(NX==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ/2+1==GRID%NGZ_rd) ) THEN
            CALL vtutor%bug("internal error 5 in FFT3D_MPI: something not properly set " // &
               str(GRID%LREAL) // " " // str(GRID%REAL2CPLX) // "\n " // str(NX) // " " // str(NY) // &
               " " // str(NZ) // "\n " // str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " // &
               str(GRID%NGZ_rd), "fft_base.F", 2130)
         ENDIF

         IF (ISN==1) THEN
! q->r FFT
            CALL FFTBRC_PLAN_MPI(C,GRID,ISN)
! concat z-lines (go from stride NZ+2 to NZ)
            CALL RESTRIDE_Q2R_CMPLX(GRID%RL%NCOL,NZ,C)
         ELSE
! z-lines (go from complex stride NZ to real stride NZ+2)
            CALL RESTRIDE_R2Q_CMPLX(GRID%RL%NCOL,NZ,C)
! r->q FFT
            CALL FFTBRC_PLAN_MPI(C,GRID,ISN)
         ENDIF
      ENDIF

      

    END SUBROUTINE FFT3D_MPI


!************************* SUBROUTINE FFTPAR_TO_SER ********************
!
! change data layout from parallel to serial data layout
! and vice versa for complex and real arrays
! operates usually in real space
!
!***********************************************************************

    SUBROUTINE FFTPAR_TO_SER(NGX, NGY, NGZ, CORIG)

!! USE moffload_struct_def

      USE prec

      INTEGER NGX, NGY, NGZ
      COMPLEX(q) :: CORIG(NGX*NGY*NGZ)
! local
      INTEGER IX, IY, IZ
      COMPLEX(q) :: C(NGX*NGY*NGZ)

!$ACC ENTER DATA CREATE(C) 
!$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(C,CORIG) 
      DO IX=0,NGX-1
         DO IY=0,NGY-1
!DIR$ IVDEP
!OCL NOVREC
            DO IZ=0,NGZ-1
! C(IX,IY,IZ)=CORIG(IZ,IX,IY)
               C(1+IX+NGX*(IY+NGY*IZ))=CORIG(1+IZ+NGZ*(IX+NGX*IY))
            ENDDO
         ENDDO
      ENDDO
!DIR$ IVDEP
!OCL NOVREC
!$ACC PARALLEL LOOP PRESENT(C,CORIG) 
      DO IX=1,NGX*NGY*NGZ
         CORIG(IX)=C(IX)
      ENDDO
!$ACC EXIT DATA DELETE(C) 
!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
    END SUBROUTINE FFTPAR_TO_SER


    SUBROUTINE FFTPAR_TO_SER_REAL(NGX, NGY, NGZ, CORIG)

!! USE moffload_struct_def

      USE prec

      INTEGER NGX, NGY, NGZ
      REAL(q) :: CORIG(NGX*NGY*NGZ)
! local
      INTEGER IX, IY, IZ
      REAL(q) :: C(NGX*NGY*NGZ)

!$ACC ENTER DATA CREATE(C) 
!$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(C,CORIG) 
      DO IX=0,NGX-1
         DO IY=0,NGY-1
!DIR$ IVDEP
!OCL NOVREC
            DO IZ=0,NGZ-1
! C(IX,IY,IZ)=CORIG(IZ,IX,IY)
               C(1+IX+NGX*(IY+NGY*IZ))=CORIG(1+IZ+NGZ*(IX+NGX*IY))
            ENDDO
         ENDDO
      ENDDO
!DIR$ IVDEP
!OCL NOVREC
!$ACC PARALLEL LOOP PRESENT(C,CORIG) 
      DO IX=1,NGX*NGY*NGZ
         CORIG(IX)=C(IX)
      ENDDO
!$ACC EXIT DATA DELETE(C) 
!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
    END SUBROUTINE FFTPAR_TO_SER_REAL


    SUBROUTINE FFTSER_TO_PAR(NGX, NGY, NGZ, CORIG)

!! USE moffload_struct_def

      USE prec

      INTEGER NGX, NGY, NGZ
      COMPLEX(q) :: CORIG(NGX*NGY*NGZ)
! local
      INTEGER IX, IY, IZ
      COMPLEX(q) :: C(NGX*NGY*NGZ)

!$ACC ENTER DATA CREATE(C) 
!$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(C,CORIG) 
      DO IX=0,NGX-1
         DO IY=0,NGY-1
!DIR$ IVDEP
!OCL NOVREC
            DO IZ=0,NGZ-1
! C(IZ,IX,IY)=CORIG(IX,IY,IZ)
               C(1+IZ+NGZ*(IX+NGX*IY))=CORIG(1+IX+NGX*(IY+NGY*IZ))
            ENDDO
         ENDDO
      ENDDO
!DIR$ IVDEP
!OCL NOVREC
!$ACC PARALLEL LOOP PRESENT(C,CORIG) 
      DO IX=1,NGX*NGY*NGZ
         CORIG(IX)=C(IX)
      ENDDO
!$ACC EXIT DATA DELETE(C) 
!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
    END SUBROUTINE FFTSER_TO_PAR


    SUBROUTINE FFTSER_TO_PAR_REAL(NGX, NGY, NGZ, CORIG)

!! USE moffload_struct_def

      USE prec

      INTEGER NGX, NGY, NGZ
      REAL(q) :: CORIG(NGX*NGY*NGZ)
! local
      INTEGER IX, IY, IZ
      REAL(q) :: C(NGX*NGY*NGZ)

!$ACC ENTER DATA CREATE(C) 
!$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(C,CORIG) 
      DO IX=0,NGX-1
         DO IY=0,NGY-1
!DIR$ IVDEP
!OCL NOVREC
            DO IZ=0,NGZ-1
! C(IZ,IX,IY)=CORIG(IX,IY,IZ)
               C(1+IZ+NGZ*(IX+NGX*IY))=CORIG(1+IX+NGX*(IY+NGY*IZ))
            ENDDO
         ENDDO
      ENDDO
!DIR$ IVDEP
!OCL NOVREC
!$ACC PARALLEL LOOP PRESENT(C,CORIG) 
      DO IX=1,NGX*NGY*NGZ
         CORIG(IX)=C(IX)
      ENDDO
!$ACC EXIT DATA DELETE(C) 
!$ACC WAIT(ACC_ASYNC_Q) IF(OFFLOAD_ON)
    END SUBROUTINE FFTSER_TO_PAR_REAL


