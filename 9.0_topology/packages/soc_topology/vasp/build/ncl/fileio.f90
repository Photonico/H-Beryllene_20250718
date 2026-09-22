# 1 "fileio.F"
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


# 2 "fileio.F" 2 
!***********************************************************************
! RCS:  $Id: fileio.F,v 1.12 2003/06/27 13:22:18 kresse Exp kresse $
!
!> collection of subroutines to perform File IO
!>
!>  not yet really supported heavily
!>  but we want to switch all IO from the main program to this module
!
!***********************************************************************
      MODULE fileio
      USE prec
      IMPLICIT NONE

!
!> unit for the density matrix GAMMA (file GAMMA)
!
      INTEGER :: IU_GAMMA=77

!*************************SUBROUTINE OUTPOT ****************************
!>   write potential     to a specified unit
!>   HEADER is currently created in the main program
!***********************************************************************

      INTERFACE
         SUBROUTINE OUTPOT(GRIDC,IU,LLONG,CVTOT)
         USE prec
         USE mgrid
         TYPE (grid_3d), INTENT(IN):: GRIDC
         INTEGER, INTENT(IN) :: IU
         LOGICAL, INTENT(IN) :: LLONG
         COMPLEX(q), INTENT(IN) :: CVTOT(GRIDC%RC%NP)
         END SUBROUTINE
      END INTERFACE

      CONTAINS


!*************************SUBROUTINE OPENWAV ***************************
!
!>   open WAVECAR file
!>
!>   optionally an three character extension can be supplied allowing the
!>   open an alternative WAVECAR file
!>
!> - the unit 12 is always connected with the WAVECAR file
!> - for the mpi chain version (repeated images) the files
!>   XXX/WAVECAR are opended
!> - the record length is read from the WAVECAR file
!
!***********************************************************************

      SUBROUTINE OPENWAV(IO, COM, EXT)

      USE base
      USE main_mpi
      USE mpimy

      USE vhdf5_base

      IMPLICIT NONE

      TYPE (in_struct) IO
      TYPE (communic) COM
      LOGICAL :: JUNK
      REAL(q) :: RDUM, RISPIN
      INTEGER :: IDUM
      INTEGER(qi8) :: WAVECAR_SIZE
      CHARACTER (LEN=3), OPTIONAL :: EXT


      JUNK=.TRUE.
      IF (PRESENT(EXT)) THEN
         INQUIRE(FILE=DIR_APP(1:DIR_LEN)//'WAVECAR.'//EXT,EXIST=junk,SIZE=WAVECAR_SIZE)
      ELSE
         INQUIRE(FILE=DIR_APP(1:DIR_LEN)//'WAVECAR',EXIST=junk,SIZE=WAVECAR_SIZE)
      ENDIF
      CALL M_bcast_l(COM, JUNK, 1 )
      IO%LFOUND=JUNK

      IF (JUNK.AND.WAVECAR_SIZE>0) THEN
        WAVECAR_FOUND=.TRUE.
      ELSE
        INQUIRE(FILE=DIR_APP(1:DIR_LEN)//'vaspwave.h5',EXIST=VASPWAVE_FOUND)
        IF (VASPWAVE_FOUND) THEN
          IO%LFOUND=.TRUE.
        ENDIF
      ENDIF


! first reopen with assumed (wrong) record length ICMPLX
      IF (PRESENT(EXT)) THEN
         OPEN(UNIT=12,FILE=DIR_APP(1:DIR_LEN)//'WAVECAR.'//EXT,ACCESS='DIRECT', &
              FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=IO%ICMPLX)
      ELSE
         OPEN(UNIT=12,FILE=DIR_APP(1:DIR_LEN)//'WAVECAR',ACCESS='DIRECT', &
              FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=IO%ICMPLX)
      ENDIF
! the first record contains the record length, get it ...
      RDUM=0._q
      READ(12,REC=1,ERR=17421) RDUM,RISPIN ; IDUM=NINT(RDUM)
! in the worst case IDUM could contain completely useless data and useless is
! all <=0 or all >10000000 (since on output we use RECL=IO%ICMPLX*MAX(NRPLWV,6)
! or RECL=(NB_TOT+2)*ICMPLX more than ten millions sounds not very realistic)
      IF (IDUM<=0) IDUM=IO%ICMPLX  ! -> error reading WAVECAR
      GOTO 17422

17421 CONTINUE
      IDUM=IO%ICMPLX  ! -> error reading WAVECAR
17422 CONTINUE
      CLOSE(12)
      IO%IRECLW=IDUM
! reopen with correct record length (clumsy all that, I know ...)
      IF (PRESENT(EXT)) THEN
         OPEN(UNIT=12,FILE=DIR_APP(1:DIR_LEN)//'WAVECAR.'//EXT,ACCESS='DIRECT', &
              FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=IO%IRECLW)
      ELSE
         OPEN(UNIT=12,FILE=DIR_APP(1:DIR_LEN)//'WAVECAR',ACCESS='DIRECT', &
              FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=IO%IRECLW)
      ENDIF

      END SUBROUTINE OPENWAV

!*************************SUBROUTINE CLOSEWAV **************************
!
!>   close unit 12 (usually connected to the  WAVECAR file)
!
!***********************************************************************

      SUBROUTINE CLOSEWAV
        CLOSE(12)
      END SUBROUTINE CLOSEWAV

!*************************SUBROUTINE INWAV  ****************************
!
!>   read wavefunctions header from file
!>
!>   just parse the old cell shape old cutoff etc.
!>   this is a new version that uses only real quantities and allows
!>   the WAVECAR file to be exchanged between different machines
!>   more easily (IEEE standard)
!
!***********************************************************************

      SUBROUTINE INWAV_HEAD(WDES, LATT_INI, LATT_CUR, ENMAXI, ISTART, IU0)
      USE prec
      USE wave
      USE lattice

      TYPE (wavedes)  WDES
      TYPE (latt)     LATT_INI,LATT_CUR
      INTEGER NODE_ME, IONODE, ISTART, I, J, IU0, NKPTSF, NBANDF, NALLOC
      REAL(q) RKPTSF, RBANDF, ENMAXI
      LOGICAL LDIFF

      NODE_ME=WDES%COMM%NODE_ME
      IONODE =WDES%COMM%IONODE


      READ(12,REC=2,ERR=100) RKPTSF,RBANDF,ENMAXI, &
     &                       ((LATT_INI%A(I,J),I=1,3),J=1,3)

      IF (IU0 >= 0) WRITE(IU0,*) 'found WAVECAR, reading the header'

      NKPTSF=NINT(RKPTSF)
      NBANDF=NINT(RBANDF)

      CALL LATTIC(LATT_INI)

      IF (ISTART==2 .AND. ENMAXI /= WDES%ENMAX) THEN
        CALL vtutor%error("ERROR: ENMAX changed please set ISTART to 1")
      ENDIF

      IF (.NOT. (NBANDF==WDES%NB_TOT .OR. (WDES%NRSPINORS==2 .AND.  NBANDF*WDES%NRSPINORS==WDES%NB_TOT))) THEN
         IF (IU0 >= 0) WRITE(IU0,'(2X,A,I6,A,I6)') &
              'number of bands has changed, file:',NBANDF,' present:',WDES%NB_TOT
         IF (IU0 >= 0) WRITE(IU0,'(2X,A,I6,A,I6)') &
              'trying to continue reading WAVECAR, but it might fail'
      ENDIF
      IF (NKPTSF/=WDES%NKPTS) THEN
         IF (IU0 >= 0) WRITE(IU0,'(2X,A,I6,A,I6)') &
              'number of k-points has changed, file:',NKPTSF,' present:',WDES%NKPTS
         IF (IU0 >= 0) WRITE(IU0,'(2X,A,I6,A,I6)') &
              'trying to continue reading WAVECAR, but it might fail'
      ENDIF
      IF (ISTART ==1) THEN

      LDIFF=.FALSE.
      DO I=1,3
        DO J=1,3
          IF (ABS(LATT_INI%A(I,J)-LATT_CUR%A(I,J)) > 1E-4) LDIFF=.TRUE.
        ENDDO
      ENDDO
      IF (ENMAXI /= WDES%ENMAX) LDIFF=.TRUE.
      IF (LDIFF) THEN
        IF (IU0>=0) &
        WRITE(IU0,*)'WAVECAR: different cutoff or change in lattice found'
      ENDIF

      ENDIF

      RETURN
!
! WAVECAR does not exist or can not be read at all
!
  100 CONTINUE
      LATT_INI%A=LATT_CUR%A
      CALL LATTIC(LATT_INI)

      IF (ISTART==3) THEN
        CALL vtutor%error("ERROR: can't restart (ISTART=3) with wavefunctions on file")
      ELSE
        ISTART=0
      ENDIF
      RETURN

      END SUBROUTINE


!*************************SUBROUTINE INWAV_FAST ************************
!
!>   read wavefunctions from file
!>
!>   the routine can read WAVECAR files generated for different
!>   lattice vectors and cutoff (both in the serial and parallel case)
!>
!>   it is also possible to restart a spin polarised calculation
!>   from a non spin polarised wavefunction file or
!>   to restart a non collinear calculation from a collinear calculation
!
!***********************************************************************

      SUBROUTINE INWAV_FAST(IO, WDES, W, GRID, LATT_CUR, LATT_INI, ISTART, EFERMI)
      USE prec
      USE base
      USE wave_high
      USE mgrid
      USE lattice
      USE string, ONLY: str
      USE tutor, ONLY: vtutor

      TYPE (in_struct)   IO
      TYPE (wavedes)  WDES
      TYPE (wavedes1) WDES1
      TYPE (wavespin) W
      TYPE (grid_3d)  GRID
      TYPE (latt)     LATT_CUR,LATT_INI
      REAL (q) :: EFERMI

      INTEGER NODE_ME, IONODE, IU0
      INTEGER ISPINREAD, IREC, ISP, IFAIL, I, J, K, IRECL, ISTART, IS
      INTEGER NALLOC, MALLOC, NBANDF, NKPTSF, NPL, NPLREAD, NCOMM
      REAL(q) RDUM, RKPTSF, RBANDF, RISPIN, RTAG, RNPL, ENMAXF

! local work arrays
      REAL(q) VKPT(3)
      INTEGER,ALLOCATABLE    :: IND(:),INDI(:)
      COMPLEX(q),ALLOCATABLE :: CW1(:),CW2(:)
      LOGICAL ::  SINGLE_PREC, SPAN_RECORDS
      COMPLEX(qs),  ALLOCATABLE :: CRD(:)
      LOGICAL, POINTER :: LCONJG(:)=>null()
      INTEGER, POINTER :: INDEX(:)=>null()
      REAL(q), ALLOCATABLE :: INBUF(:)
      CHARACTER(LEN=:), ALLOCATABLE :: MSG

      IU0=IO%IU0

      NALLOC = MAXVAL(WDES%NPLWKP_TOT)

      NODE_ME=0
      IONODE =0

      NODE_ME=WDES%COMM%NODE_ME
      IONODE =WDES%COMM%IONODE

!
! parse the header
!
        RTAG=0
        READ(12,REC=1,ERR=200) RDUM,RISPIN,RTAG
        IF (RTAG==53300) THEN
           IF (IU0>=0) WRITE(IU0,*) "VASP.5.3 WAVECAR encountered"
           SINGLE_PREC=.TRUE.
           SPAN_RECORDS=.TRUE.
        ELSE IF (RTAG==53310) THEN
           IF (IU0>=0) WRITE(IU0,*) "VASP.5.3 double precision WAVECAR encountered"
           SINGLE_PREC=.FALSE.
           SPAN_RECORDS=.TRUE.
        ELSE IF (RTAG==45200) THEN
           SINGLE_PREC=.TRUE.
           SPAN_RECORDS=.FALSE.
        ELSE
           IF (IU0>=0) WRITE(IU0,*) "double precision WAVECAR encountered, converting it"
           SINGLE_PREC=.FALSE.
           SPAN_RECORDS=.FALSE.
        ENDIF

        ISPINREAD=NINT(RISPIN)
        READ(12,REC=2,ERR=200) RKPTSF,RBANDF,ENMAXF, &
     &                         ((LATT_INI%A(I,J),I=1,3),J=1,3),EFERMI

        IREC=2
        NKPTSF=NINT(RKPTSF)
        NBANDF=NINT(RBANDF)

        CALL LATTIC(LATT_INI)
!=======================================================================
! read WAVECAR file, number of bands agree
!=======================================================================
      IF (.NOT.(WDES%NRSPINORS==2 .AND. NBANDF*WDES%NRSPINORS == WDES%NB_TOT)) THEN

        spin:    DO ISP=1, MIN(WDES%ISPIN, ISPINREAD)
        kpoints: DO K=1,WDES%NKPTS

        IF ( K> NKPTSF ) THEN
           IREC=IRECL
        ELSE
           IRECL=IREC
        ENDIF

        CALL SETWDES(WDES,WDES1,K)
# 324


        IF (NODE_ME==IONODE) THEN
        IFAIL=0
        IF (.NOT. SPAN_RECORDS) THEN
           IREC=IREC+1
           READ(12,REC=IREC,ERR=230) RNPL,VKPT, &
                (W%CELTOT(J,K,ISP),W%FERTOT(J,K,ISP),J=1,MIN(WDES%NB_TOT,NBANDF))
        ELSE
           ALLOCATE(INBUF(4+NBANDF*3))
           CALL READ_TO_BUF(12, IREC, IO%IRECLW/IO%ICMPLX*2, INBUF, IFAIL)
           RNPL   =INBUF(1)
           VKPT(1)=INBUF(2)
           VKPT(2)=INBUF(3)
           VKPT(3)=INBUF(4)
           DO J=1,MIN(WDES%NB_TOT,NBANDF)
              W%CELTOT(J,K,ISP)=CMPLX(INBUF(5+(J-1)*3),INBUF(6+(J-1)*3),q)
              W%FERTOT(J,K,ISP)=INBUF(7+(J-1)*3)
           ENDDO
           DEALLOCATE(INBUF)

        ENDIF
        NPLREAD=NINT(RNPL)
        ENDIF

        CALL M_bcast_i( WDES%COMM, IFAIL, 1)
        IF (IFAIL /=0) GOTO 230

        CALL M_bcast_i( WDES%COMM, IREC, 1)
        CALL M_bcast_i( WDES%COMM, NPLREAD, 1)
        NPL=WDES%NPLWKP_TOT(K)

        MALLOC=MAX(NPL, NPLREAD)
        ALLOCATE(CW1(MALLOC),CW2(MALLOC),CRD(MALLOC),IND(MALLOC),INDI(MALLOC))

        IF (NODE_ME==IONODE) THEN
! create an index to allow for change of cutoff or cell size
        IF (ISTART==2) THEN
! gK: fix ISTART=2 in INWAV_FAST
           CALL REPAD_INDEX_ARRAY(GRID, WDES%VKPT(:,K), VKPT, LATT_INI%B,  LATT_INI%B, &
              WDES%ENMAX, ENMAXF, NPL/WDES%NRSPINORS, NPLREAD/WDES%NRSPINORS, IND, INDI, MALLOC, IFAIL )
        ELSE
           CALL REPAD_INDEX_ARRAY(GRID, WDES%VKPT(:,K), VKPT, LATT_CUR%B,  LATT_INI%B, &
              WDES%ENMAX, ENMAXF, NPL/WDES%NRSPINORS, NPLREAD/WDES%NRSPINORS, IND, INDI, MALLOC, IFAIL )
        ENDIF
        ENDIF

        CALL M_bcast_i( WDES%COMM, IFAIL, 1)
        IF (IFAIL /=0) GOTO 220

        band: DO J=1,NBANDF
          IREC=IREC+1
          IF (J>WDES%NB_TOT) CYCLE
          IF (NODE_ME==IONODE) THEN
            IF (SINGLE_PREC) THEN
               READ(12,REC=IREC,ERR=240) (CRD(I),I=1,NPLREAD)
               CW2(1:NPLREAD)=CRD(1:NPLREAD)
            ELSE
               READ(12,REC=IREC,ERR=240) (CW2(I),I=1,NPLREAD)
            ENDIF

            CW1=0
            DO IS=1,WDES%NRSPINORS
! store the wave function coefficients according to new cutoff
               CALL REPAD_WITH_INDEX_ARRAY( MALLOC, IND, INDI, &
                  CW1((IS-1)*NPL/WDES%NRSPINORS+1), CW2((IS-1)*NPLREAD/WDES%NRSPINORS+1))
            ENDDO

          ENDIF

          IF (WDES%COMM_KINTER%NCPU.GT.1) THEN
             CALL M_bcast_z(WDES%COMM_KINTER,CW1,SIZE(CW1))
          END IF

# 401

          CALL DIS_PW_BAND(WDES1, J, CW1, W%CPTWFP(1,1,K,ISP))

        ENDDO band

# 408

        DEALLOCATE(CW1,CW2,CRD,IND,INDI)

        ENDDO kpoints

        IF (NKPTSF > WDES%NKPTS) THEN
           IREC=IREC+(NKPTSF-WDES%NKPTS)*(WDES%NB_TOT+1)
        ENDIF

! copy eigenvalues and weights to all nodes
        NCOMM=WDES%NB_TOT*WDES%NKPTS
        CALL M_bcast_z(WDES%COMM, W%CELTOT(1,1,ISP),NCOMM )
        CALL M_bcast_d(WDES%COMM, W%FERTOT(1,1,ISP),NCOMM )

        ENDDO spin

        IF (ISPINREAD > WDES%ISPIN .AND. IU0>=0 ) THEN
           WRITE(IU0,*) 'down-spin wavefunctions not read'
        ENDIF

        IF (NBANDF<WDES%NB_TOT) THEN
           IF (IU0>=0) WRITE(IU0,*) 'random initialization beyond band ',NBANDF
           CALL WFINIT(WDES, W, 1E10_q, NBANDF+1) ! ENINI=1E10 not cutoff restriction
        ENDIF

        IF (IU0 >= 0) WRITE(IU0,*) 'the WAVECAR file was read successfully'

        IF (WDES%ISPIN<=ISPINREAD) THEN
           RETURN
        ENDIF
!
!  spin down is missing
!
        IF (IU0>=0) &
        WRITE(IU0,*) 'No down-spin wavefunctions found', &
     &             ' --> setting down-spin equal up-spin ...'

        DO K=1,WDES%NKPTS
          W%CELTOT(1:WDES%NB_TOT,K,2)=W%CELTOT(1:WDES%NB_TOT,K,1)
          W%FERTOT(1:WDES%NB_TOT,K,2)=W%FERTOT(1:WDES%NB_TOT,K,1)
          NPL=WDES%NPLWKP(K)
          W%CPTWFP(1:NPL,1:WDES%NBANDS,K,2)=W%CPTWFP(1:NPL,1:WDES%NBANDS,K,1)
        ENDDO


        RETURN
!=======================================================================
! read collinear WAVECAR file for a non collinear run
!=======================================================================
      ELSE
        W%CPTWFP=0

        spin2:    DO ISP=1,ISPINREAD
        IF (IU0>=0.AND. ISP==1) &
        WRITE(IU0,*) 'reading wavefunctions of collinear run, up'
        IF (IU0>=0.AND. ISP==2) &
        WRITE(IU0,*) 'reading wavefunctions of collinear run, down'
        kpoints2: DO K=1,NKPTSF

        IF ( K> NKPTSF ) THEN
           IREC=IRECL
        ELSE
           IRECL=IREC
        ENDIF

        CALL SETWDES(WDES,WDES1,K)

        IF (NODE_ME==IONODE) THEN
        IFAIL=0
        IF (.NOT. SPAN_RECORDS) THEN
           IREC=IREC+1
           READ(12,REC=IREC,ERR=230) RNPL, VKPT(1:3), &
                       (W%CELTOT(J+NBANDF*(ISP-1),K,1),W%FERTOT(J+NBANDF*(ISP-1),K,1),J=1,NBANDF)
        ELSE
           ALLOCATE(INBUF(4+NBANDF*3))

           CALL READ_TO_BUF(12, IREC, IO%IRECLW/IO%ICMPLX*2, INBUF, IFAIL)
           RNPL   =INBUF(1)
           VKPT(1)=INBUF(2)
           VKPT(2)=INBUF(3)
           VKPT(3)=INBUF(4)
           DO J=1,NBANDF
              W%CELTOT(J+NBANDF*(ISP-1),K,1)=CMPLX(INBUF(5+(J-1)*3),INBUF(6+(J-1)*3),q)
              W%FERTOT(J+NBANDF*(ISP-1),K,1)=INBUF(7+(J-1)*3)
           ENDDO
           DEALLOCATE(INBUF)

        ENDIF
        NPLREAD=NINT(RNPL)
        ENDIF

        CALL M_bcast_i( WDES%COMM, IFAIL, 1)
        IF (IFAIL /=0) GOTO 230

        CALL M_bcast_i( WDES%COMM, IREC, 1)
        CALL M_bcast_i( WDES%COMM, NPLREAD, 1)
        NPL=WDES%NPLWKP_TOT(K)/2

        MALLOC=MAX(NPL, NPLREAD)
        ALLOCATE(CW1(2*MALLOC),CW2(2*MALLOC),CRD(2*MALLOC),IND(MALLOC),INDI(MALLOC))

        IF (NODE_ME==IONODE) THEN
        IF (ISTART==2) THEN
! gK: fix ISTART=2 in INWAV_FAST
           CALL REPAD_INDEX_ARRAY(GRID, WDES%VKPT(:,K), VKPT, LATT_INI%B,  LATT_INI%B, &
                WDES%ENMAX, ENMAXF, NPL, NPLREAD, IND, INDI, MALLOC, IFAIL )
        ELSE
           CALL REPAD_INDEX_ARRAY(GRID, WDES%VKPT(:,K), VKPT, LATT_CUR%B,  LATT_INI%B, &
                WDES%ENMAX, ENMAXF, NPL, NPLREAD, IND, INDI, MALLOC, IFAIL )
        ENDIF
        ENDIF

        CALL M_bcast_i( WDES%COMM, IFAIL, 1)
        IF (IFAIL /=0) GOTO 220

        IF (ISP==1) THEN
           DO J=1,NBANDF
              IREC=IREC+1
              IF (NODE_ME==IONODE) THEN
                IF (SINGLE_PREC) THEN
                   READ(12,REC=IREC,ERR=240) (CRD(I),I=1,NPLREAD)
                   CW2(1:NPLREAD)=CRD(1:NPLREAD)
                ELSE
                   READ(12,REC=IREC,ERR=240) (CW2(I),I=1,NPLREAD)
                ENDIF


                CW1=0
                CALL REPAD_WITH_INDEX_ARRAY( MALLOC, IND, INDI, CW1, CW2)
              ENDIF


              IF (WDES%COMM_KINTER%NCPU.GT.1) THEN
                 CALL M_bcast_z(WDES%COMM_KINTER,CW1,SIZE(CW1))
              END IF

              CALL DIS_PW_BAND(WDES1, J, CW1, W%CPTWFP(1,1,K,1))

! copy immediately to second panel
              W%CELTOT(J+NBANDF,K,1)=W%CELTOT(J,K,1)
              W%FERTOT(J+NBANDF,K,1)=W%FERTOT(J,K,1)

              CW1(NPL+1:2*NPL)=CW1(1:NPL)
              CW1(1:NPL)=0

              CALL DIS_PW_BAND(WDES1, J+NBANDF, CW1, W%CPTWFP(1,1,K,1))
           ENDDO
        ELSE
           DO J=1,NBANDF
              IREC=IREC+1
              IF (NODE_ME==IONODE) THEN
                IF (SINGLE_PREC) THEN
                   READ(12,REC=IREC,ERR=240) (CRD(I),I=1,NPLREAD)
                   CW2(1:NPLREAD)=CRD(1:NPLREAD)
                ELSE
                   READ(12,REC=IREC,ERR=240) (CW2(I),I=1,NPLREAD)
                ENDIF
                CW1=0
                CALL REPAD_WITH_INDEX_ARRAY( MALLOC, IND, INDI, CW1(NPL+1), CW2)
              ENDIF

              IF (WDES%COMM_KINTER%NCPU.GT.1) THEN
                 CALL M_bcast_z(WDES%COMM_KINTER,CW1,SIZE(CW1))
              END IF

              CALL DIS_PW_BAND(WDES1, J+NBANDF, CW1, W%CPTWFP(1,1,K,1))
           ENDDO

        ENDIF

        DEALLOCATE(CW1,CW2,CRD,IND,INDI)

        ENDDO kpoints2

        IF (NKPTSF > WDES%NKPTS) THEN
           IREC=IREC+(NKPTSF-WDES%NKPTS)*(NBANDF+1)
        ENDIF

        ENDDO spin2

! copy eigenvalues and weights to all nodes
        NCOMM=WDES%NB_TOT*WDES%NKPTS
        CALL M_bcast_z(WDES%COMM, W%CELTOT(1,1,1),NCOMM )
        CALL M_bcast_d(WDES%COMM, W%FERTOT(1,1,1),NCOMM )

        IF (IU0 >= 0) WRITE(IU0,*) 'the WAVECAR file was read successfully'
        RETURN
      ENDIF

      MSG = 'ERROR: while reading WAVECAR, file is incompatible '
      IF (WDES%ENMAX /= ENMAXF) &
         MSG = MSG // 'the energy cutoff has changed (new,old) ' // str(WDES%ENMAX) // ' ' // str(ENMAXF)
      IF (NBANDF /= WDES%NB_TOT) &
         MSG = MSG // 'the number of bands has changed (new,old) ' // str(WDES%NB_TOT) // ' ' // str(NBANDF)
      IF (NKPTSF /= WDES%NKPTS) &
         MSG = MSG // 'the number of k-points has changed (new,old) ' //      str(WDES%NKPTS) // ' ' // str(NKPTSF)
      CALL vtutor%error(MSG)
!=======================================================================
! can not do anything with WAVECAR
! hard stop, pull all breaks
!=======================================================================
  200 CONTINUE

      CALL vtutor%error("ERROR: while reading WAVECAR, header is corrupt")

  220 CONTINUE

      CALL vtutor%error("ERROR: while reading WAVECAR, plane wave coefficients changed " // str(NPL) &
         // " " // str(NPLREAD))

  230 CONTINUE

      CALL vtutor%error("ERROR: while reading eigenvalues from WAVECAR " // str(K) // " " // str(ISP))

  240 CONTINUE

      CALL vtutor%error("ERROR: while reading plane wave coef. from WAVECAR " // str(K) // " " // &
         str(ISP) // " " // str(J))

      END SUBROUTINE

!*************************SUBROUTINE INWAV_ALTERNATIVE *****************
!
!>   read wavefunctions from alternative file
!>
!>   read alternative WAVECAR file
!
!***********************************************************************

      SUBROUTINE INWAV_ALTERNATIVE(IO, WDES, W, GRID, LATT_CUR, LREAD, EXT )

      USE base
      USE wave_high
      USE lattice
      USE tutor, ONLY: vtutor, isAlert, AuxWAVECAR

      IMPLICIT NONE

      TYPE (in_struct)  IO
      TYPE (wavedes)  WDES
      TYPE (wavespin) W
      TYPE (grid_3d)  GRID
      TYPE (latt)     LATT_CUR
      LOGICAL         LREAD
      CHARACTER (LEN=3)  EXT
! local
      TYPE (latt)     LATT_INI
      INTEGER         ISTART
      REAL(q)         ENMAXI
      REAL (q)        EFERMI

      LREAD = .FALSE.

      ISTART=1
      CALL OPENWAV(IO, WDES%COMM, EXT)
      CALL INWAV_HEAD(WDES, LATT_INI, LATT_CUR, ENMAXI, ISTART, IO%IU0)

      IF (ISTART/=0) THEN
         CALL ALLOCW( WDES, W)
         CALL INWAV_FAST(IO, WDES, W, GRID, LATT_CUR, LATT_INI, ISTART, EFERMI)
         LREAD = .TRUE.
         CALL vtutor%write(isAlert, AuxWAVECAR)
      ENDIF

      CALL CLOSEWAV

      END SUBROUTINE INWAV_ALTERNATIVE

!*************************SUBROUTINE OUTWAV    ************************
!
!>   write wavefunctions to file
!>
!>   the first version writes a VASP 4 file
!>   the second version writes the newer VASP 5 file format
!>   which requires less storage if many bands are calculated
!
!***********************************************************************

      SUBROUTINE OUTWAV_4(IO, WDES, W, LATT_INI, EXT )
      USE base
      USE wave
      USE lattice
      USE main_mpi

      TYPE (in_struct)   IO
      TYPE (wavedes)  WDES
      TYPE (wavespin) W
      TYPE (latt)     LATT_INI
      CHARACTER (LEN=3), OPTIONAL :: EXT
! local work arrays
      TYPE (wavedes1) WDES1
      COMPLEX(q),ALLOCATABLE :: CPTWFP(:),EIG(:)
# 702

      COMPLEX(qs),  ALLOCATABLE :: CRD(:)

! local
      INTEGER NODE_ME, IONODE, IU0, IREC, IRECLW
      INTEGER NPL_TOT, NPL, I, J, K, ISP
      REAL(q) RDUM, RISPIN, RTAG, RNB_TOT, RNKPTS, RNPL
      LOGICAL, POINTER :: LCONJG(:)=>null()
      INTEGER, POINTER :: INDEX(:)=>null()


      NODE_ME=0
      IONODE=0

      NODE_ME=WDES%COMM%NODE_ME
      IONODE =WDES%COMM%IONODE


!
! report to  stdout
!
      IF (IO%IU6>=0) WRITE(IO%IU6,*)'writing wavefunctions'
      IF (IO%IU0>=0) WRITE(IO%IU0,*)'writing wavefunctions'
!
! open file
!
      IU0   =IO%IU0
      NPL_TOT = MAXVAL(WDES%NPLWKP_TOT)
# 732

      IO%IRECLW=MAX(MAX((NPL_TOT+1)/2,6),((WDES%NB_TOT*3+1)/2+2))*IO%ICMPLX

      IRECLW=IO%IRECLW

! open file
      IF (NODE_ME==IONODE) THEN
      IF (PRESENT(EXT)) THEN
         OPEN(UNIT=12,FILE=DIR_APP(1:DIR_LEN)//'WAVECAR.'//EXT,ACCESS='DIRECT', &
              FORM='UNFORMATTED',STATUS='REPLACE',RECL=IO%IRECLW)
      ELSE
         OPEN(UNIT=12,FILE=DIR_APP(1:DIR_LEN)//'WAVECAR',ACCESS='DIRECT', &
              FORM='UNFORMATTED',STATUS='REPLACE',RECL=IO%IRECLW)
      ENDIF
      ENDIF

      ALLOCATE(CPTWFP(NPL_TOT),CRD(NPL_TOT),EIG(WDES%NB_TOT))

      RISPIN=WDES%ISPIN
      RDUM  =IRECLW

# 755

      RTAG  =45200

      IF (NODE_ME==IONODE) WRITE(12,REC=1) RDUM,RISPIN,RTAG

      IREC=2
      IF (NODE_ME==IONODE) THEN
! in order to increase exchangeability of WAVECAR files across IEEE platforms
! avoid INTEGERS on output, write REAL(q) items instead (same below with RNPL)
      RNKPTS =WDES%NKPTS
      RNB_TOT=WDES%NB_TOT
      WRITE(12,REC=2) RNKPTS,RNB_TOT,WDES%ENMAX,((LATT_INI%A(I,J),I=1,3),J=1,3)
      ENDIF
!
! loop over spin, kpoints, bands
!
      spin: DO ISP=1,WDES%ISPIN
      kpoints: DO K=1,WDES%NKPTS


       IF (WDES%COMM_KINTER%NCPU.GT.1) THEN
          IF (WDES%COMM_KINTER%NODE_ME.NE.1) THEN
             IF (MOD(K-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) THEN
                CYCLE
             ELSE
                CALL M_send_z(WDES%COMM_KINTER, 1, W%CPTWFP(1,1,K,ISP),SIZE(W%CPTWFP,1)*SIZE(W%CPTWFP,2))
             END IF
          ELSE
             IF (MOD(K-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) THEN
                CALL M_recv_z(WDES%COMM_KINTER, MOD(K-1,WDES%COMM_KINTER%NCPU)+1, W%CPTWFP(1,1,K,ISP),SIZE(W%CPTWFP,1)*SIZE(W%CPTWFP,2))
             ENDIF
          ENDIF
       END IF

       CALL SETWDES(WDES,WDES1,K)
# 792

       NPL=WDES%NPLWKP_TOT(K)
       RNPL=NPL
       IREC=IREC+1
! write number of plane waves, k-point coordinates and all eigenvalues and
! occupation numbers for current k-point K
       IF (NODE_ME==IONODE) THEN
! write eigenvalues in real format
         DO J=1,WDES%NB_TOT
           EIG(J)=REAL(W%CELTOT(J,K,ISP),KIND=q)
         ENDDO
         WRITE(12,REC=IREC) RNPL,WDES%VKPT(1,K),WDES%VKPT(2,K), &
                       WDES%VKPT(3,K),(EIG(J),W%FERTOT(J,K,ISP),J=1,WDES%NB_TOT)
       ENDIF

       DO J=1,WDES%NB_TOT
# 810

         CALL MRG_PW_BAND(WDES1, J, CPTWFP, W%CPTWFP(1,1,K,ISP))


         IREC=IREC+1
         IF (NODE_ME==IONODE) CRD(1:NPL)=CPTWFP(1:NPL)
         IF (NODE_ME==IONODE) WRITE(12,REC=IREC) (CRD(I),I=1,NPL)
       ENDDO

# 821

      ENDDO kpoints
      ENDDO spin

      DEALLOCATE(CPTWFP,CRD,EIG)

      IF (NODE_ME==IONODE) CALL CLOSEWAV

      RETURN
      END SUBROUTINE


!*************************SUBROUTINE OUTWAV    ************************
!
!>   write wavefunctions to file
!>
!>   this version uses smaller record lenght
!>   and if required eigenvalues/ occupancies will span
!>   several records
!
!***********************************************************************

      SUBROUTINE OUTWAV(IO, WDES, W, LATT_INI, EFERMI, EXT, NBANDS_DUMP)
      USE base
      USE wave
      USE lattice
      USE main_mpi
      use tutor, only: VTUTOR

      TYPE (in_struct)   IO
      TYPE (wavedes)  WDES
      TYPE (wavespin) W
      TYPE (latt)     LATT_INI
      CHARACTER (LEN=3), OPTIONAL :: EXT
      INTEGER,INTENT(IN),  OPTIONAL :: NBANDS_DUMP
      REAL (q) EFERMI
! local work arrays
      TYPE (wavedes1) WDES1
      COMPLEX(q),ALLOCATABLE :: CPTWFP(:)
# 862

      COMPLEX(qs),  ALLOCATABLE :: CRD(:)

! local
      INTEGER IU0, NODE_ME, IONODE, IREC, IRECLW_OLD
      INTEGER NPL, NPL_TOT, I, J, K, ISP
      REAL(q) RDUM, RISPIN, RNKPTS, RNB_TOT, RTAG, RNPL
      LOGICAL, POINTER :: LCONJG(:)=>null()
      INTEGER, POINTER :: INDEX(:)=>null()
      REAL(q), ALLOCATABLE :: OUTBUF(:)
      INTEGER  NBANDS_OLD 

      NODE_ME=0
      IONODE=0

      NODE_ME=WDES%COMM%NODE_ME
      IONODE =WDES%COMM%IONODE


      NBANDS_OLD = WDES%NB_TOT 
      IF ( PRESENT( NBANDS_DUMP ) ) THEN
         IF( NBANDS_DUMP > 1 ) THEN
            WDES%NB_TOT = NBANDS_DUMP
         ELSE
            CALL VTUTOR%BUG(" silly number of bands passed "//STR(NBANDS_DUMP),"fileio.F", 886 )
         ENDIF
      ENDIF
!
! report to  stdout
!
      IF (IO%IU6>=0) WRITE(IO%IU6,*)'writing wavefunctions'
      IF (IO%IU0>=0) WRITE(IO%IU0,*)'writing wavefunctions'
!
! open file
!
      IU0   =IO%IU0
      NPL_TOT = MAXVAL(WDES%NPLWKP_TOT)
# 902

      IO%IRECLW=MAX((NPL_TOT+1)/2,7)*IO%ICMPLX
      IRECLW_OLD=MAX(MAX((NPL_TOT+1)/2,6),((WDES%NB_TOT*3+1)/2+2))*IO%ICMPLX


! open file
      IF (NODE_ME==IONODE) THEN
      IF (PRESENT(EXT)) THEN
         OPEN(UNIT=12,FILE=DIR_APP(1:DIR_LEN)//'WAVECAR.'//EXT,ACCESS='DIRECT', &
              FORM='UNFORMATTED',STATUS='REPLACE',RECL=IO%IRECLW)
      ELSE
         OPEN(UNIT=12,FILE=DIR_APP(1:DIR_LEN)//'WAVECAR',ACCESS='DIRECT', &
              FORM='UNFORMATTED',STATUS='REPLACE',RECL=IO%IRECLW)
      ENDIF
      ENDIF

      ALLOCATE(CPTWFP(NPL_TOT),CRD(NPL_TOT))

      RISPIN=WDES%ISPIN
      RDUM  =IO%IRECLW

# 925

      RTAG  =53300

! if the record lenght is identical to previous vasp version we can claim
! to write the old file format
      IF (IRECLW_OLD==IO%IRECLW) THEN
# 933

         RTAG  =45200

      ENDIF

      IF (NODE_ME==IONODE) WRITE(12,REC=1) RDUM,RISPIN,RTAG

      IREC=2
      IF (NODE_ME==IONODE) THEN
! in order to increase exchangeability of WAVECAR files across IEEE platforms
! avoid INTEGERS on output, write REAL(q) items instead (same below with RNPL)
      RNKPTS =WDES%NKPTS
      RNB_TOT=WDES%NB_TOT
      WRITE(12,REC=2) RNKPTS,RNB_TOT,WDES%ENMAX,((LATT_INI%A(I,J),I=1,3),J=1,3),EFERMI
      ENDIF
!
! loop over spin, kpoints, bands
!
      spin: DO ISP=1,WDES%ISPIN
      kpoints: DO K=1,WDES%NKPTS


       IF (WDES%COMM_KINTER%NCPU > 1) &
          CALL M_bcast_z8_from(WDES%COMM_KINTER, W%CPTWFP(1,1,K,ISP),SIZE(W%CPTWFP,1,KIND=qi8)*SIZE(W%CPTWFP,2,KIND=qi8),MOD(K-1,WDES%COMM_KINTER%NCPU)+1)

       CALL SETWDES(WDES,WDES1,K)
# 961

       NPL=WDES%NPLWKP_TOT(K)
       RNPL=NPL
! write number of plane waves, k-point coordinates and all eigenvalues and
! occupation numbers for current k-point K
! write eigenvalues in real format
       IF (NODE_ME==IONODE) THEN
       ALLOCATE(OUTBUF(4+WDES%NB_TOT*3))

       OUTBUF(1)=RNPL
       OUTBUF(2)=WDES%VKPT(1,K)
       OUTBUF(3)=WDES%VKPT(2,K)
       OUTBUF(4)=WDES%VKPT(3,K)
       DO I=1,WDES%NB_TOT
          OUTBUF(5+(I-1)*3) =REAL(W%CELTOT(I,K,ISP),q)
          OUTBUF(6+(I-1)*3) =AIMAG(W%CELTOT(I,K,ISP))
          OUTBUF(7+(I-1)*3) =W%FERTOT(I,K,ISP)
       ENDDO
       CALL WRITE_FROM_BUF(12, IREC, IO%IRECLW/IO%ICMPLX*2, OUTBUF)
       DEALLOCATE(OUTBUF)
       ENDIF

       DO J=1,WDES%NB_TOT
# 986

         CALL MRG_PW_BAND(WDES1, J, CPTWFP, W%CPTWFP(1,1,K,ISP))


         IREC=IREC+1
         IF (NODE_ME==IONODE) CRD(1:NPL)=CPTWFP(1:NPL)
         IF (NODE_ME==IONODE) WRITE(12,REC=IREC) (CRD(I),I=1,NPL)
       ENDDO

# 997

      ENDDO kpoints
      ENDDO spin

      DEALLOCATE(CPTWFP,CRD)

      IF (NODE_ME==IONODE) CALL CLOSEWAV

! restore bands to original number
      WDES%NB_TOT = NBANDS_OLD 
      RETURN
      END SUBROUTINE


      SUBROUTINE WRITE_FROM_BUF(IU, IREC, IRECL_REAL, OUTBUF)
      INTEGER :: IU         ! unit for writing
      INTEGER :: IREC       ! record to write
      INTEGER :: IRECL_REAL ! record lenght in real(q) words
      REAL(q) :: OUTBUF(:)  ! buffer to written to file

! local
      INTEGER :: I

      DO I=1,(SIZE(OUTBUF)+IRECL_REAL-1)/IRECL_REAL
         IREC=IREC+1
         WRITE(12,REC=IREC) OUTBUF((I-1)*IRECL_REAL+1 : MIN(I*IRECL_REAL, SIZE(OUTBUF)))
      ENDDO

      END SUBROUTINE

      SUBROUTINE READ_TO_BUF(IU, IREC, IRECL_REAL, INBUF, IFAIL)
      INTEGER :: IU         ! unit for reading
      INTEGER :: IREC       ! record to read from
      INTEGER :: IRECL_REAL ! record lenght in real(q) words
      REAL(q) :: INBUF(:)   ! buffer to read to
      INTEGER :: IFAIL      ! error condition

! local
      INTEGER :: I

      IFAIL=0
      DO I=1,(SIZE(INBUF)+IRECL_REAL-1)/IRECL_REAL
         IREC=IREC+1
         READ(12,REC=IREC, IOSTAT=IFAIL) INBUF((I-1)*IRECL_REAL+1 : MIN(I*IRECL_REAL, SIZE(INBUF)))
         IF (IFAIL/=0) RETURN
      ENDDO

      END SUBROUTINE

!*************************SUBROUTINE OUTCHG ****************************
!
!>   write chargedensity to a specified unit
!>
!>   HEADER is currently created in the main program
!>   charge density CHTOT must be given in reciprocal space
!
!***********************************************************************

      SUBROUTINE OUTCHG(GRIDC, IU, LLONG,CHTOT)
      USE prec
      USE mpimy
      USE mgrid
      USE string, ONLY: str
      USE tutor, ONLY: vtutor

      TYPE (grid_3d) GRIDC

      COMPLEX(q)  CHTOT(GRIDC%RC%NP)
      LOGICAL LLONG        ! long or short format
! work arrays
      COMPLEX(q),  ALLOCATABLE ::  CWORK(:)
      REAL(q),ALLOCATABLE ::  WORK (:)
      INTEGER NALLOC,NZ, NWRITE, NWRITTEN
      CHARACTER (40) FORM
      INTEGER ISTAT, NODE_ME, IONODE, IU, N, IERR
# 1075


      NALLOC=GRIDC%NGX*GRIDC%NGY

      ALLOCATE(CWORK(GRIDC%MPLWV),WORK(NALLOC),STAT=ISTAT)
      IF (ISTAT>0) RETURN ! can not write the charge immediate exit

      NODE_ME=0
      IONODE =0

      NODE_ME=GRIDC%COMM%NODE_ME
      IONODE =GRIDC%COMM%IONODE


      IF (GRIDC%NPLWV/= GRIDC%NGX*GRIDC%NGY*GRIDC%NGZ) THEN
        CALL vtutor%bug("internal ERROR: OUTCHG NPLWV is not compatibel with NGX,NGY,NGZ \n " // &
           str(GRIDC%NPLWV) // " " // str(GRIDC%NGX) // " " // str(GRIDC%NGY) // " " // str(GRIDC%NGZ), "fileio.F", 1091)
      ENDIF
! transfer to CWORK and FFT
      CALL RC_ADD(CHTOT,1.0_q,CHTOT,0.0_q,CWORK,GRIDC)
! transform to real space representation
      CALL FFT3D(CWORK,GRIDC,1)

      IF (LLONG) THEN
        FORM='(1X,E17.11)'
        NWRITE=5
      ELSE
        FORM='(1X,G11.5)'
        NWRITE=10
      ENDIF

      IF (NODE_ME==IONODE) WRITE(IU,'(3I5)') GRIDC%NGX,GRIDC%NGY,GRIDC%NGZ

      NWRITTEN=0
      DO NZ=1,GRIDC%NGZ
         CALL MRG_GRID_RL_PLANE(GRIDC, WORK, CWORK, NZ)
         IF (NODE_ME==IONODE) THEN

         DO N=1,NALLOC
            NWRITTEN=NWRITTEN+1
            IF ( MOD(NWRITTEN,NWRITE)==0 ) THEN
               WRITE(IU,FORM) WORK(N)
            ELSE
               WRITE(IU,FORM,ADVANCE='NO') WORK(N)
            ENDIF
         ENDDO
# 1141

         ENDIF
      ENDDO
      IF ( MOD(NWRITTEN,NWRITE)/=0 ) WRITE(IU,*)' '

      DEALLOCATE(CWORK,WORK)

      RETURN
      END SUBROUTINE

!*************************SUBROUTINE OUTCHG_RGRID **********************
!
!>   write chargedensity to a specified unit
!>
!>   simple copy from OUTCHG, however
!>   no FFT to real space is 1._q, so CHTOT must be provided on
!>   real-space grid
!
!***********************************************************************

      SUBROUTINE OUTCHG_RGRID(GRID, IU, LLONG,SV)
         USE poscar_struct_def
         USE prec
         USE mpimy
         USE mgrid
         USE string, ONLY: str
         USE tutor, ONLY: vtutor

         TYPE (grid_3d) GRID

         COMPLEX(q)  SV(GRID%MPLWV)
         LOGICAL LLONG        ! long or short format
! work arrays
         REAL(q),ALLOCATABLE ::  WORK (:)
         INTEGER NALLOC,NZ, NWRITE, NWRITTEN
         CHARACTER (40) FORM
         INTEGER ISTAT, NODE_ME, IONODE, IU, N, IERR

         NALLOC=GRID%NGX*GRID%NGY

         ALLOCATE(WORK(NALLOC),STAT=ISTAT)
         IF (ISTAT>0) RETURN ! can not write the charge immediate exit

         NODE_ME=0
         IONODE =0

         NODE_ME=GRID%COMM%NODE_ME
         IONODE =GRID%COMM%IONODE


         IF (GRID%NPLWV/= GRID%NGX*GRID%NGY*GRID%NGZ) THEN
           CALL vtutor%bug("OUTCHG_RGRID NPLWV is not compatibel with NGX,NGY,NGZ \n " // str(GRID%NPLWV) // " " // str(GRID%NGX) // " " // str(GRID%NGY) // " " // str(GRID%NGZ), "fileio.F", 1192 )
         ENDIF
! original subroutine stores CHTOT (complex) to SV (real)
! and transforms SV to real-space
! note that SV has 2 times the dimension of CHTOT
! so that Re and Im part of CHTOT fits into SV
! this requires that RC_ADD to reside outside of a module
! for the fortran compiler to understand this
!        ! transfer to SV and FFT
!        CALL RC_ADD(CHTOT,1.0_q,CHTOT,0.0_q,SV,GRID)
!        ! transform to real space representation
!        CALL FFT3D(SV,GRID,1)

         IF (LLONG) THEN
           FORM='(1(1X,E17.11))'
           NWRITE=5
         ELSE
           FORM='(1(1X,G11.5))'
           NWRITE=10
         ENDIF

         IF (NODE_ME==IONODE) WRITE(IU,'(3I5)') GRID%NGX,GRID%NGY,GRID%NGZ

         NWRITTEN=0
         DO NZ=1,GRID%NGZ
            CALL MRG_GRID_RL_PLANE(GRID, WORK, SV, NZ)
            IF (NODE_ME==IONODE) THEN
            DO N=1,NALLOC
               NWRITTEN=NWRITTEN+1
               IF ( MOD(NWRITTEN,NWRITE)==0 ) THEN
                  WRITE(IU,FORM) WORK(N)
               ELSE
                  WRITE(IU,FORM,ADVANCE='NO') WORK(N)
               ENDIF
            ENDDO
            ENDIF
         ENDDO
         IF ( MOD(NWRITTEN,NWRITE)/=0 ) WRITE(IU,*)' '

         DEALLOCATE(WORK)

         RETURN
   END SUBROUTINE OUTCHG_RGRID


!*************************SUBROUTINE  INCHG ****************************
!>   read chargedensity form a specified unit
!>
!>   HEADER must have been read previously
!***********************************************************************

      SUBROUTINE INCHG(IU,CHTOT,GRID,IERR)
      USE prec
      USE charge
      USE mpimy
      USE mgrid
      USE string, ONLY: str
      USE tutor, ONLY: vtutor

      TYPE (grid_3d) GRID
      COMPLEX(q)   CHTOT(GRID%RC%NP)
! work array
      REAL(q),ALLOCATABLE ::  WORK(:)
      INTEGER NALLOC,NZ, NWRITE, NWRITTEN
      INTEGER ISTAT, NODE_ME, IONODE, IU, IERR
      INTEGER NGXFIL, NGYFIL, NGZFIL, N

      NALLOC=GRID%NGX*GRID%NGY
      NODE_ME=0
      IONODE =0

      NODE_ME=GRID%COMM%NODE_ME
      IONODE =GRID%COMM%IONODE


      ALLOCATE(WORK(NALLOC),STAT=ISTAT)
      IF (ISTAT>0) THEN
         IERR=2
         RETURN                 ! can not read potential, immediate exit
      ENDIF

      IERR=0
      IF (GRID%NPLWV/= GRID%NGX*GRID%NGY*GRID%NGZ) THEN
        CALL vtutor%bug("internal ERROR: INCHG GRID%NPLWV is not compatible with GRID%NGX,GRID%NGY,&
           &NGZC \n " // str(GRID%NPLWV) // " " // str(GRID%NGX) // " " // str(GRID%NGY) // " " // &
           str(GRID%NGZ), "fileio.F", 1277)
      ENDIF

      IERR=0
      NGXFIL=0; NGYFIL=0; NGZFIL=0


      IF (NODE_ME==IONODE) THEN

        READ(IU,*,ERR=120,END=120) NGXFIL,NGYFIL,NGZFIL
 120    CONTINUE
        IF (NGXFIL==0 .AND. NGYFIL==0 .AND. NGZFIL==0) THEN
           IERR=2
        ELSE IF (NGXFIL/=GRID%NGX .OR. NGYFIL/=GRID%NGY .OR. NGZFIL/=GRID%NGZ ) THEN
          IERR=1
        ENDIF

      ENDIF


      CALL M_bcast_i( GRID%COMM, IERR, 1)
      IF (IERR /=0 ) GOTO 200

      NWRITE=5
      NWRITTEN=0
      DO NZ=1,GRID%NGZ
         IF (NODE_ME==IONODE) THEN
         DO N=1,NALLOC
            NWRITTEN=NWRITTEN+1
! non advancing read is not the same on all machine
! this is the only version that works on all tested platforms
            IF ( MOD(NWRITTEN,NWRITE)==0 .OR. (N==NALLOC .AND. NZ==GRID%NGZ)) THEN
               READ(IU,'(1(1X,E17.11))',ERR=100,END=100) WORK(N)
            ELSE
               READ(IU,'(1(1X,E17.11))',ADVANCE='NO',ERR=100,END=100) WORK(N)
            ENDIF
         ENDDO
         ENDIF
         CALL DIS_GRID_RL_PLANE(GRID, WORK, CHTOT, .TRUE., NZ )
      ENDDO

      CALL FFT_RC_SCALE(CHTOT,CHTOT,GRID)
      CALL SETUNB_COMPAT(CHTOT,GRID)

  200 CONTINUE
      DEALLOCATE(WORK)
      RETURN

  100 CONTINUE
      IERR=2
      DEALLOCATE(WORK)
      RETURN
      END SUBROUTINE

!*************************SUBROUTINE READNI ****************************
!>   reads total number of ions from a line consisting of numbers
!>   for several species
!***********************************************************************

      SUBROUTINE READNI(ZPARSE,NIONSF)
      USE prec
      USE tutor, ONLY: vtutor
      CHARACTER (80) ZPARSE
      CHARACTER (80) ZITEM
      CHARACTER (80) ZWORK
      CHARACTER (15) ZFORM
      CHARACTER (1)  ZCHAR
      INTEGER IDATA, NDATA, NIONSF, NI, LFORM
      INTEGER, EXTERNAL :: NITEMS

!  get number of data items
      NDATA=NITEMS(ZPARSE,ZWORK,.TRUE.,'I')

!  parse and read list:
      NIONSF=0

      DO 100 IDATA=1,NDATA
         CALL SUBWRD(ZPARSE,ZITEM,IDATA,1)
         CALL CHKINT(ZITEM,ZWORK,ZCHAR,ZFORM)
         IF (ZCHAR=='Y') THEN
            ZWORK='('//ZFORM//')'
            CALL STRIP(ZWORK,LFORM,'A')
            READ(ZITEM,ZWORK(1:LFORM)) NI
         ELSE
            CALL vtutor%error("Fatal error in READNI: Invalid data found ...")
         ENDIF

         NIONSF=NIONSF+NI

  100 CONTINUE

      RETURN
      END SUBROUTINE

!*************************SUBROUTINE READCH ****************************
!
!>  This routine reads in a chargedensity from file IU
!>
!>  When the chargedensity of the file is incompatible
!>  i.e. NGX,Y,Z differ from that used in the programm a warning is
!>  reported and ICHARG is set to 0
!>  if the magnetisation charge density could not be read in
!>  ICHARG is set to -1
!>  on entry to the routine:
!>     CHTOT  may be uninitialised
!>     RHOLM  must be already correctly initialised according to MAGMOM
!>  on exit:
!>     CHTOT  as read from file
!>     RHOLM  as read from file
!>  the routine is complicated by several special cases:
!>  o if the magnetisation charge density could not be read in
!>     ICHARG is set to -1, CHTOT(:,2:4)=0
!>  o non collinear case: if only mz could be read
!>     the x and y components are set to 0 in CHTOT and RHOLM
!
!***********************************************************************

      SUBROUTINE READCH(GRIDC, LOVERL, T_INFO, CHTOT, RHOLM, ICHARG, ISPIN, &
           LATT_CUR, P, CSTRF, IU, IU0)
      USE prec
      USE lattice
      USE poscar
      USE mgrid
      USE pseudo
      USE charge
      USE pawm

      TYPE (grid_3d)     GRIDC
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P (T_INFO%NTYP)
      TYPE (latt)        LATT_CUR
      COMPLEX(q)   CSTRF(GRIDC%MPLWV,T_INFO%NTYP)
      INTEGER ISPIN, IU, IU0, ICHARG
      REAL(q)         :: RHOLM(:,:)
! local work arrays
      TYPE (type_info)   T_INFO_OLD
      REAL(q)         :: TMP(T_INFO%NIONS)
      COMPLEX(q) CHTOT(GRIDC%MPLWV,ISPIN)
      CHARACTER (40) TITEL
      CHARACTER (80) ZPARSE
      CHARACTER (1)  CHARAC
      COMPLEX(q), ALLOCATABLE :: CD(:),CSTRF_OLD(:,:)
      INTEGER I, NIONSF, ISP, N
      INTEGER NODE_ME, IONODE, IERR
      LOGICAL LOVERL

      NODE_ME=0
      IONODE =0

      NODE_ME=GRIDC%COMM%NODE_ME
      IONODE =GRIDC%COMM%IONODE


!=======================================================================
! read in old ionic-positions and header
!=======================================================================
      T_INFO_OLD=T_INFO

      NULLIFY(T_INFO_OLD%POSION)
      ALLOCATE(T_INFO_OLD%POSION(3,T_INFO_OLD%NIONS))

      READ(IU,*,ERR=1000,END=1000) TITEL
      DO I=1,4
        READ(IU,*,ERR=1000,END=1000)
      ENDDO
      READ(IU,'(A80)') ZPARSE
      READ(ZPARSE,*) CHARAC
      IF (.NOT.(CHARAC>='0' .AND. CHARAC<='9')) READ(IU,'(A80)') ZPARSE

      CALL READNI(ZPARSE,NIONSF)
      IF (NIONSF/=T_INFO%NIONS) THEN
        IF (IU0>=0) &
        WRITE(IU0,*) 'WARNING: number of atoms are different on CHGCAR file'
        ICHARG=0
        RETURN
      ENDIF
      READ(IU,*,ERR=1000,END=1000)

      DO  I=1,NIONSF
        READ(IU,*) T_INFO_OLD%POSION(:,I)
      ENDDO
!=======================================================================
! read in the charge-density
!=======================================================================
      CALL INCHG(IU,CHTOT(1,1),GRIDC,IERR)
      IF (IERR==1) THEN
        IF (IU0>=0) &
        WRITE(*,*) 'WARNING: dimensions on CHGCAR file are different'
        ICHARG=0
        RETURN
      ENDIF
      IF (IERR==2) THEN
        IF (IU0>=0) &
        WRITE(IU0,*) 'WARNING: chargedensity file is incomplete'
        ICHARG=0
        RETURN
      ENDIF
!=======================================================================
! now subtract charge density according to old position and add
! that (1._q,0._q) according to new positions
!=======================================================================
      ALLOCATE(CD(GRIDC%MPLWV),CSTRF_OLD(GRIDC%MPLWV,T_INFO%NTYP))

!---- subtract atomic charge-denisty for old positions
      CALL STUFAK(GRIDC,T_INFO_OLD,CSTRF_OLD)
      CALL RHOATO_WORK(.FALSE.,.FALSE.,GRIDC,T_INFO_OLD,LATT_CUR%B,P,CSTRF_OLD,CD)
      DO N=1,GRIDC%RC%NP
        CHTOT(N,1)= CHTOT(N,1)- CD(N)
      ENDDO

!---- add atomic charge-denisty for new positions
      CALL RHOATO_WORK(.FALSE.,.FALSE.,GRIDC,T_INFO,LATT_CUR%B,P,CSTRF,CD)
      DO N=1,GRIDC%RC%NP
        CHTOT (N,1)= CHTOT(N,1)+ CD(N)
      ENDDO

      DEALLOCATE(CD,CSTRF_OLD,T_INFO_OLD%POSION)

! just in case initialize the magnetization density to 0
      DO ISP=2,ISPIN
        CALL RC_ADD(CHTOT(1,ISP),0.0_q,CHTOT(1,ISP),0.0_q,CHTOT(1,ISP),GRIDC)
      ENDDO

      CALL RD_RHO_PAW(P, T_INFO, LOVERL, RHOLM(:,1), GRIDC%COMM, IU, IERR )
      IF (IERR/=0) THEN
        IF (IU0>=0) &
        WRITE(IU0,*) 'WARNING: PAW occupancies are missing on CHGCAR '
        ICHARG=0
        RETURN
      ENDIF
!=======================================================================
! set/read spin density
!=======================================================================
      IF (IU0>=0) &
        WRITE(IU0,*) 'charge-density read from file: ',TITEL

      DO ISP=2,ISPIN
! read in the spin-density
        IERR=0
! in vasp.4.4 the ATOMOM array was read from the CHGCAR file
! reading the magnetic moments from the file really makes no sense
! and spoiles non collinear calculations
        IF (NODE_ME==IONODE) READ(IU,*,IOSTAT=IERR) (TMP(I),I=1,T_INFO%NIONS)
        CALL M_bcast_i( GRIDC%COMM, IERR, 1)

        IF (IERR==0) CALL INCHG(IU,CHTOT(1,ISP),GRIDC,IERR)
        IF (IERR==0) CALL RD_RHO_PAW(P, T_INFO, LOVERL, RHOLM(:,ISP), GRIDC%COMM, IU, IERR )

! error occured
        IF (IERR/=0) THEN
           IF (ISP==2) THEN
! we could not read any entry in the magnetisation density return ICHARG=-1
! (magnetisation according to overlapping atoms)
              ICHARG=-1
              RETURN
           ELSE
! non collinear case we could read at least m_z(r)
! copy that to the right place and initialise everything else to 0
              CALL RC_ADD(CHTOT(1,2),1.0_q,CHTOT(1,2),0.0_q,CHTOT(1,4),GRIDC) ! copy
              CALL RC_ADD(CHTOT(1,2),0.0_q,CHTOT(1,2),0.0_q,CHTOT(1,2),GRIDC) ! =0
              CALL RC_ADD(CHTOT(1,3),0.0_q,CHTOT(1,3),0.0_q,CHTOT(1,3),GRIDC) ! =0

              RHOLM(:,4)=RHOLM(:,2)
              RHOLM(:,2)=0
              RHOLM(:,3)=0

              RETURN
           ENDIF
        ENDIF

        IF (IU0>=0) WRITE(IU0,'(A,I2)') ' magnetization density read from file',ISP-1

      ENDDO
      RETURN

      ICHARG=-1
      RETURN

 1000 CONTINUE
      IF (IU0>=0) &
       WRITE(IU0,*) 'WARNING: chargedensity file is incomplete'
      ICHARG=0
      RETURN
      END SUBROUTINE


!*************************SUBROUTINE READPOT ***************************
!
!>  This routine reads in a chargedensity from file IU
!>
!>  When the chargedensity of the file is incompatible
!>  i.e. NGX,Y,Z differ from that used in the programm a warning is
!>  reported and ICHARG is set to 0
!>  if the magnetisation charge density could not be read in
!>  ICHARG is set to -1
!>  on entry to the routine:
!>     CHTOT  may be uninitialised
!>     RHOLM  must be already correctly initialised according to MAGMOM
!>  on exit:
!>     CHTOT  as read from file
!>     RHOLM  as read from file
!>  the routine is complicated by several special cases:
!>  o if the magnetisation charge density could not be read in
!>     ICHARG is set to -1, CHTOT(:,2:4)=0
!>  o non collinear case: if only mz could be read
!>     the x and y components are set to 0 in CHTOT and RHOLM
!
!***********************************************************************

      SUBROUTINE READPOT(GRIDC, T_INFO, CVTOT, ICHARG, WDES, &
              LATT_CUR, P, LMDIM, CDIJ, N_MIX_PAW, IU, IU0)

      USE prec
      USE lattice
      USE poscar
      USE mgrid
      USE pseudo
      USE charge
      USE pawm
      USE wave

      TYPE (grid_3d)     GRIDC
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P (T_INFO%NTYP)
      TYPE (latt)        LATT_CUR
      TYPE (wavedes)     WDES
      INTEGER LMDIM
      COMPLEX(q)  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      INTEGER N_MIX_PAW
      INTEGER IU, IU0, ICHARG
! local work arrays
      REAL(q)         :: DLM(N_MIX_PAW,WDES%NCDIJ)
      TYPE (type_info)   T_INFO_OLD
      REAL(q)         :: TMP(T_INFO%NIONS)
      COMPLEX(q) CVTOT(GRIDC%MPLWV,WDES%NCDIJ)
      CHARACTER (40) TITEL
      CHARACTER (80) ZPARSE
      CHARACTER (1)  CHARAC
      INTEGER I, NIONSF, N, ISP
      INTEGER NODE_ME, IONODE, IERR

      NODE_ME=0
      IONODE =0

      NODE_ME=GRIDC%COMM%NODE_ME
      IONODE =GRIDC%COMM%IONODE


!=======================================================================
! read in old ionic-positions and header
!=======================================================================
      T_INFO_OLD=T_INFO

      NULLIFY(T_INFO_OLD%POSION)
      ALLOCATE(T_INFO_OLD%POSION(3,T_INFO_OLD%NIONS))

      READ(IU,*,ERR=1000,END=1000) TITEL
      DO I=1,4
        READ(IU,*,ERR=1000,END=1000)
      ENDDO
      READ(IU,'(A80)') ZPARSE
      READ(ZPARSE,*) CHARAC
      IF (.NOT.(CHARAC>='0' .AND. CHARAC<='9')) READ(IU,'(A80)') ZPARSE

      CALL READNI(ZPARSE,NIONSF)
      IF (NIONSF/=T_INFO%NIONS) THEN
        IF (IU0>=0) &
        WRITE(IU0,*) 'WARNING: number of atoms are different on CHGCAR file'
        ICHARG=0
        RETURN
      ENDIF
      READ(IU,*,ERR=1000,END=1000)

      DO  I=1,NIONSF
        READ(IU,*) T_INFO_OLD%POSION(:,I)
      ENDDO
!=======================================================================
! read in the charge-density
!=======================================================================
      CALL INCHG(IU,CVTOT(1,1),GRIDC,IERR)
      CALL FFT3D(CVTOT(1,1),GRIDC,1)

      IF (IERR==1) THEN
        IF (IU0>=0) &
        WRITE(*,*) 'WARNING: dimensions on LOCPOT file are different'
        ICHARG=0
        RETURN
      ENDIF
      IF (IERR==2) THEN
        IF (IU0>=0) &
        WRITE(IU0,*) 'WARNING: potential file is incomplete'
        ICHARG=0
        RETURN
      ENDIF

      DEALLOCATE(T_INFO_OLD%POSION)

! just in case initialize the magnetization density to 0
      DO ISP=2,WDES%NCDIJ
        CALL RC_ADD(CVTOT(1,ISP),0.0_q,CVTOT(1,ISP),0.0_q,CVTOT(1,ISP),GRIDC)
      ENDDO

      DLM(:,1)=0
      CALL RD_RHO_PAW(P, T_INFO, WDES%LOVERL, DLM(:,1), GRIDC%COMM, IU, IERR )
      IF (IERR/=0) THEN
        IF (IU0>=0) &
        WRITE(IU0,*) 'WARNING: PAW occupancies are missing on POT file'
        ICHARG=0
        RETURN
      ENDIF
!=======================================================================
! set/read spin density
!=======================================================================
      IF (IU0>=0) &
        WRITE(IU0,*) 'local potential read from file: ',TITEL

      DO ISP=2,WDES%NCDIJ
! read in the spin-density
        IERR=0
! in vasp.4.4 the ATOMOM array was read from the CHGCAR file
! reading the magnetic moments from the file really makes no sense
! and spoiles non collinear calculations
        IF (NODE_ME==IONODE) READ(IU,*,IOSTAT=IERR) (TMP(I),I=1,T_INFO%NIONS)
        CALL M_bcast_i( GRIDC%COMM, IERR, 1)

        IF (IERR==0) CALL INCHG(IU,CVTOT(1,ISP),GRIDC,IERR)
        CALL FFT3D(CVTOT(1,ISP),GRIDC,1)

        IF (IERR==0) CALL RD_RHO_PAW(P, T_INFO, WDES%LOVERL, DLM(:,ISP), GRIDC%COMM, IU, IERR )

! error occured
        IF (IERR/=0) THEN
           IF (ISP==2) THEN
! we could not read any entry in the magnetisation density return ICHARG=-1
! (magnetisation according to overlapping atoms)
              ICHARG=-1
              RETURN
           ELSE
! non collinear case we could read at least m_z(r)
! copy that to the right place and initialise everything else to 0
              CALL RC_ADD(CVTOT(1,2),1.0_q,CVTOT(1,2),0.0_q,CVTOT(1,4),GRIDC) ! copy
              CALL RC_ADD(CVTOT(1,2),0.0_q,CVTOT(1,2),0.0_q,CVTOT(1,2),GRIDC) ! =0
              CALL RC_ADD(CVTOT(1,3),0.0_q,CVTOT(1,3),0.0_q,CVTOT(1,3),GRIDC) ! =0

              DLM(:,4)=DLM(:,2)
              DLM(:,2)=0
              DLM(:,3)=0

              GOTO 2000
           ENDIF
        ENDIF

        IF (IU0>=0) WRITE(IU0,'(A,I2)') ' magnetization potential read from file',ISP-1

      ENDDO

 2000 CONTINUE
! (1._q,0._q) centre corrections are passed through mixer and added now
      CALL RETRIVE_RHO_PAW(WDES, P, T_INFO, WDES%LOVERL, LMDIM, CDIJ, DLM)

      RETURN

 1000 CONTINUE
      IF (IU0>=0) &
       WRITE(IU0,*) 'WARNING: potential file is incomplete'
      ICHARG=0
      RETURN
      END SUBROUTINE


!*************************SUBROUTINE OPENGAMMA *************************
!
! the followin subroutines open and close the file GAMMA for reading
! and writing
! furthermore they write and read the header of the file
!
!***********************************************************************

!> open the file GAMMA for reading and writing
      SUBROUTINE OPENGAMMA
        USE main_mpi
        OPEN(UNIT=IU_GAMMA  ,FILE=DIR_APP(1:DIR_LEN)//'GAMMA',STATUS='UNKNOWN')

      END SUBROUTINE OPENGAMMA


!> Close the file GAMMA
      SUBROUTINE CLOSEGAMMA

        CLOSE(UNIT=IU_GAMMA)

      END SUBROUTINE CLOSEGAMMA

!> Write the header to the file GAMMA
      SUBROUTINE WRITEGAMMA_HEAD( NKPTS , NBANDS, IO)
        USE base
        IMPLICIT NONE
        TYPE (in_struct)   IO
        INTEGER :: NKPTS   !< number of k-points
        INTEGER :: NBANDS  !< number of bands

        IF (IO%IU0>=0) THEN
           WRITE(IO%IU0,*) 'writing the density matrix to GAMMA'
        ENDIF
        IF (IO%IU6>=0) THEN
           WRITE(IU_GAMMA,'(2I10," ! number of k-points and bands")') NKPTS, NBANDS
        ENDIF

      END SUBROUTINE WRITEGAMMA_HEAD


!> Write the header to the file GAMMA
      SUBROUTINE WRITEGAMMA_HEAD_HAM( NKPTS , NBANDS, IO)
        USE base
        IMPLICIT NONE
        TYPE (in_struct)   IO
        INTEGER :: NKPTS   !< number of k-points
        INTEGER :: NBANDS  !< number of bands

        IF (IO%IU0>=0) THEN
           WRITE(IO%IU0,*) 'writing the density matrix to GAMMA'
        ENDIF
        IF (IO%IU6>=0) THEN
           WRITE(IU_GAMMA,'("Hamiltonian",2I10," ! number of k-points and bands")') NKPTS, NBANDS
        ENDIF

      END SUBROUTINE WRITEGAMMA_HEAD_HAM

!> Read the header from the file GAMMA
      SUBROUTINE READGAMMA_HEAD( IMODE, NKPTS , NBANDS, IO)
        USE base
        USE tutor, ONLY: vtutor, isError, Gamma
        IMPLICIT NONE
        TYPE (in_struct)   IO
        INTEGER :: NKPTS    !< number of k-points
        INTEGER :: NBANDS   !< number of bands
        CHARACTER (LEN=80) ::STRING ! header of Gamma file
        CHARACTER (LEN=1)  ::MODE   !  mode selector
        INTEGER :: IMODE    !< mode: IMODE=1 density matrix ; IMODE=2 Hamiltonian
! local
        INTEGER :: NKPTS_   ! number of k-points
        INTEGER :: NBANDS_  ! number of bands
        INTEGER :: IERR

! read first line
        READ(IU_GAMMA,'(A)') STRING
! try old format first
        IMODE=1
        READ(STRING,*, IOSTAT = IERR) NKPTS_, NBANDS_

        IF (IERR/=0) THEN
! try new format
           READ(STRING,*, IOSTAT = IERR) MODE, NKPTS_, NBANDS_
           IF (MODE=='H' .OR. MODE=='h') IMODE=2
        ENDIF

! if you do not care, (1._q,0._q) may simply supply -1 in the file
        IF (NKPTS_==-1)  NKPTS_=NKPTS
        IF (NBANDS_==-1) NBANDS_=NBANDS

        IF (IERR/= 0 .OR. NBANDS /= NBANDS_ .OR. NKPTS_ /= NKPTS) THEN
           CALL vtutor%write(isError, Gamma)
        ENDIF

        IF (IO%IU0>=0) THEN
          WRITE(IO%IU0,*) 'reading the density matrix from GAMMA'
        ENDIF

      END SUBROUTINE READGAMMA_HEAD


!*************************SUBROUTINE WRITEGAMMA *************************
!
!> write the change of the density matrix for (1._q,0._q) k-point
!>
!> the matrix MUST not be stored in a distributed fashion, and
!> only (1._q,0._q) core must call the routine
!
!***********************************************************************

      SUBROUTINE WRITEGAMMA( NK, NB_TOT, NBDIM, CHAM , LROT)
        USE prec
        IMPLICIT NONE
        INTEGER :: NK    !< k-point index
        INTEGER :: NB_TOT!< total number of bands
        INTEGER :: NBDIM !< leading dimension of CHAM
        COMPLEX(q) :: CHAM(NBDIM, NB_TOT)
!> (1._q,0._q) can also supply an approximate unitary rotation matrix
!> for instance obtained by Loewdin perturbation theory
        LOGICAL :: LROT
        INTEGER :: N1, N2
        INTEGER :: NB_START, NB_END

! change sign in the upper triangle (convert unitary matrix to density matrix)
        IF (LROT) THEN
           DO N1=1, NB_TOT
              DO N2=N1+1, NB_TOT
                 CHAM(N1, N2)= -CHAM(N1, N2)
              ENDDO
              CHAM(N1, N1)=0
           ENDDO
        ENDIF

! write starting band, and final band
        NB_START = 1
        NB_END = NB_TOT
        WRITE(IU_GAMMA, '(3I10," ! k-point index, start and end band")') NK, NB_START, NB_END

        DO N1=NB_START, NB_END
           WRITE(IU_GAMMA,'(16E15.7)') CHAM(N1, NB_START:NB_END)
        ENDDO

! change back to rotation matrix)
        IF (LROT) THEN
           DO N1=1, NB_TOT
              DO N2=N1+1, NB_TOT
                 CHAM(N1, N2)= -CHAM(N1, N2)
              ENDDO
              CHAM(N1, N1)=1
           ENDDO
        ENDIF

      END SUBROUTINE WRITEGAMMA


!> read the change of the density matrix for (1._q,0._q) k-point
!>
!> READGAMMA may be either called by all cores, or better
!> called only by (1._q,0._q) core with a subsequent broadcast of the
!> data to all cores
      SUBROUTINE READGAMMA( NK, NB_TOT, NBDIM, CHAM, IO )
        USE prec
        USE base
        USE tutor, ONLY: vtutor, isError, Gamma2
        IMPLICIT NONE
        TYPE (in_struct)   IO
        INTEGER :: NK    !< k-point index
        INTEGER :: NB_TOT!< total number of bands
        INTEGER :: NBDIM !< leading dimension of CHAM
        COMPLEX(q) :: CHAM(NBDIM, NB_TOT)
! local
        INTEGER :: NK_
        INTEGER :: N1, N2
        INTEGER :: NB_START, NB_END
        REAL(q) :: BUFFER(NBDIM*2)

! read starting band, and final band
        NB_START = 1
        NB_END = NB_TOT
        READ(IU_GAMMA, *) NK_, NB_START, NB_END

        NB_START=MAX(1, NB_START)
        NB_END  =MIN(NB_TOT, NB_END)
        IF (NK_==-1) NK_=NK

        IF (NK /= NK_) THEN
            CALL vtutor%write(isError, Gamma2)
        ENDIF

        CHAM=0
        DO N1=NB_START,NB_END
# 1940

! not to read into real buffer (since complex entries must be in parantheses)
           READ(IU_GAMMA,*) BUFFER(1:(NB_END-NB_START+1)*2)
! store as complex
           DO N2=NB_START,NB_END
              CHAM(N1, N2)=CMPLX(BUFFER((N2-NB_START)*2+1),BUFFER((N2-NB_START)*2+2),q)
           ENDDO


        ENDDO
      END SUBROUTINE READGAMMA

      END MODULE


!*************************SUBROUTINE OUTPOT ****************************
!>   write potential     to a specified unit
!>   HEADER is currently created in pot.F
!>   z is slowest, x is fastest index
!***********************************************************************

      SUBROUTINE OUTPOT(GRIDC, IU, LLONG, CVTOT)
      USE prec
      USE mgrid
      TYPE(grid_3d), INTENT(IN) :: GRIDC              !< grid for potential
      INTEGER, INTENT(IN)       :: IU                 !< unit
      LOGICAL, INTENT(IN)       :: LLONG              !< TRUE: E17.11, FALSE: G11.5
      COMPLEX(q), INTENT(IN)         :: CVTOT(GRIDC%RL%NP) !< real-valued potential
! local variables
      INTEGER NALLOC,NZ, NWRITE, NWRITTEN
      REAL(q),ALLOCATABLE ::  WORK(:)
      INTEGER ISTAT
      CHARACTER (40) FORM

      NODE_ME=0
      IONODE =0

      NODE_ME=GRIDC%COMM%NODE_ME
      IONODE =GRIDC%COMM%IONODE


      NALLOC=GRIDC%NGX*GRIDC%NGY

      ALLOCATE(WORK(NALLOC),STAT=ISTAT)
      IF (ISTAT>0) RETURN ! can not write the potential immediate exit

      IF (GRIDC%NPLWV/= GRIDC%NGX*GRIDC%NGY*GRIDC%NGZ) THEN
        CALL vtutor%bug("internal ERROR: OUTPOT GRIDC%NPLWV is not compatibel with GRIDC%NGX,&
           &GRIDC%NGY,NGZC \n " // str(GRIDC%NPLWV) // " " // str(GRIDC%NGX) // " " // str(GRIDC%NGY) &
           // " " // str(GRIDC%NGZ), "fileio.F", 1989)
      ENDIF

      IF (LLONG) THEN
        FORM='(1(1X,E17.11))'
        NWRITE=5
      ELSE
        FORM='(1G11.5)'
        NWRITE=10
      ENDIF

      IF (NODE_ME==IONODE) WRITE(IU,'(3I5)') GRIDC%NGX,GRIDC%NGY,GRIDC%NGZ

      NWRITTEN=0
      DO NZ=1,GRIDC%NGZ
         CALL MRG_GRID_RL_PLANE(GRIDC, WORK, CVTOT, NZ)
         IF (NODE_ME==IONODE) THEN
         DO N=1,NALLOC
            NWRITTEN=NWRITTEN+1
            IF ( MOD(NWRITTEN,NWRITE)==0 ) THEN
               WRITE(IU,FORM) WORK(N)
            ELSE
               WRITE(IU,FORM,ADVANCE='NO') WORK(N)
            ENDIF
         ENDDO
         ENDIF
      ENDDO
      IF ( MOD(NWRITTEN,NWRITE)/=0 ) WRITE(IU,*)' '

      DEALLOCATE(WORK)

      RETURN
      END SUBROUTINE



!***********************************************************************
!>  write out initial header for PCDAT
!***********************************************************************

      SUBROUTINE PCDAT_HEAD(IU, T_INFO, LATT_CUR, DYN, PACO, SZNAM1)
      USE prec
      USE lattice
      USE poscar
      USE base
      IMPLICIT NONE

      INTEGER IU
      TYPE (type_info) :: T_INFO
      TYPE (latt)::       LATT_CUR
      TYPE (dynamics)  :: DYN
      TYPE (paco_struct)  PACO
      CHARACTER (40) SZNAM1
! local variables
      INTEGER I
      REAL(q) AOMEGA
      AOMEGA=LATT_CUR%OMEGA/T_INFO%NIONS

      WRITE(IU,'(4I4,2E15.7)')1,T_INFO%NIONS,1,0,AOMEGA,DYN%TEMP
      WRITE(IU,*) ' CAR '
      WRITE(IU,*) SZNAM1
      WRITE(IU,'(3I4)') 0,0,0
      WRITE(IU,'(2I4)') 1,DYN%NBLOCK
      WRITE(IU,'(3I4)') PACO%NPACO,PACO%NPACO,PACO%NPACO
      WRITE(IU,'(1I4)') PACO%NPACO
      WRITE(IU,'(1E15.7)') 1E-10
      WRITE(IU,'(1E15.7)') PACO%APACO*1E-10/PACO%NPACO
      WRITE(IU,'(1I4)') DYN%NSW/DYN%NBLOCK/DYN%KBLOCK
      WRITE(IU,'(4E15.7)') DYN%POTIM*1E-15,((LATT_CUR%ANORM(I)*1E-10),I=1,3)

      RETURN
      END SUBROUTINE

!***********************************************************************
!>  write out temperature  to PCDAT file
!***********************************************************************

      SUBROUTINE PCDAT_WRITE(IU, PACO, TMEAN0, TMEAN)
      USE prec
      USE base
      IMPLICIT NONE

      INTEGER IU
      TYPE (paco_struct)  PACO
      REAL(q) :: TMEAN0, TMEAN
! local variables
      INTEGER :: I

      WRITE(IU,'(3E15.7)') TMEAN0,TMEAN

! more then (1._q,0._q) species write all data
      IF (UBOUND(PACO%SIPACO,2)>1) THEN
         DO  I=0,PACO%NPACO-1
            WRITE(IU,'(100F7.3)') PACO%SIPACO(I,:)
         ENDDO
      ELSE
! write only total (entry 0 equal entry 1)
         DO  I=0,PACO%NPACO-1
            WRITE(IU,'(100F7.3)') PACO%SIPACO(I,0)
         ENDDO
      ENDIF
      CALL WFORCE(IU)

      END SUBROUTINE

!***********************************************************************
!
! handling of vasp.lock  file
! this file created (open and close)
! there are three functions
! CREATE_VASP_LOCK creates the vasp.lock file on the root node
! REMOVE_VASP_LOCK removes the vasp.lock file on the root node
! WAIT_VASP_LOCK   waits for the file to be created by an external process
!
!***********************************************************************

!> creates the vasp.lock file on the root node
      SUBROUTINE CREATE_VASP_LOCK(MYCOMM)
        USE mpimy
        USE main_mpi
        TYPE (communic) :: MYCOMM


        IF (MYCOMM%NODE_ME==MYCOMM%IONODE) THEN

           OPEN(UNIT=23,FILE=DIR_APP(1:DIR_LEN)//'vasp.lock',STATUS='UNKNOWN')
           WRITE(23,*)'vasp is running happily'
           CLOSE(UNIT=23)

        ENDIF

      END SUBROUTINE


!>  waits for the file to be created by an external process
      SUBROUTINE WAIT_VASP_LOCK(MYCOMM)
        USE mpimy
        USE main_mpi
        TYPE (communic) :: MYCOMM

! local
        INTEGER :: ISTAT

        DO
           ISTAT=0
! try opening file on root node

           IF (MYCOMM%NODE_ME==MYCOMM%IONODE) THEN

              OPEN(UNIT=23,FILE=DIR_APP(1:DIR_LEN)//'vasp.lock',STATUS='OLD',IOSTAT=ISTAT)

           ENDIF

! "propagate" error status to all nodes
           CALL M_sum_i(MYCOMM, ISTAT,1)
           IF (ISTAT==0) EXIT
! all nodes sleep for 2 seconds
           CALL SLEEP(1)
! all nodes wait for root node
           CALL M_barrier(MYCOMM)
        ENDDO

        CLOSE(UNIT=23)

      END SUBROUTINE WAIT_VASP_LOCK

!******************** SUBROUTINE FILE_EXISTS ***************************
!> @brief Check if a file exists in the current directory
!***********************************************************************
      SUBROUTINE FILE_EXISTS(MYCOMM, FILENAME, RES)
        USE mpimy
        USE main_mpi
        TYPE (communic) :: MYCOMM
        CHARACTER(LEN=*), INTENT(IN) :: FILENAME
        LOGICAL, INTENT(OUT) :: RES
! local
        LOGICAL :: FILE_PRESENT
        INTEGER :: FILE_INT


           IF (MYCOMM%NODE_ME==MYCOMM%IONODE) THEN

           INQUIRE(FILE=TRIM(FILENAME), EXIST=FILE_PRESENT)
! cast to int to make broadcast easier
           FILE_INT = FILE_PRESENT

           ENDIF

         CALL M_bcast_i(COMM%MPI_COMM, FILE_INT , 1 )
! cast back to logical for easier evaluation
         RES = FILE_INT

      END SUBROUTINE FILE_EXISTS

!> removes the vasp.lock file on the root node
      SUBROUTINE REMOVE_VASP_LOCK(MYCOMM)
        USE mpimy
        USE main_mpi
        TYPE (communic) :: MYCOMM

! local
        INTEGER :: ISTAT

        IF (MYCOMM%NODE_ME==MYCOMM%IONODE) THEN

           OPEN(UNIT=23,FILE=DIR_APP(1:DIR_LEN)//'vasp.lock',STATUS='UNKNOWN')
           WRITE(23,*) 'file will be removed'
           CLOSE(UNIT=23,STATUS='DELETE')

        ENDIF

! all nodes sleep for 2 seconds to be sure NFS has propagated this change
        CALL SLEEP(1)
! all nodes wait for root node
        CALL M_barrier(MYCOMM)

      END SUBROUTINE REMOVE_VASP_LOCK
