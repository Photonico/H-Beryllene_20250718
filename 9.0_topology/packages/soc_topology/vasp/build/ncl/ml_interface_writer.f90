# 1 "ml_interface_writer.F"
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


# 2 "ml_interface_writer.F" 2 

MODULE ml_interface_writer
   USE mpimy
   USE prec
   USE tutor, ONLY: vtutor, isError, RandomSeed, argument
   
   CONTAINS

!=======================================================================
!
! This subroutine does some post processing after the energy, forces
! and stress was predicted in a previous step. This routine also does
! the writing to the OUTCAR file.
!
!=======================================================================


   SUBROUTINE ML_TO_VASP_WRITING(DYN, GRIDC, IO, LATT_CUR, LREMOVE_DRIFT, LWRITE_FORCE, &
                                 NSTEP, T_INFO, TOTEN, TIFOR, TSIF)
      USE base
      USE classicfields
      USE constant
      USE ini
      USE lattice
      USE poscar
      USE mgrid
      USE mpimy
      USE msymmetry
      USE poscar
      USE ML_FF_STRUCT, ONLY: ML_SUPER_TYPE, ML_IO_WRITE
      USE ML_FF_CONSTANT, ONLY: MUNIT, TUNIT
! Input and output variables
      TYPE (dynamics)               :: DYN
      TYPE (grid_3d)                :: GRIDC             ! grid for potentials/charge
      TYPE (in_struct)              :: IO
      TYPE (latt), INTENT(IN)       :: LATT_CUR          ! Bravais lattice information
      TYPE (type_info), INTENT(IN)  :: T_INFO            ! type and position information
      LOGICAL, INTENT(IN)           :: LREMOVE_DRIFT     
      LOGICAL, INTENT(IN)           :: LWRITE_FORCE      ! If LWRITE_FORCE=.TRUE., write forces in OUTCAR file.
      INTEGER, INTENT(IN)           :: NSTEP             ! MD step
      REAL (q), INTENT(OUT)         :: TOTEN             ! total energy
      REAL (q), INTENT(OUT)         :: TIFOR(:,:)        !(1:3,1:NIOND) ! forces on all atoms
      REAL (q), INTENT(OUT)         :: TSIF(:,:)         !(1:3,1:3)      ! stress tensor
! Local variables
      REAL(q)   :: EV2KB
      REAL(q)   :: FAC
      INTEGER   :: I
      INTEGER   :: INIONS
      INTEGER   :: IXYZ
      INTEGER   :: JXYZ
      INTEGER   :: NOFFS
      INTEGER   :: NI
      INTEGER   :: NT
      REAL(q)   :: OFIELD_E                        ! energy contribution from order field
      REAL(q)   :: OFIELD_FOR(1:3,T_INFO%NIONS)
      REAL(q)   :: PRESS
      REAL(q)   :: pressure
      REAL(q)   :: TMP(1:6)
      REAL(q)   :: VEL(1:3)
      REAL(q)   :: VTMP(1:3)


! Constants.
      EV2KB = EVTOJ*1E+22_q/LATT_CUR%OMEGA
      FAC  = AMTOKG* &
             1E5_q* &
             1E5_q* &
             1E30_q* &
             1E-8_q
! Add the external order field, if necessary
      CALL OFIELD(GRIDC%COMM, IO, LATT_CUR, DYN, T_INFO, OFIELD_FOR, OFIELD_E)
      TIFOR=TIFOR+OFIELD_FOR
      TOTEN=TOTEN+OFIELD_E
! Remove drift.
      IF (DYN%IBRION/=0) THEN
! remove drift from the forces
         IF (LREMOVE_DRIFT) CALL SYMVEC(T_INFO%NIONS,TIFOR)
      ENDIF
! Only write at every DYN%NBLOCK steps
      IF (MOD(NSTEP,DYN%NBLOCK)==0) THEN
! Output to OUTCAR file.
         IF (IO%IU6>=0) THEN
            PRESS=(TSIF(1,1)+TSIF(2,2)+TSIF(3,3))/3._q- &
                  DYN%PSTRESS/(EVTOJ*1E22_q)*LATT_CUR%OMEGA
            WRITE(IO%IU6,1)       (TSIF(IXYZ,IXYZ),IXYZ=1,3),      TSIF (1,2),     TSIF(2,3),     TSIF (3,1), &
                                  (EV2KB*TSIF(IXYZ,IXYZ),IXYZ=1,3),EV2KB*TSIF(1,2),EV2KB*TSIF(2,3),EV2KB*TSIF(3,1), &
                                   PRESS*EV2KB,DYN%PSTRESS
            TMP = 0._q
            NOFFS = 0
            DO NT=1,T_INFO%NTYP
               DO NI=1+NOFFS,T_INFO%NITYP(NT)+NOFFS
                  VEL(1) = DYN%VEL(1,NI)/DYN%POTIM
                  VEL(2) = DYN%VEL(2,NI)/DYN%POTIM
                  VEL(3) = DYN%VEL(3,NI)/DYN%POTIM
                  CALL  DIRKAR( 1, VEL, LATT_CUR%A )
                  TMP(1) = TMP(1) + VEL(1)*VEL(2) * T_INFO%POMASS(NT)
                  TMP(2) = TMP(2) + VEL(2)*VEL(3) * T_INFO%POMASS(NT)
                  TMP(3) = TMP(3) + VEL(3)*VEL(1) * T_INFO%POMASS(NT)
                  TMP(4) = TMP(4) + VEL(1)*VEL(1) * T_INFO%POMASS(NT)
                  TMP(5) = TMP(5) + VEL(2)*VEL(2) * T_INFO%POMASS(NT)
                  TMP(6) = TMP(6) + VEL(3)*VEL(3) * T_INFO%POMASS(NT)
               ENDDO
               NOFFS = NOFFS + T_INFO%NITYP(NT)
            ENDDO
            TMP = TMP/LATT_CUR%OMEGA*FAC         ! NOW TMP IS IN KB
            WRITE(IO%IU6,'(''  kinetic pressure (ideal gas correction) = '',F9.2,'' kB'')') (TMP(4)+TMP(5)+TMP(6))/3
            pressure=(TSIF(1,1)+TSIF(2,2)+TSIF(3,3))*EV2KB
            pressure=pressure + (TMP(4)+TMP(5)+TMP(6))
            pressure=pressure/3.
            WRITE(IO%IU6,'(''  total pressure  = '',F9.2,'' kB'')') pressure
            WRITE(IO%IU6,'(''  Total+kin. '',F9.3,5F12.3)') TSIF(1,1)*EV2KB + TMP(4),  &
                                                            TSIF(2,2)*EV2KB + TMP(5),  &
                                                            TSIF(3,3)*EV2KB + TMP(6),  &
                                                            TSIF(1,2)*EV2KB + TMP(1),  &
                                                            TSIF(2,3)*EV2KB + TMP(2),  &
                                                            TSIF(3,1)*EV2KB + TMP(3)
            WRITE(IO%IU6,2) LATT_CUR%OMEGA, &
                 ((LATT_CUR%A(IXYZ,JXYZ),IXYZ=1,3),(LATT_CUR%B(IXYZ,JXYZ),IXYZ=1,3),JXYZ=1,3), &
                 (LATT_CUR%ANORM(IXYZ),IXYZ=1,3),(LATT_CUR%BNORM(IXYZ),IXYZ=1,3)
            IF(LWRITE_FORCE) THEN
               WRITE(IO%IU6,3)
               DO INIONS=1, T_INFO%NIONS
                  VTMP(1)=T_INFO%POSION(1,INIONS)
                  VTMP(2)=T_INFO%POSION(2,INIONS)
                  VTMP(3)=T_INFO%POSION(3,INIONS)
                  CALL  DIRKAR(1,VTMP(1),LATT_CUR%A(1,1))
                  WRITE(IO%IU6,4) (VTMP(I),I=1,3),(TIFOR(I,INIONS),I=1,3)
               ENDDO
               VTMP(1:3)=0.0_q
               DO INIONS=1, T_INFO%NIONS
                  DO IXYZ=1, 3
                     VTMP(IXYZ)=VTMP(IXYZ)+TIFOR(IXYZ,INIONS)
                  ENDDO
               ENDDO
               WRITE(IO%IU6,5) VTMP
            ENDIF
            WRITE(IO%IU6,6)
            WRITE(IO%IU6,7) '  ML FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)',TOTEN,TOTEN,TOTEN
         ENDIF
      ENDIF
! Convert to enthalpy
      IF (DYN%PSTRESS/=0) THEN
         TOTEN=TOTEN+DYN%PSTRESS/(EVTOJ*1E22_q)*LATT_CUR%OMEGA
! Again printing only at every NBLOCK step
         IF (MOD(NSTEP,DYN%NBLOCK)==0) THEN
            IF (IO%IU6>=0) WRITE(IO%IU6,8) TOTEN,DYN%PSTRESS/(EVTOJ*1E22_q)*LATT_CUR%OMEGA
         ENDIF
      ENDIF
! Format
1        FORMAT(/'  ML FORCE on cell =-STRESS in cart. coord. units (eV/cell)'/ &
                 '  Direction',4X,'XX', 10X,'YY', 10X,'ZZ', 10X,'XY', 10X,'YZ', 10X,'ZX'/ &
                 '  --------------------------------------------------------------------------------------'/ &
                 '  Total: ',6F12.5/ &
                 '  in kB  ',6F12.5/, &
                 '  external pressure = ',F11.2,' kB', &
                 '  Pullay stress = ',F11.2,' kB'/)
2        FORMAT( '  volume of cell :  ',F10.2/ &
                 '      direct lattice vectors',17X,'reciprocal lattice vectors'/ &
                 3(2(3X,3F13.9)/)/ &
                 '  length of vectors'/ &
                  (2(3X,3F13.9)/) /)
3        FORMAT( '  POSITION    ',35X,'TOTAL-FORCE (eV/Angst) (ML)'/ &
                 ' -----------------------------------------------------------------------------------')
4        FORMAT((3F13.5,3X,3F14.6))
5        FORMAT( ' ----------------------------------------------', &
                 '-------------------------------------',/ &
                 '    total drift:      ',20X,3F14.6)
6        FORMAT (5X, //, &
                '----------------------------------------------------', &
                '----------------------------------------------------'//)
7        FORMAT(/ &
                A/ &
                '  ---------------------------------------------------'/ &
                '  free  energy ML TOTEN  = ',F18.8,' eV'// &
                '  ML energy  without entropy=',F18.8, &
                '  ML energy(sigma->0) =',F18.8/)
8        FORMAT ('  enthalpy is ML TOTEN    = ',F18.8,' eV   P V=',F18.8/)
   END SUBROUTINE ML_TO_VASP_WRITING

END MODULE ml_interface_writer
