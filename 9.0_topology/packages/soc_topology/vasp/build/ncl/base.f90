# 1 "base.F"
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


# 2 "base.F" 2 
!************************************************************************
! RCS:  $Id: base.F,v 1.2 2001/02/20 14:44:56 kresse Exp $
!
!> this module contains some control data structures for VASP
!
!***********************************************************************
      MODULE PREC
      INTEGER, PARAMETER :: q =SELECTED_REAL_KIND(10)

      INTEGER, PARAMETER :: qd=SELECTED_REAL_KIND(30)
# 14

      INTEGER, PARAMETER :: qs=SELECTED_REAL_KIND(5)
      INTEGER, PARAMETER :: qi4=SELECTED_INT_KIND(8)
      INTEGER, PARAMETER :: qi8=SELECTED_INT_KIND(15)
!
!> this parameter controls the step width in for numerical differentiation
!> in some VASP routines
!> 1E-5 is very reliable yielding at least 7 digits in all contributions
!> to the forces
!> 1E-4, however, is better suited for second derivatives
!> for reasons of consistentcy with previous versions 1E-5
      REAL(q), PARAMETER :: fd_displacement=1E-5
      END MODULE

      MODULE BASE
      USE prec
      IMPLICIT NONE
!
!> header type information
!
      TYPE PLUGINS_SETTINGS
!> Configuration of the interface to plugins
        CHARACTER(LEN=:), ALLOCATABLE :: MODE !< are the plugins executed in serial or parallel
        LOGICAL :: FORCE_AND_STRESS = .FALSE. !< call the force and stress plugins
        LOGICAL :: LOCAL_POTENTIAL = .FALSE. !< call the local potential plugins
        LOGICAL :: MACHINE_LEARNING = .FALSE.  !< augment the VASP run by input from an external ML
        LOGICAL :: OCCUPANCIES = .FALSE. !< allow modification of occupation related parameters
        LOGICAL :: STRUCTURE = .FALSE. !< call the structure update plugins
        INTEGER :: OUTBLOCK = 1 !< Step size of the
        INTEGER :: OUTPUT_MODE = 0 !< Controls whether PCDAT etc. files are calculated and written or not
      END TYPE PLUGINS_SETTINGS

      TYPE info_struct
! only INFO
! mode information
        LOGICAL LREAL                !< real space projection/ reciprocal space proj.
        LOGICAL LOVERL               !< vanderbilt type PP read in ?
        LOGICAL LCORE                !< any partial core read in ?
        LOGICAL LCHCON               !< charge density constant during run
        LOGICAL LCHCOS               !< allways same as above in cur. impl.
        LOGICAL LONESW               !< use all bands simultaneous
        LOGICAL LONESW_AUTO          !< switch automatically between LONESW and DIIS
        LOGICAL LPRECONDH            !< precondition the subspace rotation matrix
        LOGICAL LDAVID               !< use block davidson
        LOGICAL LDAVID_FULL          !< use non-blocked davidson
        LOGICAL LEXACT_DIAG          !< use exact diagonalization
        LOGICAL LRMM                 !< use RMM-DIIS algorithm
        LOGICAL LORTHO               !< orthogonalize
        LOGICAL LABORT,LSOFT         !< soft / hard stop
        LOGICAL LSTOP                !< stop with this iteration
        LOGICAL LPOTOK               !< local potential ok ?
        LOGICAL LMIX                 !< mixing 1._q
        LOGICAL LASPH                !< aspherical radial PAW
! exchange correlation spin
        LOGICAL LXLDA                !< calculate LDA exchange only in POTLOK
        INTEGER ISPIN                !< spin 1=no 2 =yes
        REAL(q) RSPIN                !< 2 for spinpolar. 1 for non spinpol.
! electronic relaxation
        INTEGER ISTART               !< how to start up
        INTEGER ICHARG               !< initial charge density
        INTEGER INIWAV               !< how to initialize wavefunctions
        INTEGER INICHG               !< how to initialize charge
        INTEGER NELM                 !< maximal number of el-steps in SCC
        INTEGER NELMALL              !< maximal number of el-steps in direct optimization
        INTEGER NELMIN               !< minimal number of el-steps
        INTEGER NELMDL               !< number of delay el-steps
        INTEGER IALGO                !< algorithm for el-relax
        INTEGER IALGO_COMPAT         !< compatibility flag for IALGO
!< MOD(IALGO_COMPAT,2) = 1 use old DIIS algorithm
        INTEGER NDAV                 !< number of steps in RMM-DIIS
        INTEGER NKRYLOV              !< maximum size of reduced basis in the non-blocked Davidson
        INTEGER DAVITER              !< maximum number of iterations in the non-blocked Davidson
        INTEGER NSTRICT              !< the number of orbitals that will be converged tightly
        REAL(q) FBREAK               !< break-off criterion w.r.t. residual size
        REAL(q) TIME                 !< timestep for el (IALGO>50)
        REAL(q) WEIMIN,EBREAK,DEPER  !< control tags for elm
        REAL(q) EDIFF                !< accuracy for electronic relaxation
        LOGICAL LDIAG                !< level reordering allowed or not
        LOGICAL LSUBROT              !< sub space rotation to optimize rotation matrix
        LOGICAL LPDIAG               !< sub space rotation before iterat. diag.
        LOGICAL LCDIAG               !< recalculate eigenvalues after iterat. diag.
! cutoff information
        REAL(q) ENMAX                !< cutoff for calculations
        REAL(q) ENINI                !< cutoff during delay
        REAL(q) ENAUG                !< cutoff for augmentation charges
! some important things
        REAL(q) EALLAT               !< total energy of all atoms
        REAL(q) NELECT               !< number of electrons
        REAL(q) NUP_DOWN             !< spin multiplicity
        INTEGER NBANDTOT             !< total number of bands
        INTEGER MCPU                 !< max number of proc (dimensioned)
        INTEGER NCPU                 !< actual number of proc
        LOGICAL LCORR                !< correction to forces
        INTEGER IQUASI               !< eigenvalue/occupation number corrections
        INTEGER TURBO                !< turbo mode
        INTEGER IFOLD                !< folded eigenproblem: 1=LFOLD,2=LFOLDHpsi
        INTEGER IRESTART             !< whether to restart: 2=restart with 2 optimized vectors
        INTEGER IHARMONIC            !< harmonic Ritz values
        INTEGER NREBOOT              !< number of reboots
        INTEGER NMIN                 !< reboot dimension
        REAL(q) EREF                 !< reference energy to select bands
        LOGICAL NLSPLINE             !< use spline interpolation to construct projection operator
! GW algorithms
        INTEGER ICHIREAL             !< which GW algorithm , > 0 selects space time algos
        LOGICAL LGW                  !< perform GW calculations
        LOGICAL LCHI                 !< calculate response functions
        LOGICAL LscQPGW              !< iterate wavefunctions (selfconsistent GW)
        LOGICAL LGW0                 !< leave initial wavefunctions unmodified for calculations of W
        LOGICAL LG0W0                !< leave initial wavefunctions unmodified for calculations of W
        LOGICAL LCRPA                !< swithes on CRPA in GW routines
        LOGICAL LBSE                 !< solve BSE or Cassida equations (TD-DFT)
        LOGICAL LACFDT               !< ACFDT calculation (adiabatic connection, fluctuation dissipation theorem)
        LOGICAL LOEP                 !< OEP calculation: calculate the RPA-OEP potential
        LOGICAL LEXX                 !< OEP calculation: calculate the EXX-OEP potential
        LOGICAL LHFCORRECT           !< HFCORRECT calculation (adiabatic connection, fluctuation dissipation theorem)
        INTEGER EXXOEP               !< Energy cutoff settings for OEP/LHF methods
        LOGICAL LGWNO                !< GW natural orbitals (optimized for description of Coloumb hole)
        LOGICAL LCORBSE              !< solve BSE or Cassida equations and calculate correlation energy
        LOGICAL LQPBSE               !< solve BSE or Cassida equations and calculate self-energy
        LOGICAL L2E4W                !< Calculate 2-electron 4-wannier-orbital integrals
        LOGICAL L2E4W_ALL            !< all terms not just (1._q,0._q)-center terms
        TYPE(PLUGINS_SETTINGS) PLUGIN  !< which quantities may be altered by plugins
! characters allways last
        CHARACTER*40 SZNAM1          !< header of INCAR
        CHARACTER*12 SZPREC          !< precision information
      END TYPE

!> This type controls which potential are written
      TYPE write_pot
        LOGICAL :: LVTOT = .FALSE.    !< legacy flag to write total local potential
        LOGICAL :: LVHAR = .FALSE.    !< legacy flag to write Hartree potential
        LOGICAL :: TOTAL = .FALSE.    !< write total local potential
        LOGICAL :: HARTREE = .FALSE.  !< write hartree potential
        LOGICAL :: IONIC = .FALSE.    !< write ionic potential
        LOGICAL :: XC = .FALSE.       !< write xc potential
        LOGICAL :: PAW = .FALSE.      !< write paw strengths
      END TYPE write_pot

      TYPE in_struct
! only IO
        LOGICAL LOPEN                !< files open at startup
        INTEGER IU0                  !< unit for error
        INTEGER IU6                  !< unit for stdout
        INTEGER IU5                  !< unit for stdin
        INTEGER NWRITE               !< how much information is written out
        INTEGER IDIOT                !< how much information is written out
        INTEGER ICMPLX               !< size of a complex item upon IO
        INTEGER MRECL                !< maximal size of record length
        LOGICAL LREALD               !< no LREAL read in
        LOGICAL LMUSIC               !< jF (just a joke)
        LOGICAL LFOUND               !< WAVECAR exists ?
        LOGICAL LWAVE                !< write WAVECAR
        LOGICAL LCHARG               !< write CHGCAR
        LOGICAL LPDENS               !< write partial density (charge density for (1._q,0._q) band)
        INTEGER LORBIT               !< write orbit/dos
        LOGICAL LELF                 !< write elf
        LOGICAL LOPTICS              !< calculate/write optical matrix elements
        LOGICAL LPETIM               !< timing information
        INTEGER IUVTOT               !< unit for local potential
        LOGICAL INTERACTIVE          !< vasp runs interactive
        INTEGER IRECLW               !< record lenght for WAVECAR
        LOGICAL LDOWNSAMPLE          !< read WAVECAR of denser k-grid
        LOGICAL LWAVEDERF            !< write file WAVEDERF
        LOGICAL :: DRY_RUN = .FALSE. !< execute dry run of VASP
        LOGICAL :: VELOCITY = .FALSE.!< write velocities to HDF5 file
        LOGICAL LWAP                 !< activate writting of the potential for electron-phonon calculation
# 182

        LOGICAL :: LH5               !< redirect restart information to vaspwave
        LOGICAL :: LCHARGH5          !< write the density to vaspwave.h5
        LOGICAL :: LPARCHGH5         !< write partial charge densities to vaspout.h5
        LOGICAL :: LSYNCH5 = .FALSE. !< synchronize the vaspout.h5 during run to allow SWMR access
        TYPE(write_pot) WRT_POTENTIAL !< which potentials are written
      END TYPE

      TYPE mixing
! only MIX
        INTEGER IUBROY               !< unit for broyden mixer
        REAL(q) AMIX                 !< mixing parameter A
        REAL(q) BMIX                 !< mixing parameter B
        REAL(q) AMIX_MAG             !< mixing parameter A for magnetization
        REAL(q) BMIX_MAG             !< mixing parameter B for magnetization
        REAL(q) AMIN                 !< minimal mixing parameter A
        REAL(q) WC                   !< weight factor for Johnsons method
        INTEGER IMIX                 !< type of mixing
        INTEGER INIMIX               !< initial mixing matrix
        INTEGER MIXPRE               !< form of metric for mixing
        LOGICAL LRESET               !< reset mixer on next call (set when ions move)
        LOGICAL HARD_RESET           !< force hard reset of mixer (force full reset regardless of MAXMIX)
        INTEGER MAXMIX               !< maximum number of mixing steps (if positive LRESET does not apply)
        INTEGER NEIG                 !< number of eigenvalues
        INTEGER MREMOVE              !< how many vectors are removed once the iteration depth is reached
        REAL(q) EIGENVAL(512)        !< eigenvalues of dielectric matrix
        REAL(q) AMEAN                !< mean eigenvalue
        LOGICAL MIXFIRST             !< mix before diagonalization (or after)
      END TYPE

      TYPE symmetry
!only SYMM
        INTEGER, POINTER :: ROTMAP(:,:,:) !
        REAL(q), POINTER :: TAU(:,:)      ! jF
        REAL(q), POINTER :: TAUROT(:,:)   ! jF
        REAL(q), POINTER :: WRKROT(:)     ! jF
        REAL(q), POINTER :: PTRANS(:,:)   ! jF
        REAL(q), POINTER :: MAGROT(:,:)   ! jF
        INTEGER, POINTER :: INDROT(:)     ! jF
        INTEGER ISYMOP(3,3,48)            !< Space group symmetry operation
        INTEGER IGRPOP(3,3,48)            !< Rotation part of the sapce group symmetry operation
        REAL(q) GTRANS(3,48)              !< Translation part of the space group symmetry operation
        INTEGER INVMAP(48)                !< Map to the inverse elements of each group element
        REAL(q) AP(3,3)                   !< Lattice parameters of the primitive cell
        INTEGER NPCELL                    !< Number of primitive cells
        INTEGER ISYM                      !< symmetry on/of
        INTEGER NROT                      !< number of rotations
        INTEGER NPTRANS                   !< number of primitive translations
      END TYPE

      TYPE prediction
!only PRED
        INTEGER IWAVPR               !< prediction of wavefunctions
        INTEGER INIPRE               !< initialized yes/no
        INTEGER IPRE                 !< what was 1._q in wavefunction predic.
        INTEGER IUDIR                !< unit for prediction of wavefunction
        INTEGER ICMPLX               !< size of complex word
        REAL(q)  ALPHA,BETA
      END TYPE

      TYPE dipol
!only DIP
        INTEGER IDIPCO               !< direction (0 no dipol corrections)
        LOGICAL LCOR_DIP             !< correct potential
        REAL(q) POSCEN(3)            !< position of center
        REAL(q) DIPOLC(3)            !< calculated dipol
        REAL(q) QUAD                 !< trace of quadrupol
        INTEGER INDMIN(3)            !< position of minimum
        REAL(q) EDIPOL,EMONO,E_ION_EXTERN
        REAL(q), POINTER :: FORCE(:,:)
        REAL(q) VACUUM(2)            !< vacuum level
      END TYPE

      TYPE smear_struct
!only SMEAR_LOOP
        INTEGER              :: ISMCNT          !
        REAL(q), allocatable :: SMEARS(:)       !< (ismear, sigma, ...)
      END TYPE


      TYPE paco_struct
!only PACO
        INTEGER NPACO                !< number of grid points for pair corr.
        REAL(q) APACO                !< cutoff
        REAL(q),ALLOCATABLE :: SIPACO(:,:) ! accumulated partial pair correlation function, second index is flattened
        INTEGER SMEANP               !< number of configurations stored
        REAL(q) OMEGA                !< average volume times number of accumulated timesteps
      END TYPE

      TYPE energy
        REAL(q)    :: TOTENASPH  = 0.0_q  !< total energy for aspherical GGA
        REAL(q)    :: EBANDSTR   = 0.0_q  !< bandstructure energy
        REAL(q)    :: DENC       = 0.0_q  !< -1/2 hartree (d.c.)
        REAL(q)    :: XCENC      = 0.0_q  !< -V(xc)+E(xc) (d.c.)
        REAL(q)    :: EXCG       = 0.0_q  !< E(xc) (LDA+GGA)
        REAL(q)    :: EXLDA      = 0.0_q  !< LDA excchange energy
        REAL(q)    :: ECLDA      = 0.0_q  !< LDA correlation energy
        REAL(q)    :: EXGGA      = 0.0_q  !< GGA exchange energy
        REAL(q)    :: ECGGA      = 0.0_q  !< GGA correlation energy
        REAL(q)    :: EXHF       = 0.0_q  !< Hartree-Fock exchange energy
        REAL(q)    :: EXHF_ACFDT = 0.0_q  !< difference between HF energy, and exchange energy in ACFDT
        REAL(q)    :: EDOTP      = 0.0_q  !< Electric field \dot Polarization
        REAL(q)    :: TEWEN      = 0.0_q  !< Ewald energy
        REAL(q)    :: PSCENC     = 0.0_q  !< alpha Z (V(q->0) Z)
        REAL(q)    :: EENTROPY   = 0.0_q  !< Entropy term
        REAL(q)    :: PAWPS      = 0.0_q  !< paw double counting corrections
        REAL(q)    :: PAWAE      = 0.0_q  !< paw double counting corrections
        REAL(q)    :: PAWPSG     = 0.0_q  !< paw xc energies (LDA+GGA)
        REAL(q)    :: PAWAEG     = 0.0_q  !< paw xc energies (LDA+GGA)
        REAL(q)    :: PAWCORE    = 0.0_q  !< exchange correlation energy of core (LDA+GGA)
        REAL(q)    :: PAWPSAS    = 0.0_q  !< paw xc energies (aspherical)
        REAL(q)    :: PAWAEAS    = 0.0_q  !< paw xc energies (aspherical)
        COMPLEX(q) :: CVZERO     = (0.0_q, 0.0_q)  !< average local potential
        REAL(q)    :: ECGWGM     = 0.0_q  !< GW Galitskii-Migdal correlation energy
        REAL(q)    :: EKLGG0     = 0.0_q  !< Klein-contribution: Tr( G_0^-1 . G - 1 )
        REAL(q)    :: EKLLOG     = 0.0_q  !< Klein-contribution: Tr[Ln( G . G_0^-1 )]
        REAL(q)    :: ETRGHF     = 0.0_q  !< trace of HF hamiltonian with denstiy matrix
        REAL(q)    :: ELOG1G     = 0.0_q  !< Tr Log ( 1 - G S )
        REAL(q)    :: ELOGG0     = 0.0_q  !< Tr Log ( G_0 )
        REAL(q)    :: ERPA_CUT   = 0.0_q  !< RPA corr. energy with highest cutoff caculated
        REAL(q)    :: ERPA_INF   = 0.0_q  !< RPA corr. energy extrapolated to infinte basis set cutoff
        REAL(q)    :: ESCPC      = 0.0_q  !< self consistent potential correction double counting
        REAL(q)    :: EPLUGINS   = 0.0_q  !< contributions from plugins (python)
        REAL(q)    :: POTHARZERO = 0.0_q  !< g=0 component for the hartree potential
        REAL(q)    :: POTIONZERO = 0.0_q  !< g=0 component for the ion potential
      END TYPE

      CHARACTER (LEN=10) :: INCAR='INCAR'

!> Create a pointer to complex array reinterpreting it as a real array
      INTERFACE CAST_COMPLEX_TO_REAL
         MODULE PROCEDURE C_D_CONVERT1, C_D_CONVERT2
      END INTERFACE CAST_COMPLEX_TO_REAL

!> Create a pointer to complex array reinterpreting it as an COMPLEX(q) array
      INTERFACE CAST_COMPLEX_TO_RGRID
         MODULE PROCEDURE C_RG_CONVERT1, C_RG_CONVERT2
      END INTERFACE

!> Create a pointer to COMPLEX(q) array reinterpreting it as a complex
      INTERFACE CAST_RGRID_TO_COMPLEX
         MODULE PROCEDURE RG_C_CONVERT1
      END INTERFACE

      PRIVATE RG_C_CONVERT1, C_RG_CONVERT1, C_RG_CONVERT2, C_D_CONVERT1, C_D_CONVERT2

      CONTAINS
!
!> small subroutine which tries to give good dimensions for 1 dimension
!
      SUBROUTINE MAKE_STRIDE (N)
      INTEGER N,NEW
      INTEGER, PARAMETER :: NGOOD=16

      NEW=(N+NGOOD)/NGOOD
      NEW=NEW*NGOOD+1
      N=NEW

      END SUBROUTINE

!****************************************************************************
! DESCRIPTION:
!>returns a double-precision real from a possible emulated quad precision
!>real
!> @param[in]  A  (argument)
!****************************************************************************
   FUNCTION TOREAL( A )
      USE prec
# 352

      REAL(q) :: TOREAL
      REAL(qd) :: A
      TOREAL = A
   END FUNCTION TOREAL

      SUBROUTINE C_D_CONVERT1(C, R)
         USE ISO_C_BINDING, ONLY: C_F_POINTER, C_LOC
         COMPLEX(Q), INTENT(IN), TARGET :: C(:)
         REAL(Q), INTENT(INOUT), POINTER :: R(:)

         CALL C_F_POINTER(C_LOC(C), R, SHAPE=[2*SIZE(C)])

      END SUBROUTINE C_D_CONVERT1

      SUBROUTINE C_D_CONVERT2(C, R)
         USE ISO_C_BINDING, ONLY: C_F_POINTER, C_LOC
         COMPLEX(Q), INTENT(IN), TARGET :: C(:,:)
         REAL(Q), INTENT(INOUT), POINTER :: R(:,:)

         CALL C_F_POINTER(C_LOC(C), R, SHAPE=[2*SIZE(C,1), SIZE(C,2)])

      END SUBROUTINE C_D_CONVERT2


      SUBROUTINE C_RG_CONVERT1(C, RG)
         COMPLEX(Q), INTENT(IN), TARGET :: C(:)
         COMPLEX(q), INTENT(INOUT), POINTER :: RG(:)
# 382

         RG => C

      END SUBROUTINE C_RG_CONVERT1

      SUBROUTINE C_RG_CONVERT2(C, RG)
         COMPLEX(Q), INTENT(IN), TARGET :: C(:,:)
         COMPLEX(q), INTENT(INOUT), POINTER :: RG(:,:)
# 392

         RG => C

      END SUBROUTINE C_RG_CONVERT2


      SUBROUTINE RG_C_CONVERT1(R, C)
         USE ISO_C_BINDING, ONLY: C_F_POINTER, C_LOC
         COMPLEX(q), INTENT(IN), TARGET :: R(:)
         COMPLEX(Q), INTENT(INOUT), POINTER :: C(:)
# 404

         C => R

      END SUBROUTINE RG_C_CONVERT1

      END MODULE
