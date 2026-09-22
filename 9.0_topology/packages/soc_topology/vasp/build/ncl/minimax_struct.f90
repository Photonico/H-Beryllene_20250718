# 1 "minimax_struct.F"
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


# 2 "minimax_struct.F" 2 
MODULE minimax_struct
   USE prec
   USE base, ONLY: TOREAL
# 7

   USE mpimy
   IMPLICIT NONE
!
!> general descriptor for a quadrature type
!
   TYPE quadrature_handle
!>  number that identifies basis functions used for fit
      INTEGER          :: ITYPE=0
!> grid id
      INTEGER          :: GRID_ID=0
!> number of quadrature points
      INTEGER          :: N=0
!> scaled interval \f$[1,R=B/A]\f$ (or \f$[0,B]\f$ ) for which quadrature error is minimized
      REAL(qd)          :: A, B
!> scaling factor
      REAL(qd)          :: SCALING
!> L_infty norm or max error of quadature
      REAL(qd)          :: INFINITYNORM
!> quadrature coefficients, typically of size 2N+1
      REAL(qd), POINTER :: C(:) => NULL()
!> alternant of error function, typically of size 2N+1
      REAL(qd), POINTER :: X0(:) => NULL()
     
!> function used to fit object function
      PROCEDURE( basis_function ), POINTER, NOPASS :: PHI2 => NULL()
!> first derivative of basis function used for fit w.r.t. first argument
      PROCEDURE( basis_function ), POINTER, NOPASS :: D_PHI2_DZ => NULL()
!> first derivative of basis function used for fit w.r.t. second argument
      PROCEDURE( basis_function ), POINTER, NOPASS :: D_PHI2_DL => NULL()
!> second derivative of basis function used for fit w.r.t. first argument
      PROCEDURE( basis_function ), POINTER, NOPASS :: D2_PHI2_DZ2 => NULL()

!> function \f$ f(z) \f$ that is approximated by basis functions
      PROCEDURE( object_function ), POINTER, NOPASS :: F => NULL()
!> first derivative of function \f$ f(z) \f$ that is approximated by basis functions
      PROCEDURE( object_function ), POINTER, NOPASS :: D_F_DZ => NULL()
!> second derivative of function \f$ f(z) \f$ that is approximated by basis functions
      PROCEDURE( object_function ), POINTER, NOPASS :: D2_F_DZ2 => NULL()

!> subroutines that set initial quadrature values for R=1E3
      PROCEDURE( init_quad_subroutine ), POINTER, NOPASS :: R1E3 => NULL()
!> subroutines that set initial quadrature values for R=1E4
      PROCEDURE( init_quad_subroutine ), POINTER, NOPASS :: R1E4 => NULL()
!> subroutines that set initial quadrature values for R=1E5
      PROCEDURE( init_quad_subroutine ), POINTER, NOPASS :: R1E5 => NULL()
!> subroutines that set initial quadrature values for R=1E6
      PROCEDURE( init_quad_subroutine ), POINTER, NOPASS :: R1E6 => NULL()
!> tabulated coefficients from minimax_ini
      REAL(qd), POINTER :: CTAB(:) => NULL()

!> basis function \f$\phi(z,\lambda) \f$ used to fit function and derivatives
      PROCEDURE( basis_function ), POINTER, NOPASS :: PHI => NULL()
!> dual basis function \f$\psi(\lambda,z) \f$ used to fit function and derivatives
      PROCEDURE( basis_function ), POINTER, NOPASS :: PSI => NULL()
!> auxilary function used for transformation matrix for small frequency points
      PROCEDURE( object_function ), POINTER, NOPASS :: TRANS_TAYLOR => NULL()

!> conjugate basis function
      PROCEDURE( basis_function ), POINTER, NOPASS :: PHI_CONJG => NULL()
!> conjugate basis function
      PROCEDURE( basis_function ), POINTER, NOPASS :: PSI_CONJG => NULL()
!> auxilary function used for transformation matrix for small frequency points
      PROCEDURE( object_function ), POINTER, NOPASS :: TRANS_TAYLOR_CONJG => NULL()

!> vectorial function that fits the minimax coefficients and alternatnt for requested interval length
      PROCEDURE( fit_coeff_function ), POINTER, NOPASS :: FITCOF => NULL()
      REAL(qd)          :: R1 !< interval, for which fit is taylored
      REAL(qd)          :: R2 !

!> communicator used in internal minimax routines for 1 parallelization
      TYPE( communic ) :: COMM

   END TYPE quadrature_handle
!
!> unified imaginary grid descriptor for time and frequency domain
!> contains communicators and other local as well as global indices
!
   TYPE loop_des
      INTEGER :: NPOINTS                           !< total number of time/freq points tau
      INTEGER :: NPOINTSC                          !< loop counter, indexing the current grid  point handled
      INTEGER :: NPOINTS_IN_GROUP                  !< number of points in group
      INTEGER :: NPOINTS_IN_ROOT_GROUP             !< number of points in first group (always the largest)
      REAL(q) :: POINT_CURRENT                     !< imaginary time treated currently (only (1._q,0._q) at each step)
      REAL(q) :: POINTMIN_CURRENT                  !< current smallest tau
      INTEGER :: NPOINT_GLOBAL                     !< global grid point
      LOGICAL :: LDO_POINT_LOCAL                   !< skip GEMM calls?
 
      TYPE (communic) :: COMM_BETWEEN_GROUPS       !< communicator between tau
!< processor groups
      TYPE (communic) :: COMM_IN_GROUP             !< communicator witin tau
!< processor group
      TYPE (communic) :: COMM                      !< global communicator
 
      REAL(q), POINTER :: POINTS_LOCAL(:)=>NULL()  !< local points treated by (1._q,0._q) processor group

      INTEGER, POINTER :: DISTRIBUTION(:,:)=>NULL()!< information how points are distributed

      INTEGER, POINTER :: GROUP_OF_NODE(:)=>NULL() !< processor group id
!< based on global id
!> inverse temperature (in eV), as default take the (1._q,0._q) for 300 K
      REAL(q) :: BETA = 38.68149168875351_q
   END TYPE loop_des
   
!> imaginary grid type used in GW subroutines (high level)
   TYPE imag_grid_handle
!> number of grid points (NOMEGA)
      INTEGER          :: NOMEGA=0
      
!> decides if finite temperature grids are calculated
      LOGICAL          :: LFINITE_TEMPERATURE = .FALSE.
!> inverse temperature (in eV), as default take the (1._q,0._q) for 300 K
      REAL(q)          :: BETA = 38.68149168875351_q
!> minimization interval (left boundary)
      REAL(q)          :: X1=0
!> minimization interval (right boundary)
      REAL(q)          :: X2=0
!> minimization interval (scaled right boundary )
      REAL(q)          :: R=0

!> this is used for CRPA, where set to -1 for subtracting correlated part of polarizability
      REAL(q) :: FACTOR = 1._q

!> integration error for Galitskii-Migdal correlation energy \f$ \int_0^\beta d\tau G(\tau) \Sigma_c(-\tau) \f$
      REAL(q)          :: GM_INT_ERROR=0
      
!> imaginary time grid points
      REAL(q), POINTER :: TAU(:) => NULL( )
!> imaginary time integration weights points
      REAL(q), POINTER :: TAU_WEIGHT(:)=> NULL( )
!> error of imaginary time integration quadrature
      REAL(q)          :: TAU_ERROR=0
!> transformation weights from FER_RE to TAU=0+ point
      REAL(q), POINTER :: TO_TAU0(:,:) => NULL( )
!> transformation error from FER_RE to TAU=0+ point
      REAL(q), POINTER :: TO_TAU0_ERROR(:) => NULL( )
!> transformation weights from FER_RE to TAU=0+ point
      REAL(q), POINTER :: TO_TAU0_CONJG(:,:) => NULL( )
!> transformation error from FER_IM to TAU=0+ point
      REAL(q), POINTER :: TO_TAU0_CONJG_ERROR(:) => NULL( )
!> transformation weights to TAU=0- are same as above with minus sign

!> imaginary frequency points for real part of bosonic functions
      REAL(q), POINTER :: BOS_RE(:) => NULL( )
!> weights for frequency integration of real part of bosonic functions
      REAL(q), POINTER :: BOS_RE_WEIGHT(:)=> NULL( )
!> corresponding error
      REAL(q)          :: BOS_RE_ERROR=0
!> transformation matrix from TAU to BOS_RE grid
      REAL(q), POINTER :: TO_BOS_RE(:,:) => NULL( )
!> error of transformation from TAU to BOS_RE grid points
      REAL(q), POINTER :: TO_BOS_RE_ERROR(:) => NULL( )
!> transformation matrix from TAU to BOS_RE grid ( conjugate basis )
      REAL(q), POINTER :: TO_BOS_RE_CONJG(:,:) => NULL( )
!> error of transformation from TAU to BOS_RE grid points
      REAL(q), POINTER :: TO_BOS_RE_CONJG_ERROR(:) => NULL( )

!> imaginary frequency points for imaginary part of bosonic functions
      REAL(q), POINTER :: BOS_IM(:) => NULL( )
!> weights for frequency integration of imaginary part of bosonic functions
      REAL(q), POINTER :: BOS_IM_WEIGHT(:)=> NULL( )
!> corresponding error
      REAL(q)          :: BOS_IM_ERROR=0
!> transformation matrix from TAU to BOS_RE grid
      REAL(q), POINTER :: TO_BOS_IM(:,:) => NULL( )
!> error of transformation from TAU to BOS_RE grid points
      REAL(q), POINTER :: TO_BOS_IM_ERROR(:) => NULL( )
!> transformation matrix from TAU to BOS_RE grid (conjugate basis )
      REAL(q), POINTER :: TO_BOS_IM_CONJG(:,:) => NULL( )
!> error of transformation from TAU to BOS_RE grid points
      REAL(q), POINTER :: TO_BOS_IM_CONJG_ERROR(:) => NULL( )

!> imaginary frequency points for real part of fermionic functions
      REAL(q), POINTER :: FER_RE(:) => NULL( )
!> weights for frequency integration of real part of fermionic functions
      REAL(q), POINTER :: FER_RE_WEIGHT(:)=> NULL( )
!> corresponding error
      REAL(q)          :: FER_RE_ERROR=0
!> transformation matrix from TAU to FER_RE grid
      REAL(q), POINTER :: TO_FER_RE(:,:) => NULL( )
!> error of transformation from TAU to FER_RE grid points
      REAL(q), POINTER :: TO_FER_RE_ERROR(:) => NULL( )
!> transformation matrix from TAU to FER_RE grid (conjugate basis)
      REAL(q), POINTER :: TO_FER_RE_CONJG(:,:) => NULL( )
!> error of transformation from TAU to FER_RE grid points
      REAL(q), POINTER :: TO_FER_RE_CONJG_ERROR(:) => NULL( )

!> imaginary frequency points for imaginary part of fermionic functions
      REAL(q), POINTER :: FER_IM(:) => NULL( )
!> weights for frequency integration of imaginary part of fermionic functions
      REAL(q), POINTER :: FER_IM_WEIGHT(:)=> NULL( )
!> corresponding error
      REAL(q)          :: FER_IM_ERROR=0
!> transformation matrix from TAU to FER_RE grid
      REAL(q), POINTER :: TO_FER_IM(:,:) => NULL( )
!> error of transformation from TAU to FER_RE grid points
      REAL(q), POINTER :: TO_FER_IM_ERROR(:) => NULL( )
!> transformation matrix from TAU to FER_RE grid (conjugate basis)
      REAL(q), POINTER :: TO_FER_IM_CONJG(:,:) => NULL( )
!> error of transformation from TAU to FER_RE grid points
      REAL(q), POINTER :: TO_FER_IM_CONJG_ERROR(:) => NULL( )
 
!> quadrature handles, 5 distinct grids are possible
      TYPE( quadrature_handle ) :: TIME
!> grid for real part of bosonic correlation functions
      TYPE( quadrature_handle ) :: FREQ_BOS_RE
!> grid for imaginary part of bosonic correlation functions
      TYPE( quadrature_handle ) :: FREQ_BOS_IM
!> grid for real part of fermionic correlation functions
      TYPE( quadrature_handle ) :: FREQ_FER_RE
!> grid for imaginary part of fermionic correlation functions
      TYPE( quadrature_handle ) :: FREQ_FER_IM

!> loop descriptors for time grid
      TYPE(loop_des) :: T
!> loop descriptors for bosonic grid
      TYPE(loop_des) :: B
!> loop descriptors for fermionic grid
      TYPE(loop_des) :: F

!> global communicator
      TYPE( communic ) :: COMM
!> inter communicator between quadratures
      TYPE( communic ) :: COMM_INTER

   END TYPE imag_grid_handle

!> Interface declarations for module procedures used in derived types

   ABSTRACT INTERFACE 
!> interface for basis function call
      FUNCTION OBJECT_FUNCTION( Z )
         USE prec
# 242

         IMPLICIT NONE
         REAL(qd),INTENT(IN)    :: Z         !< argument
         REAL(qd)               :: OBJECT_FUNCTION 

      END FUNCTION OBJECT_FUNCTION

!> interface for basis function call
      FUNCTION BASIS_FUNCTION( Z, L )
         USE prec
# 254

         IMPLICIT NONE
         REAL(qd),INTENT(IN)    :: Z         !< argument
         REAL(qd),INTENT(IN)    :: L         !< abszissa parameter that is found by fit
         REAL(qd)               :: BASIS_FUNCTION

      END FUNCTION BASIS_FUNCTION

!> interface for initialization of quadrature coefficients
!> for various minimization interval lengths R
      SUBROUTINE INIT_QUAD_SUBROUTINE()
         USE prec
# 268

         IMPLICIT NONE

      END SUBROUTINE INIT_QUAD_SUBROUTINE

!> interface for function that calculates minimax coefficients and alternant
!> from a polynomial that is stored in minimax_dependence
      FUNCTION FIT_COEFF_FUNCTION( N, R )
         USE prec
         IMPLICIT NONE
         INTEGER,INTENT(IN)    :: N         !< requested fit order
         REAL(q),INTENT(IN)    :: R         !< requested interval length
         REAL(q)               :: FIT_COEFF_FUNCTION(4*N+4) !< coefficients and alternant
!> polynomial coefficients of fit describing
!> R-dependence of every quadrature coefficient and alternant
         REAL(q)               :: ALPHA( 24,100, 6 ) 

      END FUNCTION FIT_COEFF_FUNCTION

   END INTERFACE 

   CONTAINS

!***********************************************************************
! DESCRIPTION:
!>calculates the error function for which the minimax solution is found
!> \f$
!>     \eta(z) = f(z) - \sum_{i=1}^N c_{i+N} \phi^2(z,c_i)
!> \f$
!> here \f$ f \f$ is the object function function and \f$ \phi \f$ the
!> basis function.
!> @param[in]  QDT quadrature_handle
!> @param[in]  Z  argument of error function
!***********************************************************************
   FUNCTION ERROR_FUNCTION( QDT , Z ) 
      TYPE( quadrature_handle ) :: QDT
      REAL(qd) :: Z
      REAL(qd) :: ERROR_FUNCTION
! local
      INTEGER :: I
      

      ERROR_FUNCTION = QDT%F( Z )
      DO I = 1 , QDT%N 
         ERROR_FUNCTION = ERROR_FUNCTION - &
           QDT%C( I + QDT%N )*QDT%PHI2( Z, QDT%C(I) )
      ENDDO

      
   END FUNCTION ERROR_FUNCTION

!***********************************************************************
! DESCRIPTION:
!> calculates the frist derivative of the error function for which the
!> minimax solution is found, w.r.t. to \f$ z\f$.
!> \f$
!>     \eta(z) = f(z) - \sum_{i=1}^N c_{i+N} \phi^2(z,c_i)
!> \f$
!> here \f$ f \f$ is the object function function and \f$ \phi \f$ the
!> basis function.
!> @param[in]  QDT quadrature_handle
!> @param[in]  Z  argument of error function
!***********************************************************************
   FUNCTION D1_ERROR_FUNCTION( QDT , Z ) 
      TYPE( quadrature_handle ) :: QDT
      REAL(qd) :: Z
      REAL(qd) :: D1_ERROR_FUNCTION
! local
      INTEGER :: I
      

      D1_ERROR_FUNCTION = QDT%D_F_DZ( Z )
      DO I = 1 , QDT%N 
         D1_ERROR_FUNCTION = D1_ERROR_FUNCTION - &
           QDT%C( I + QDT%N )*QDT%D_PHI2_DZ( Z, QDT%C(I) )
      ENDDO

      
   END FUNCTION D1_ERROR_FUNCTION

!***********************************************************************
! DESCRIPTION:
!> calculates the second derivative of the error function for which the
!> minimax solution is found, w.r.t. to \f$ z\f$.
!> \f$
!>     \eta(z) = f(z) - \sum_{i=1}^N c_{i+N} \phi^2(z,c_i)
!> \f$
!> here \f$ f \f$ is the object function function and \f$ \phi \f$ the
!> basis function.
!> @param[in]  QDT quadrature_handle
!> @param[in]  Z  argument of error function
!***********************************************************************
   FUNCTION D2_ERROR_FUNCTION( QDT , Z ) 
      TYPE( quadrature_handle ) :: QDT
      REAL(qd) :: Z
      REAL(qd) :: D2_ERROR_FUNCTION
! local
      INTEGER :: I
      

      D2_ERROR_FUNCTION = QDT%D2_F_DZ2( Z )

      DO I = 1 , QDT%N 
         D2_ERROR_FUNCTION = D2_ERROR_FUNCTION - &
           QDT%C( I + QDT%N )*QDT%D2_PHI2_DZ2( Z, QDT%C(I) )
      ENDDO

      
   END FUNCTION D2_ERROR_FUNCTION

!***********************************************************************
! DESCRIPTION:
!> calculates the gradient of the error function for which the
!> minimax solution is found, w.r.t. to \f$ \vec c\f$.
!> \f$
!>     \eta(z) = f(z) - \sum_{i=1}^N c_{i+N} \phi^2(z,c_i)
!> \f$
!> here \f$ f \f$ is the object function function and \f$ \phi \f$ the
!> basis function.
!> @param[in]  QDT quadrature_handle
!> @param[in]  Z  argument of error function
!***********************************************************************
   FUNCTION GRAD_ERROR_FUNCTION( QDT , Z, J ) 
      TYPE( quadrature_handle ) :: QDT
      REAL(qd) :: Z
      REAL(qd) :: GRAD_ERROR_FUNCTION
      INTEGER :: J
      

      IF ( J <= QDT%N ) THEN
         GRAD_ERROR_FUNCTION = &
           -QDT%C( J + QDT%N )*QDT%D_PHI2_DL( Z, QDT%C(J) )
      ELSE
         GRAD_ERROR_FUNCTION = -QDT%PHI2( Z, QDT%C( J-QDT%N ) )
      ENDIF 

      
   END FUNCTION GRAD_ERROR_FUNCTION

!*******************************************************************************
!> helper routine, prints quadrature coefficients to file quad_IT.dat
!*******************************************************************************
   SUBROUTINE PRINT_QUADRATURE( QUADRATURE, IO, IT, FNAME )
      USE base
      USE prec
      TYPE( quadrature_handle )  :: QUADRATURE
      TYPE( in_struct )           :: IO
! local
      INTEGER,INTENT(IN),OPTIONAL :: IT    ! add this number to file name
      CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: FNAME ! changes file name
! local
      INTEGER                     :: I
      CHARACTER(LEN=4)            :: APP 

      

      IF ( IO%IU0 >= 0 ) THEN
      IF (PRESENT(IT)) THEN
         WRITE(APP,'(I4)')IT
         APP=ADJUSTL(APP)
         IF ( PRESENT( FNAME ) ) THEN
            OPEN(UNIT=90,FILE=FNAME//TRIM(APP)//'.dat',ACTION='WRITE',STATUS='REPLACE')
         ELSE
            OPEN(UNIT=90,FILE='quad'//TRIM(APP)//'.dat',ACTION='WRITE',STATUS='REPLACE')
         ENDIF
      ELSE
         OPEN(UNIT=90,FILE='quad.dat',ACTION='WRITE',STATUS='REPLACE')
      ENDIF
      WRITE(90,'(2I4,3E13.6)')QUADRATURE%ITYPE,QUADRATURE%N,&
         QUADRATURE%A, QUADRATURE%B, QUADRATURE%SCALING
   
!choose abscissas from coefficients
      DO I = 1, 2*QUADRATURE%N+1
         WRITE(90,'(2F45.32)')QUADRATURE%C(I),QUADRATURE%X0(I)
      ENDDO

      CLOSE(90)
      ENDIF
      
      RETURN
   END SUBROUTINE PRINT_QUADRATURE

!*******************************************************************************
!> helper routine, dumps quadrature coefficients to file
!*******************************************************************************
   SUBROUTINE DUMP_QUADRATURE( QUADRATURE, IT, FNAME )
      USE prec
      USE base
      TYPE( quadrature_handle )  :: QUADRATURE
! local
      INTEGER                    :: IT    
      CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: FNAME ! changes file name
! local
      INTEGER                     :: I
      INTEGER                     :: IUNIT
      CHARACTER(LEN=4)            :: APP 

      

      WRITE(APP,'(I4)')IT
      APP=ADJUSTL(APP)
      IF ( PRESENT( FNAME ) ) THEN
         OPEN(UNIT=IT,FILE=FNAME//TRIM(APP)//'.dat',ACTION='WRITE',STATUS='REPLACE')
      ELSE
         OPEN(UNIT=IT,FILE='quad'//TRIM(APP)//'.dat',ACTION='WRITE',STATUS='REPLACE')
      ENDIF
      WRITE(IT,'(2I4,3E13.6)')QUADRATURE%ITYPE,QUADRATURE%N,&
         QUADRATURE%A, QUADRATURE%B, QUADRATURE%SCALING
   
!choose abscissas from coefficients
      DO I = 1, 2*QUADRATURE%N+1
         WRITE(IT,'(2F45.32)')QUADRATURE%C(I),QUADRATURE%X0(I)
      ENDDO

      CLOSE(IT)
      
      RETURN
   END SUBROUTINE DUMP_QUADRATURE

!*******************************************************************************
!>helper routine, dumps grids to file
!*******************************************************************************
   SUBROUTINE DUMP_IGRID( IGRID, FNAME, INODE )
      USE prec
      TYPE( imag_grid_handle )  :: IGRID
! local
      CHARACTER(LEN=*),INTENT(IN) :: FNAME ! changes file name
      INTEGER                     :: INODE
! local
      INTEGER                     :: I,J
      INTEGER                     :: IT

      IF ( INODE == 1 ) THEN
         IT=100
         OPEN(UNIT=IT,FILE=FNAME//'.dat',ACTION='WRITE',STATUS='REPLACE')
         WRITE(IT,'(I4,L4,6E13.6)')IGRID%NOMEGA,IGRID%LFINITE_TEMPERATURE,&
            IGRID%BETA,IGRID%X1,IGRID%X2,IGRID%R,IGRID%FACTOR,IGRID%GM_INT_ERROR
      
         IF ( ASSOCIATED( IGRID%TAU ) ) THEN
!dump tau grid
         WRITE(IT,'(A,E13.6)')'# TAU-GRID, Error=',IGRID%TAU_ERROR
         DO I = 1, IGRID%NOMEGA
            WRITE(IT,'(2F25.16)')IGRID%TAU(I),IGRID%TAU_WEIGHT(I)
         ENDDO
         ENDIF

!dump BOS_RE grid
         IF ( ASSOCIATED( IGRID%BOS_RE ) ) THEN
         WRITE(IT,'(A,E13.6)')'# BOS_RE-GRID, Error=',IGRID%BOS_RE_ERROR
         DO I = 1, IGRID%NOMEGA
            WRITE(IT,'(2F25.16)')IGRID%BOS_RE(I),IGRID%BOS_RE_WEIGHT(I)
         ENDDO
         ENDIF
         IF ( ASSOCIATED( IGRID%TO_BOS_RE ) ) THEN
!TO_BOS_RE transformation
         IF ( IGRID%NOMEGA>6 ) THEN
            WRITE(IT,'(A)')'# TO_BOS_RE and TO_BOS_RE_CONJG transformation: '
            DO I = 1, IGRID%NOMEGA
            DO J = 1, IGRID%NOMEGA
               WRITE(IT,'(2I3,2F25.16)')I,J,IGRID%TO_BOS_RE(I,J),IGRID%TO_BOS_RE_CONJG(I,J)
            ENDDO
            ENDDO
         ELSE
            WRITE(IT,'(A)')'# TO_BOS_RE transformation: '
            DO I = 1, IGRID%NOMEGA
               WRITE(IT,'(6F25.16)')IGRID%TO_BOS_RE(I,:)
            ENDDO

            WRITE(IT,'(A)')'# TO_BOS_RE_CONJ transformation: '
            DO I = 1, IGRID%NOMEGA
               WRITE(IT,'(6F25.16)')IGRID%TO_BOS_RE_CONJG(I,:)
            ENDDO
         ENDIF
         WRITE(IT,'(A)')'# TO_BOS_RE_ERROR and TO_BOS_RE_CONJG_ERROR:'
         DO I = 1, IGRID%NOMEGA
            WRITE(IT,'(I3,2F25.16)')I,IGRID%TO_BOS_RE_ERROR(I),IGRID%TO_BOS_RE_CONJG_ERROR(I)
         ENDDO
         ENDIF

         IF ( ASSOCIATED( IGRID%TO_BOS_IM ) ) THEN
!dump BOS_IM grid
         WRITE(IT,'(A,E13.6)')'# BOS_IM-GRID, Error=',IGRID%BOS_IM_ERROR
         DO I = 1, IGRID%NOMEGA
            WRITE(IT,'(2F25.16)')IGRID%BOS_IM(I),IGRID%BOS_IM_WEIGHT(I)
         ENDDO
!TO_BOS_IM transformation
         IF ( IGRID%NOMEGA>6 ) THEN
            WRITE(IT,'(A)')'# TO_BOS_IM and TO_BOS_IM_CONJG transformation: '
            DO I = 1, IGRID%NOMEGA
            DO J = 1, IGRID%NOMEGA
               WRITE(IT,'(2I3,2F25.16)')I,J,IGRID%TO_BOS_IM(I,J),IGRID%TO_BOS_IM_CONJG(I,J)
            ENDDO
            ENDDO
         ELSE
            WRITE(IT,'(A)')'# TO_BOS_IM transformation: '
            DO I = 1, IGRID%NOMEGA
               WRITE(IT,'(6F25.16)')IGRID%TO_BOS_IM(I,:)
            ENDDO

            WRITE(IT,'(A)')'# TO_BOS_IM_CONJ transformation: '
            DO I = 1, IGRID%NOMEGA
               WRITE(IT,'(6F25.16)')IGRID%TO_BOS_IM_CONJG(I,:)
            ENDDO
         ENDIF
         WRITE(IT,'(A)')'# TO_BOS_IM_ERROR and TO_BOS_IM_CONJG_ERROR:'
         DO I = 1, IGRID%NOMEGA
            WRITE(IT,'(I3,2F25.16)')I,IGRID%TO_BOS_IM_ERROR(I),IGRID%TO_BOS_IM_CONJG_ERROR(I)
         ENDDO
         ENDIF

         IF ( ASSOCIATED( IGRID%TO_FER_RE ) ) THEN
!dump FER_RE grid
         WRITE(IT,'(A,E13.6)')'# FER_RE-GRID, Error=',IGRID%FER_RE_ERROR
         DO I = 1, IGRID%NOMEGA
            WRITE(IT,'(2F25.16)')IGRID%FER_RE(I),IGRID%FER_RE_WEIGHT(I)
         ENDDO
!TO_FER_RE transformation
         IF ( IGRID%NOMEGA>6 ) THEN
            WRITE(IT,'(A)')'# TO_FER_RE and TO_FER_RE_CONJG transformation: '
            DO I = 1, IGRID%NOMEGA
            DO J = 1, IGRID%NOMEGA
               WRITE(IT,'(2I3,2F25.16)')I,J,IGRID%TO_FER_RE(I,J),IGRID%TO_FER_RE_CONJG(I,J)
            ENDDO
            ENDDO
         ELSE
            WRITE(IT,'(A)')'# TO_FER_RE transformation: '
            DO I = 1, IGRID%NOMEGA
               WRITE(IT,'(6F25.16)')IGRID%TO_FER_RE(I,:)
            ENDDO

            WRITE(IT,'(A)')'# TO_FER_RE_CONJ transformation: '
            DO I = 1, IGRID%NOMEGA
               WRITE(IT,'(6F25.16)')IGRID%TO_FER_RE_CONJG(I,:)
            ENDDO
         ENDIF
         WRITE(IT,'(A)')'# TO_FER_RE_ERROR and TO_FER_RE_CONJG_ERROR:'
         DO I = 1, IGRID%NOMEGA
            WRITE(IT,'(I3,2F25.16)')I,IGRID%TO_FER_RE_ERROR(I),IGRID%TO_FER_RE_CONJG_ERROR(I)
         ENDDO
         ENDIF

         IF ( ASSOCIATED( IGRID%TO_FER_IM ) ) THEN
!dump FER_IM grid
         WRITE(IT,'(A,E13.6)')'# FER_IM-GRID, Error=',IGRID%FER_IM_ERROR
         DO I = 1, IGRID%NOMEGA
            WRITE(IT,'(2F25.16)')IGRID%FER_IM(I),IGRID%FER_IM_WEIGHT(I)
         ENDDO
!TO_FER_IM transformation
         IF ( IGRID%NOMEGA>6 ) THEN
            WRITE(IT,'(A)')'# TO_FER_IM and TO_FER_IM_CONJG transformation: '
            DO I = 1, IGRID%NOMEGA
            DO J = 1, IGRID%NOMEGA
               WRITE(IT,'(2I3,2F25.16)')I,J,IGRID%TO_FER_IM(I,J),IGRID%TO_FER_IM_CONJG(I,J)
            ENDDO
            ENDDO
         ELSE
            WRITE(IT,'(A)')'# TO_FER_IM transformation: '
            DO I = 1, IGRID%NOMEGA
               WRITE(IT,'(6F25.16)')IGRID%TO_FER_IM(I,:)
            ENDDO

            WRITE(IT,'(A)')'# TO_FER_IM_CONJ transformation: '
            DO I = 1, IGRID%NOMEGA
               WRITE(IT,'(6F25.16)')IGRID%TO_FER_IM_CONJG(I,:)
            ENDDO
         ENDIF
         WRITE(IT,'(A)')'# TO_FER_IM_ERROR and TO_FER_IM_CONJG_ERROR:'
         DO I = 1, IGRID%NOMEGA
            WRITE(IT,'(I3,2F25.16)')I,IGRID%TO_FER_IM_ERROR(I),IGRID%TO_FER_IM_CONJG_ERROR(I)
         ENDDO
         ENDIF

      ENDIF

      CLOSE(IT)
   END SUBROUTINE DUMP_IGRID
   
!*******************************************************************************
!>helper routine, prints the nonlinear error function and its first two
! derivatives w.r.t. argument X
!*******************************************************************************
   SUBROUTINE PRINT_ERROR_FUNCTION( QUADRATURE, IO, IT, FNAME )
      USE base
      USE prec
      TYPE( quadrature_handle )  :: QUADRATURE
      TYPE( in_struct )           :: IO
! local
      INTEGER,INTENT(IN),OPTIONAL :: IT    ! add this number to file name
      CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: FNAME ! changes file name
! local
      INTEGER                     :: I,J,N
      REAL(qd)                     :: X, X2, X1
      REAL(qd)                     :: Y,D1Y,D2Y
      CHARACTER(LEN=4)            :: APP 
      INTEGER, PARAMETER          :: IMAX = 100

      
   
      IF ( IO%IU0 >= 0 ) THEN
      IF (PRESENT(IT)) THEN
         WRITE(APP,'(I4)')IT
         APP=ADJUSTL(APP)
         IF ( PRESENT( FNAME ) ) THEN
            OPEN(UNIT=91,FILE=FNAME//TRIM(APP)//'.dat',ACTION='WRITE',STATUS='REPLACE')
         ELSE
            OPEN(UNIT=91,FILE='erf'//TRIM(APP)//'.dat',ACTION='WRITE',STATUS='REPLACE')
         ENDIF
      ELSE
         OPEN(UNIT=91,FILE='erf.dat',ACTION='WRITE',STATUS='REPLACE')
      ENDIF
      WRITE(91,*)'#',QUADRATURE%A,QUADRATURE%B
      N = QUADRATURE%N
      DO J = 1, 2*N+1
         WRITE(91,'(A,2F50.32)')'# ',QUADRATURE%C(J),QUADRATURE%X0(J)
      ENDDO
   
!choose abscissas from coefficients
      X1 = QUADRATURE%A
      DO J = 1, N
         X2 = X1+QUADRATURE%C(J)*QUADRATURE%B/QUADRATURE%C(N)
         DO I=1,IMAX  !points in between two abscissas
!equidistant abscissas between coefficients
            X=X1+((I-1)*(X2-X1))/IMAX
            IF( X > QUADRATURE%B) THEN
               X = QUADRATURE%B
               Y = ERROR_FUNCTION( QUADRATURE, X )
               D1Y = D1_ERROR_FUNCTION( QUADRATURE, X )
               D2Y = D2_ERROR_FUNCTION( QUADRATURE, X )
               WRITE(91,'(4E40.32)')X,Y,D1Y,D2Y
               EXIT    
            ENDIF
            Y = ERROR_FUNCTION( QUADRATURE, X )
            D1Y = D1_ERROR_FUNCTION( QUADRATURE, X )
            D2Y = D2_ERROR_FUNCTION( QUADRATURE, X )
            WRITE(91,'(4E40.32)')X,Y,D1Y,D2Y
         ENDDO
         X1 = X2
      ENDDO

      CLOSE(91)
      ENDIF 
      
   END SUBROUTINE PRINT_ERROR_FUNCTION
   
!*******************************************************************************
!>helper routine, prints the nonlinear error function and its first two
! derivatives w.r.t. argument X
!*******************************************************************************
   SUBROUTINE PRINT_TRANS_ERROR_FUNCTION( QUAD_A, QUAD_B, ALPHA, OMEGA, &
      LCONJG, NODE_ME, IT )
      USE prec
      USE base
      TYPE( quadrature_handle )  :: QUAD_A
      TYPE( quadrature_handle )  :: QUAD_B
      REAL(qd)                     :: ALPHA(:)
      REAL(qd)                     :: OMEGA
      INTEGER                     :: NODE_ME 
      LOGICAL                     :: LCONJG
      INTEGER,INTENT(IN),OPTIONAL :: IT    ! add this number to file name
! local
      INTEGER                     :: I, J, K
      REAL(qd)                     :: X, X2, X1
      REAL(qd)                     :: Y
      CHARACTER(LEN=4)            :: APP 
      INTEGER, PARAMETER          :: IMAX = 100

      
   
      IF (PRESENT(IT)) THEN
         WRITE(APP,'(I4)')IT
         APP=ADJUSTL(APP)
         IF ( LCONJG ) THEN
            OPEN(UNIT=100+NODE_ME,FILE='erf_w'//TRIM(APP)//'_conjg.dat',ACTION='WRITE',STATUS='REPLACE')
         ELSE
            OPEN(UNIT=100+NODE_ME,FILE='erf_w'//TRIM(APP)//'.dat',ACTION='WRITE',STATUS='REPLACE')
         ENDIF
      ELSE
         IF ( LCONJG ) THEN
            OPEN(UNIT=100+NODE_ME,FILE='erf_w_conjg.dat',ACTION='WRITE',STATUS='REPLACE')
         ELSE
            OPEN(UNIT=100+NODE_ME,FILE='erf_w.dat',ACTION='WRITE',STATUS='REPLACE')
         ENDIF
      ENDIF
      WRITE(100+NODE_ME,'(A,5E14.6)')'#',QUAD_A%A,QUAD_A%B,&
        QUAD_B%A,QUAD_B%B, OMEGA
   
!choose abscissas from coefficients
      X1 = QUAD_B%A
      DO J = 1, QUAD_B%N
         IF ( QUAD_B%N > 1 .AND. ABS( QUAD_B%C(QUAD_B%N) ) > REAL(1.E-16,KIND=qd) ) THEN
            X2 = X1+QUAD_B%C(J)*QUAD_B%B/QUAD_B%C(QUAD_B%N)
         ELSE
            X2 = X1+QUAD_B%B
         ENDIF
         DO I=1,IMAX  !points in between two abscissas
!equidistant abscissas between coefficients
            X=X1+((I-1)*(X2-X1))/IMAX
            IF( X > QUAD_B%B) THEN
               X = QUAD_B%B
               IF ( LCONJG ) THEN
                  Y = QUAD_B%PHI_CONJG( X, OMEGA )
                  DO K = 1, QUAD_A%N
                     Y = Y - ALPHA( K )*QUAD_B%PSI_CONJG( X, QUAD_A%C( K ) )
                  ENDDO
               ELSE
                  Y = QUAD_B%PHI( X, OMEGA )
                  DO K = 1, QUAD_A%N
                     Y = Y - ALPHA( K )*QUAD_B%PSI( X, QUAD_A%C( K ) )
                  ENDDO
               ENDIF
               WRITE(100+NODE_ME,'(2E40.32)')X,Y
               EXIT    
            ENDIF
            IF ( LCONJG ) THEN
               Y = QUAD_B%PHI_CONJG( X, OMEGA )
               DO K = 1, QUAD_A%N
                  Y = Y - ALPHA( K )*QUAD_B%PSI_CONJG( X, QUAD_A%C( K ) )
               ENDDO
            ELSE
               Y = QUAD_B%PHI( X, OMEGA )
               DO K = 1, QUAD_A%N
                  Y = Y - ALPHA( K )*QUAD_B%PSI( X, QUAD_A%C( K ) )
               ENDDO
            ENDIF
            WRITE(100+NODE_ME,'(2E40.32)')X,Y
         ENDDO
         X1 = X2
      ENDDO
      CLOSE(100+NODE_ME)
      
   END SUBROUTINE PRINT_TRANS_ERROR_FUNCTION

!*******************************************************************************
!>helper routine, prints the nonlinear error function and its first two
! derivatives w.r.t. argument X
!*******************************************************************************
   SUBROUTINE PRINT_TRANSFORMATION_ERROR( OMEGA, ERRORS, N, LCONJG, NODE_ME, IT )
      USE prec
      USE base
      REAL(q)                     :: OMEGA(:)
      REAL(q)                     :: ERRORS(:)
      INTEGER                     :: N ! quadrature points come to separate file
      LOGICAL                     :: LCONJG
      INTEGER                     :: NODE_ME 
      INTEGER,INTENT(IN),OPTIONAL :: IT    ! add this number to file name
! local
      INTEGER                     :: I
      CHARACTER(LEN=4)            :: APP 

      
      
      IF (PRESENT(IT)) THEN
         WRITE(APP,'(I4)')IT
         APP=ADJUSTL(APP)
         IF ( LCONJG ) THEN
            OPEN(UNIT=1000+NODE_ME,FILE='erf_trans'//TRIM(APP)//'_conjg.dat',ACTION='WRITE',STATUS='REPLACE')
            OPEN(UNIT=2000+NODE_ME,FILE='erf_trans_opt'//TRIM(APP)//'_conjg.dat',ACTION='WRITE',STATUS='REPLACE')
         ELSE
            OPEN(UNIT=1000+NODE_ME,FILE='erf_trans'//TRIM(APP)//'.dat',ACTION='WRITE',STATUS='REPLACE')
            OPEN(UNIT=2000+NODE_ME,FILE='erf_trans_opt'//TRIM(APP)//'.dat',ACTION='WRITE',STATUS='REPLACE')
         ENDIF
      ELSE
         IF ( LCONJG ) THEN
            OPEN(UNIT=1000+NODE_ME,FILE='erf_trans_conjg.dat',ACTION='WRITE',STATUS='REPLACE')
            OPEN(UNIT=2000+NODE_ME,FILE='erf_trans_opt_conjg.dat',ACTION='WRITE',STATUS='REPLACE')
         ELSE
            OPEN(UNIT=1000+NODE_ME,FILE='erf_trans.dat',ACTION='WRITE',STATUS='REPLACE')
            OPEN(UNIT=2000+NODE_ME,FILE='erf_trans_opt.dat',ACTION='WRITE',STATUS='REPLACE')
         ENDIF
      ENDIF
!choose abscissas from coefficients
      DO I = 1, SIZE( ERRORS )-N
         WRITE(1000+NODE_ME,'( 2E20.10 )' )OMEGA(I),ERRORS(I)
      ENDDO
      CLOSE(1000+NODE_ME)

!choose abscissas from coefficients
      DO I = SIZE(ERRORS)-N+1, SIZE( ERRORS )
         WRITE(2000+NODE_ME,'( 2E20.10 )' ) OMEGA(I), ERRORS(I)
      ENDDO
      CLOSE(2000+NODE_ME)

      

   END SUBROUTINE PRINT_TRANSFORMATION_ERROR

!*******************************************************************************
!>helper routine, prints Transformation Matrix to file
!*******************************************************************************
   SUBROUTINE PRINT_TRANSFORMATION_MATRIX( QUAD_A, QUAD_B, ALPHA, LCONJG, NODE_ME, IT )
      USE prec
      USE base
      TYPE( quadrature_handle )   :: QUAD_A
      TYPE( quadrature_handle )   :: QUAD_B
      REAL(q)                     :: ALPHA(:,:)
      INTEGER                     :: NODE_ME 
      LOGICAL                     :: LCONJG
      INTEGER,INTENT(IN),OPTIONAL :: IT    ! add this number to file name
! local
      INTEGER                     :: I, J
      CHARACTER(LEN=4)            :: APP 
      INTEGER, PARAMETER          :: IMAX = 100

      
   
      IF (PRESENT(IT)) THEN
         WRITE(APP,'(I4)')IT
         APP=ADJUSTL(APP)
         IF ( LCONJG ) THEN
            OPEN(UNIT=100+NODE_ME,FILE='trafo'//TRIM(APP)//'_conjg.dat',ACTION='WRITE',STATUS='REPLACE')
         ELSE
            OPEN(UNIT=100+NODE_ME,FILE='trafo'//TRIM(APP)//'.dat',ACTION='WRITE',STATUS='REPLACE')
         ENDIF
      ELSE
         IF ( LCONJG ) THEN
            OPEN(UNIT=100+NODE_ME,FILE='trafo_conjg.dat',ACTION='WRITE',STATUS='REPLACE')
         ELSE
            OPEN(UNIT=100+NODE_ME,FILE='trafo.dat',ACTION='WRITE',STATUS='REPLACE')
         ENDIF
      ENDIF
      WRITE(100+NODE_ME,'(A,4E14.6,L4)')'#',QUAD_A%A,QUAD_A%B,&
        QUAD_B%A,QUAD_B%B, LCONJG
      WRITE(100+NODE_ME,'(A)')'#     I      J              ALPHA(I,J)       '
   
!choose abscissas from coefficients
      DO I = 1, QUAD_A%N
      DO J = 1, QUAD_B%N
         WRITE(100+NODE_ME,'(2I7,F25.16)')I,J,ALPHA(I,J)
      ENDDO
      ENDDO
      CLOSE(100+NODE_ME)
      
   END SUBROUTINE PRINT_TRANSFORMATION_MATRIX

# 1425

END MODULE minimax_struct
