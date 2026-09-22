# 1 "minimax_functions1D.F"
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


# 2 "minimax_functions1D.F" 2 
!***********************************************************************
!
! MODULE: minimax_functions
!
!> @author
!> Merzuk Kaltak, VASP Software GmbH
!
! DESCRIPTION:
!> contains all (1._q,0._q)-dimensional functions
!>
!***********************************************************************
MODULE minimax_functions1D
   USE prec
# 20

   IMPLICIT NONE 
!truncation threshold for Taylor expansion for finite temperature error functions
   REAL(qd), PARAMETER :: ACC3=1.E-8_qd   
   REAL(qd) ,PARAMETER :: PIQ=3.1415926535897932384626433832795028841971693993751_qd

! maximum exponent to base 2 is obtained from intrinsic function
! MAXEXPONENT, to base e the maximum exponent is MAXEXPONENT(X) * LN 2
! for QD MAXEXPONENT return usually 16384* LN 2 = 11356
! we use a 1/2.5 part of this (functions are being squared) ->4500
   INTEGER, PARAMETER :: MAXE=100

   CONTAINS

!***********************************************************************
! DESCRIPTION:
!>calculates following function:
!> \f$
!>     \frac{1}{z}
!> \f$
!> @param[in]  Z  (argument)
!***********************************************************************
  FUNCTION INVERSE_Z(Z)
      REAL(qd), INTENT(IN) :: Z
      REAL(qd)             :: INVERSE_Z
       
      IF (ABS(Z)<REAL(1.E-9,KIND=qd)) THEN
         INVERSE_Z=0
      ELSE
         INVERSE_Z=REAL(1,KIND=qd)/Z
      ENDIF
   END FUNCTION INVERSE_Z

!***********************************************************************
! DESCRIPTION:
!>calculates following function:
!> \f$
!>     -\frac{1}{z^2}
!> \f$
!> @param[in]  Z  (argument)
!***********************************************************************
  FUNCTION D1_INVERSE_Z(Z)
      REAL(qd), INTENT(IN) :: Z
      REAL(qd)             :: D1_INVERSE_Z
       
      IF (ABS(Z)<REAL(1.E-9,KIND=qd)) THEN
         D1_INVERSE_Z=0
      ELSE
         D1_INVERSE_Z=-REAL(1,KIND=qd)/(Z*Z)
      ENDIF
   END FUNCTION D1_INVERSE_Z

!***********************************************************************
! DESCRIPTION:
!>calculates following function:
!> \f$
!>     \frac{2}{z^3}
!> \f$
!> @param[in]  Z  (argument)
!***********************************************************************
  FUNCTION D2_INVERSE_Z(Z)
      REAL(qd), INTENT(IN) :: Z
      REAL(qd)             :: D2_INVERSE_Z
       
      IF (ABS(Z)<REAL(1.E-9,KIND=qd)) THEN
         D2_INVERSE_Z=0
      ELSE
         D2_INVERSE_Z=REAL(2,KIND=qd)/(Z**3)
      ENDIF
   END FUNCTION D2_INVERSE_Z

!****************************************************************************
!
! DESCRIPTION:
!>calculates the hyperbolic tangent tanh
!
!> @param[in]  Z  (argument)
!
!****************************************************************************

   FUNCTION TH( Z ) 
      USE prec
# 104

      REAL(qd) :: Z
      REAL(qd) :: TH
   
      
! error is smaller than 7.4*10^{-44}
      IF ( ABS( Z ) >  50) THEN
         TH = SIGN_QD(Z)
      ELSE IF ( ABS( Z) < REAL(1E-8,KIND=qd) )THEN
         TH = Z**2/REAL(2,KIND=qd)-Z**3/REAL(24,KIND=qd)+Z**5/REAL(240,KIND=qd)
      ELSE
         TH = TANH( Z ) 
      ENDIF
      
   END FUNCTION TH

!****************************************************************************
!
! DESCRIPTION:
!>calculates the derivative of hyperbolic tangent tanh
!
!> @param[in]  Z  (argument)
!
!****************************************************************************

   FUNCTION DTH( Z ) 
      USE prec
# 133

      REAL(qd) :: Z
      REAL(qd) :: DTH
   
      
! error is smaller than 7.4*10^{-44}
      IF ( ABS( Z ) >  50) THEN
         DTH = 0
      ELSE
         DTH = REAL(1,KIND=qd)-TH(Z)**2
      ENDIF
      
   END FUNCTION DTH

!****************************************************************************
!
! DESCRIPTION:
!>calculates:
!> \f$
!>  \frac12\tanh\frac{z}{2}
!> \f$
!
!> @param[in]  Z  (argument)
!
!****************************************************************************

   FUNCTION F_TANH( Z ) 
      USE prec
      REAL(qd), INTENT(IN) :: Z
      REAL(qd) :: F_TANH
      
      F_TANH = TH( Z/REAL(2,KIND=qd) ) /REAL(2,KIND=qd)
      
   END FUNCTION F_TANH

!****************************************************************************
!
! DESCRIPTION:
!>calculates first derivative of:
!> \f$
!>  \frac12\tanh\frac{z}{2}
!> \f$
!
!> @param[in]  Z  (argument)
!
!****************************************************************************

   FUNCTION D1_F_TANH( Z ) 
      USE prec
      REAL(qd), INTENT(IN) :: Z
      REAL(qd) :: D1_F_TANH
      
      D1_F_TANH = (REAL(1,KIND=qd)-REAL(4,KIND=qd)*F_TANH(Z)**2)/REAL(4,KIND=qd)
      
   END FUNCTION D1_F_TANH

!****************************************************************************
!
! DESCRIPTION:
!>calculates second derivative of:
!> \f$
!>  \frac12\tanh\frac{z}{2}
!> \f$
!
!> @param[in]  Z  (argument)
!
!****************************************************************************

   FUNCTION D2_F_TANH( Z ) 
      USE prec
      REAL(qd), INTENT(IN) :: Z
      REAL(qd) :: T
      REAL(qd) :: D2_F_TANH
      
      T = F_TANH(Z)
      D2_F_TANH = -REAL(2,KIND=qd)*T*(REAL(1,KIND=qd)-REAL(4,KIND=qd)*T**2)/REAL(4,KIND=qd)
      
   END FUNCTION D2_F_TANH

!***********************************************************************
! DESCRIPTION:
!>calculates following function:
!> \f$
!>     \frac{\tanh\frac{z}2}{z}
!> \f$
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION F_THZ_OVER_Z( Z ) 
      REAL(qd)               :: F_THZ_OVER_Z
      REAL(qd) , INTENT(IN)  :: Z 
      

! use a taylor expansion for small arguments
      IF ( ABS( Z ) < REAL(1E-8,KIND=qd) ) THEN
         F_THZ_OVER_Z=REAL(0.5,KIND=qd)-Z**2/REAL(24,KIND=qd)+Z**4/REAL(240,KIND=qd)&
                    -REAL(17,KIND=qd)*Z**6/REAL(40320,KIND=qd)
      ELSE
         F_THZ_OVER_Z =TH( Z/REAL(2,KIND=qd) )/Z
      ENDIF
      
   END FUNCTION F_THZ_OVER_Z

!***********************************************************************
! DESCRIPTION:
!>calculates first derivative of following function:
!> \f$
!>     \frac{\tanh\frac{z}2}{z}
!> \f$
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION D1_THZ_OVER_Z( Z ) 
      REAL(qd)               :: D1_THZ_OVER_Z
      REAL(qd) , INTENT(IN)  :: Z 
      REAL(qd)               :: THZ
      

! use a taylor expansion for small arguments
      IF ( ABS( Z ) < REAL(1E-8,KIND=qd) ) THEN
         D1_THZ_OVER_Z= -Z/REAL(12,KIND=qd)+Z**3/REAL(60,KIND=qd)-REAL(17,KIND=qd)*Z**5/REAL(6720,KIND=qd)
      ELSE
         THZ = TH(Z/REAL(2,KIND=qd))
         D1_THZ_OVER_Z= (REAL(1,KIND=qd)-THZ**2)/(REAL(2,KIND=qd)*Z) - THZ/Z**2
      ENDIF
      
   END FUNCTION D1_THZ_OVER_Z

!***********************************************************************
! DESCRIPTION:
!>calculates second derivative of following function:
!> \f$
!>     \frac{\tanh\frac{z}2}{z}
!> \f$
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION D2_THZ_OVER_Z( Z ) 
      REAL(qd)               :: D2_THZ_OVER_Z
      REAL(qd) , INTENT(IN)  :: Z 
      REAL(qd)               :: THZ
      

! use a taylor expansion for small arguments
      IF ( ABS( Z ) < REAL(1E-8,KIND=qd) ) THEN
         D2_THZ_OVER_Z= -REAL(1,KIND=qd)/REAL(12,KIND=qd)                    &
                        +Z**2/REAL(20,KIND=qd)                       &
                        -REAL(17,KIND=qd)*Z**4/REAL(1344,KIND=qd)            &
                        +REAL(31,KIND=qd)*Z**6/REAL(12960,KIND=qd)
      ELSE
         THZ = TH(Z/REAL(2,KIND=qd))
         D2_THZ_OVER_Z= REAL(2,KIND=qd)*THZ/Z**3                   &
                      - (REAL(1,KIND=qd)-THZ**2)/Z**2              &
                      - THZ*(REAL(1,KIND=qd)-THZ**2)/(REAL(2,KIND=qd)*Z)   
      ENDIF
      
   END FUNCTION D2_THZ_OVER_Z

!***********************************************************************
! DESCRIPTION:
!>calculates following function:
!> \f$
!>     \frac{1}{8}\left(1-\tanh^2\frac{z}{2}+\frac{2}{z}\tanh\frac{z}{2}\right)
!> \f$
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION F_BETA_LRS( Z ) 
      REAL(qd)               :: F_BETA_LRS
      REAL(qd) , INTENT(IN)  :: Z 
      REAL(qd)               :: ARG 
      

      ARG = Z/REAL(2,KIND=qd)
! use a taylor expansion for small arguments
      IF ( ABS( Z ) < REAL(1E-8,KIND=qd) ) THEN
         F_BETA_LRS=REAL(0.25,KIND=qd)-Z**2/REAL(24,KIND=qd)+Z**4/REAL(160,KIND=qd)&
                       -REAL(17,KIND=qd)*Z**6/REAL(20160,KIND=qd)
      ELSE
         F_BETA_LRS =( REAL(0.125,KIND=qd)-REAL(0.125,KIND=qd)*TH( ARG )**2 )&
                        +( REAL(2,KIND=qd)*TH( ARG )/(REAL(8,KIND=qd)*Z) )
      ENDIF
      
   END FUNCTION F_BETA_LRS

!***********************************************************************
! DESCRIPTION:
!>calculates following function:
!> \f$
!>     \frac{1}{8}\left(1-\tanh^2\frac{z}{2}+\frac{2}{z}\tanh\frac{z}{2}\right)
!>     -\frac{1}{z^2}\tanh^2\frac{z}{2}
!> \f$
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION F_BETA_LRS_REG( Z ) 
      REAL(qd)               :: F_BETA_LRS_REG
      REAL(qd) , INTENT(IN)  :: Z 
      REAL(qd)               :: ARG 
      

      ARG = Z/REAL(2,KIND=qd)
! use a taylor expansion for small arguments
      IF ( ABS( Z ) < REAL(1E-8,KIND=qd) ) THEN
         F_BETA_LRS_REG=Z**4/REAL(2880,KIND=qd)-Z**6/REAL(13440,KIND=qd)
      ELSE
         ARG = TH( ARG ) 
         F_BETA_LRS_REG =( REAL(0.125,KIND=qd)-REAL(0.125,KIND=qd)*ARG**2 )&
                        +( REAL(2,KIND=qd)*ARG/(REAL(8,KIND=qd)*Z) ) &
                        -( ARG/Z )**2
      ENDIF
      
   END FUNCTION F_BETA_LRS_REG

!***********************************************************************
! DESCRIPTION:
!>calculates first two derivatives of following function:
!> \f$
!>     \frac{1}{8}\left(1-\tanh^2\frac{z}{2}+\frac{2}{z}\tanh\frac{z}{2}\right)
!> \f$
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION F_BETA_LRS_DER( Z ) 
      REAL(qd)               :: F_BETA_LRS_DER(2)
      REAL(qd) , INTENT(IN)  :: Z 
      REAL(qd)               :: ZP(2),T
      

      T = TH( Z/REAL(2,KIND=qd) )
! For small arguments use taylor expansion
      IF ( ABS( Z ) < ACC3 ) THEN
! for |Z|<1.E-6 the error of expansion and exact value is < 10^-39

! first derivative  w.r.t. Z
         ZP(1)= -Z/REAL(12,KIND=qd) + Z**3/REAL(40,KIND=qd) - (REAL(17,KIND=qd)*Z**5)/REAL(3360,KIND=qd)
               
! second derivative  w.r.t. Z
         ZP(2)= -REAL(1,KIND=qd)/REAL(12,KIND=qd) + REAL(3,KIND=qd)*Z**2/REAL(40,KIND=qd) - (REAL(17,KIND=qd)*Z**4)/REAL(672,KIND=qd)
      ELSE
! write functions in terms of Tanh[ Z/2 ]
! first derivative  w.r.t. Z
         ZP(1)= ( REAL(1,KIND=qd)/( REAL(8,KIND=qd)*Z  ) ) &
              - ( T/REAL(8,KIND=qd) )          &
              - T*(REAL(1,KIND=qd)/( REAL(2,KIND=qd)*Z) )**2  &
              - ( T**2/( REAL(8,KIND=qd)*Z ) ) &
              + ( T**3/( REAL(8,KIND=qd) ) )

! second derivative  w.r.t. Z
         ZP(2)= -REAL(0.0625,KIND=qd)                     &
                -( REAL(1,KIND=qd)/(REAL(2,KIND=qd)*Z))**2        &
                +( T/REAL(2,KIND=qd) )*( REAL(1,KIND=qd)/Z )**3   &
                -( T/(REAL(8,KIND=qd)*Z) )                &
                +( T/REAL(2,KIND=qd) )**2                 &
                +( T/(REAL(2,KIND=qd)*Z) )**2             &
                +( T**3/(REAL(8,KIND=qd)*Z) )             &
                -( (REAL(3,KIND=qd)*T**4)/REAL(16,KIND=qd) )     
      ENDIF 

      F_BETA_LRS_DER(1) = ZP(1)
      F_BETA_LRS_DER(2) = ZP(2)

      
   END FUNCTION F_BETA_LRS_DER

!***********************************************************************
! DESCRIPTION:
!>calculates first derivative of following function:
!> \f$
!>     \frac{1}{8}(1-\tanh^2\frac{z}{2}+\frac{2}{z}\tanh\frac{z}{2})
!> \f$
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION D1_F_BETA_LRS( Z ) 
      REAL(qd)               :: D1_F_BETA_LRS
      REAL(qd) , INTENT(IN)  :: Z 
      REAL(qd)               :: T,DT
      

      T = TH( Z/REAL(2,KIND=qd) )
! For small arguments use taylor expansion
      IF ( ABS( Z ) < ACC3 ) THEN
! for |Z|<1.E-6 the error of expansion and exact value is < 10^-39

! first derivative  w.r.t. Z
         D1_F_BETA_LRS= -Z/REAL(12,KIND=qd) + Z**3/REAL(40,KIND=qd) - (REAL(17,KIND=qd)*Z**5)/REAL(3360,KIND=qd)
      ELSE
! write functions in terms of Tanh[ Z/2 ]
! first derivative  w.r.t. Z
         D1_F_BETA_LRS= ( REAL(1,KIND=qd)/( REAL(8,KIND=qd)*Z  ) ) &
              - ( T/REAL(8,KIND=qd) )          &
              - T*(REAL(1,KIND=qd)/( REAL(2,KIND=qd)*Z) )**2  &
              - ( T**2/( REAL(8,KIND=qd)*Z ) ) &
              + ( T**3/( REAL(8,KIND=qd) ) )
      ENDIF 

      
   END FUNCTION D1_F_BETA_LRS

!***********************************************************************
! DESCRIPTION:
!>calculates second derivative of following function:
!> \f$
!>     \frac{1}{8}\left(1-\tanh^2\frac{z}{2}+\frac{2}{z}\tanh\frac{z}{2}\right)
!> \f$
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION D2_F_BETA_LRS( Z ) 
      REAL(qd)               :: D2_F_BETA_LRS
      REAL(qd) , INTENT(IN)  :: Z 
      REAL(qd)               :: T
      

      T = TH( Z/REAL(2,KIND=qd) )
! For small arguments use taylor expansion
      IF ( ABS( Z ) < ACC3 ) THEN
! for |Z|<1.E-6 the error of expansion and exact value is < 10^-39
! second derivative  w.r.t. Z
         D2_F_BETA_LRS= -REAL(1,KIND=qd)/REAL(12,KIND=qd) + REAL(3,KIND=qd)*Z**2/REAL(40,KIND=qd) - (REAL(17,KIND=qd)*Z**4)/REAL(672,KIND=qd)
      ELSE
! second derivative  w.r.t. Z
         D2_F_BETA_LRS= -REAL(0.0625,KIND=qd)             &
                -( REAL(1,KIND=qd)/(REAL(2,KIND=qd)*Z))**2        &
                +( T/REAL(2,KIND=qd) )*( REAL(1,KIND=qd)/Z )**3   &
                -( T/(REAL(8,KIND=qd)*Z) )                &
                +( T/REAL(2,KIND=qd) )**2                 &
                +( T/(REAL(2,KIND=qd)*Z) )**2             &
                +( T**3/(REAL(8,KIND=qd)*Z) )             &
                -( (REAL(3,KIND=qd)*T**4)/REAL(16,KIND=qd) )     
      ENDIF 

      
   END FUNCTION D2_F_BETA_LRS

!***********************************************************************
! DESCRIPTION:
!>calculates frst two dervatives of following function:
!> \f$
!>     \frac{1}{8}\left(1-\tanh^2\frac{z}{2}+\frac{2}{z}\tanh\frac{z}{2}\right)
!>     -\frac{1}{z^2}\tanh^2\frac{z}{2}
!> \f$
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION F_BETA_LRS_REG_DER( Z ) 
      REAL(qd)               :: F_BETA_LRS_REG_DER(2)
      REAL(qd) , INTENT(IN)  :: Z 
      REAL(qd)               :: ZP(2),T,TZ
      

! For small arguments use taylor expansion
      IF ( ABS( Z ) < ACC3 ) THEN
! first derivative  w.r.t. Z
         ZP(1)= Z**3/REAL(720,KIND=qd)-Z**5/REAL(2240,KIND=qd)
! second derivative  w.r.t. Z
         ZP(2)=Z**2/REAL(240,KIND=qd)-Z**4/REAL(448,KIND=qd)+REAL(7,KIND=qd)*Z**6/REAL(10800,KIND=qd)
      ELSE
         T = TH( Z/REAL(2,KIND=qd) )
         TZ = T/Z
! write functions in terms of Tanh[ Z/2 ]
! first derivative  w.r.t. Z
         ZP(1)=                                   &
              - ( T/REAL(8,KIND=qd) )                   &
              + ( T**3/REAL(8,KIND=qd) )                &
              + ( REAL(2,KIND=qd)*TZ**2/Z  )              &
              - ( REAL(5,KIND=qd)/REAL(4,KIND=qd)*TZ/Z )          &
              + ( T*TZ**2 )                       &
              + ( REAL(1,KIND=qd)/( REAL(8,KIND=qd)*Z  ) )    &
              - ( T*TZ/REAL(8,KIND=qd) )    

! second derivative  w.r.t. Z
         ZP(2)=                                   &
              -  REAL(0.0625,KIND=qd)                     &
              + ( T**2/REAL(4,KIND=qd) )                &
              - (REAL(3,KIND=qd)/REAL(16,KIND=qd))*T**4           &
              - (REAL(6,KIND=qd)*(TZ**2/Z**2 ))           &
              + (REAL(9,KIND=qd)*TZ/(REAL(2,KIND=qd)*Z**2 ))      &
              - (REAL(4,KIND=qd)*TZ**3)                   &
              - (REAL(3,KIND=qd)/(REAL(4,KIND=qd)*Z**2 ))         &
              + (REAL(9,KIND=qd)/REAL(4,KIND=qd)*TZ**2)           &
              - (REAL(3,KIND=qd)/REAL(2,KIND=qd)*T**2*TZ**2)      &
              - (TZ/REAL(8,KIND=qd) )                     &
              + (TZ*T**2/REAL(8,KIND=qd) )      
      ENDIF 

      F_BETA_LRS_REG_DER(1) = ZP(1)
      F_BETA_LRS_REG_DER(2) = ZP(2)
      
   END FUNCTION F_BETA_LRS_REG_DER

!***********************************************************************
! DESCRIPTION:
!>calculates first dervatives of following function:
!> \f$
!>     \frac{1}{8}\left(1-\tanh^2\frac{z}{2}+\frac{2}{z}\tanh\frac{z}{2}\right)
!>     -\frac{1}{z^2}\tanh^2\frac{z}{2}
!> \f$
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION D1_F_BETA_LRS_REG( Z ) 
      REAL(qd)               :: D1_F_BETA_LRS_REG
      REAL(qd) , INTENT(IN)  :: Z 
      REAL(qd)               :: ZP(2),T,TZ
      

! For small arguments use taylor expansion
      IF ( ABS( Z ) < ACC3 ) THEN
! first derivative  w.r.t. Z
         D1_F_BETA_LRS_REG= Z**3/REAL(720,KIND=qd)-Z**5/REAL(2240,KIND=qd)
      ELSE
         T = TH( Z/REAL(2,KIND=qd) )
         TZ = T/Z
! write functions in terms of Tanh[ Z/2 ]
! first derivative  w.r.t. Z
         D1_F_BETA_LRS_REG=                       &
              - ( T/REAL(8,KIND=qd) )                   &
              + ( T**3/REAL(8,KIND=qd) )                &
              + ( REAL(2,KIND=qd)*TZ**2/Z  )              &
              - ( REAL(5,KIND=qd)/REAL(4,KIND=qd)*TZ/Z )          &
              + ( T*TZ**2 )                       &
              + ( REAL(1,KIND=qd)/( REAL(8,KIND=qd)*Z  ) )    &
              - ( T*TZ/REAL(8,KIND=qd) )    
      ENDIF 
      
   END FUNCTION D1_F_BETA_LRS_REG

!***********************************************************************
! DESCRIPTION:
!>calculates second dervatives of following function:
!> \f$
!>     \frac{1}{8}\left(1-\tanh^2\frac{z}{2}+\frac{2}{z}\tanh\frac{z}{2}\right)
!>     -\frac{1}{z^2}\tanh^2\frac{z}{2}
!> \f$
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION D2_F_BETA_LRS_REG( Z ) 
      REAL(qd)               :: D2_F_BETA_LRS_REG
      REAL(qd) , INTENT(IN)  :: Z 
      REAL(qd)               :: T,TZ
      

! For small arguments use taylor expansion
      IF ( ABS( Z ) < ACC3 ) THEN
! first derivative  w.r.t. Z
         D2_F_BETA_LRS_REG = Z**2/REAL(240,KIND=qd)-Z**4/REAL(448,KIND=qd)+REAL(7,KIND=qd)*Z**6/REAL(10800,KIND=qd)
      ELSE
         T = TH( Z/REAL(2,KIND=qd) )
         TZ = T/Z
! write functions in terms of Tanh[ Z/2 ]
! first derivative  w.r.t. Z
         D2_F_BETA_LRS_REG=                       &
              -  REAL(0.0625,KIND=qd)                     &
              + ( T**2/REAL(4,KIND=qd) )                &
              - (REAL(3,KIND=qd)/REAL(16,KIND=qd))*T**4           &
              - (REAL(6,KIND=qd)*(TZ**2/Z**2 ))           &
              + (REAL(9,KIND=qd)*TZ/(REAL(2,KIND=qd)*Z**2 ))      &
              - (REAL(4,KIND=qd)*TZ**3)                   &
              - (REAL(3,KIND=qd)/(REAL(4,KIND=qd)*Z**2 ))         &
              + (REAL(9,KIND=qd)/REAL(4,KIND=qd)*TZ**2)           &
              - (REAL(3,KIND=qd)/REAL(2,KIND=qd)*T**2*TZ**2)      &
              - (TZ/REAL(8,KIND=qd) )                     &
              + (TZ*T**2/REAL(8,KIND=qd) )      
      ENDIF 
      
   END FUNCTION D2_F_BETA_LRS_REG

!***********************************************************************
! DESCRIPTION:
!>calculates following function:
!> \f$
!>     \frac{1}{4}\left( \frac{ \tanh\frac{z}{2} }{z}
!>                 -\frac12( 1-\tanh^2\frac{z}{2} ) \right)
!> \f$
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION F_BETA_RRS( Z ) 
      REAL(qd)               :: F_BETA_RRS
      REAL(qd) , INTENT(IN)  :: Z 
      

! use a taylor expansion for small arguments
      IF ( ABS( Z ) < REAL(1E-8,KIND=qd) ) THEN
         F_BETA_RRS = Z**2/REAL(48,KIND=qd) - Z**4/REAL(240,KIND=qd) + REAL(17,KIND=qd)*Z**6/REAL(26880,KIND=qd)
      ELSE
         F_BETA_RRS = ( REAL(2,KIND=qd)*F_THZ_OVER_Z( Z ) &
                    - ( REAL(1,KIND=qd)-TH( Z/REAL(2,KIND=qd) )**2 ) )/REAL(8,KIND=qd)
      ENDIF
      
   END FUNCTION F_BETA_RRS

!***********************************************************************
! DESCRIPTION:
!>calculates first derivative of following function:
!> \f$
!>     \frac{1}{4}\left( \frac{ \tanh\frac{z}{2} }{z}
!>                 -\frac12\left(1-\tanh^2\frac{z}{2}\right) \right)
!> \f$
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION D1_F_BETA_RRS( Z ) 
      REAL(qd)              :: D1_F_BETA_RRS
      REAL(qd), INTENT(IN)  :: Z 
! local
      REAL(qd)              :: TZ, T

      

! use a taylor expansion for small arguments
      IF ( ABS( Z ) < REAL(1E-8,KIND=qd) ) THEN
         D1_F_BETA_RRS=Z/REAL(24,KIND=qd)-Z**3/REAL(60,KIND=qd)+REAL(17,KIND=qd)*Z**5/REAL(4480,KIND=qd)
      ELSE
         TZ = F_THZ_OVER_Z( Z )
         T = TH( Z/REAL(2,KIND=qd) ) 
!T/8 - T^3/8 - T/(4 z^2) + 1/(8 z) - T^2/(8 z)
         D1_F_BETA_RRS = (T-T**3)/REAL(8,KIND=qd)-(REAL(2,KIND=qd)*TZ/Z+T*TZ)/REAL(8,KIND=qd) &
                       + REAL(1,KIND=qd)/(REAL(8,KIND=qd)*Z)
      ENDIF
      
   END FUNCTION D1_F_BETA_RRS

!***********************************************************************
! DESCRIPTION:
!>calculates second derivative of following function:
!> \f$
!>     \frac{1}{4}\left( \frac{ \tanh\frac{z}{2} }{z}
!>                 -\frac12(1-\tanh^2\frac{z}{2}) \right)
!> \f$
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION D2_F_BETA_RRS( Z ) 
      REAL(qd)              :: D2_F_BETA_RRS
      REAL(qd), INTENT(IN)  :: Z 
! local
      REAL(qd)              :: TZ, T

      

! use a taylor expansion for small arguments
      IF ( ABS( Z ) < REAL(1E-8,KIND=qd) ) THEN
         D2_F_BETA_RRS = REAL(1,KIND=qd)/REAL(24,KIND=qd)-Z**2/REAL(20,KIND=qd) &
                       + REAL(17,KIND=qd)*Z**4/REAL(896,KIND=qd) &
                       - REAL(31,KIND=qd)*Z**6/REAL(6480,KIND=qd)
      ELSE
         TZ = F_THZ_OVER_Z( Z )
         T = TH( Z/REAL(2,KIND=qd) ) 
!1/16 - T^2/4 + (3 T^4)/16 + T/(2 z^3) - 1/(4 z^2) + T^2/(4 z^2) - T/( 8 z) + T^3/(8 z)
         D2_F_BETA_RRS = REAL(1,KIND=qd)/REAL(16,KIND=qd)       &
                       - T**2/REAL(4,KIND=qd)         &
                       + REAL(3,KIND=qd)*T**4/REAL(16,KIND=qd)  &
                       + TZ/(REAL(2,KIND=qd)*Z**2)      &
                       - REAL(1,KIND=qd)/(REAL(4,KIND=qd)*Z**2) &
                       + TZ**2/REAL(4,KIND=qd)          &
                       - TZ/REAL(8,KIND=qd)             &
                       + T**2*TZ/REAL(8,KIND=qd)      
      ENDIF
      
   END FUNCTION D2_F_BETA_RRS

   FUNCTION F_BETA_RRS_DER( Z ) 
      REAL(qd)               :: F_BETA_RRS_DER(2)
      REAL(qd) , INTENT(IN)  :: Z 
      REAL(qd)               :: T, TZ 
      

      IF ( ABS( Z ) < REAL(1E-8,KIND=qd) ) THEN
         F_BETA_RRS_DER(1)=Z/REAL(24,KIND=qd)-Z**3/REAL(60,KIND=qd)+REAL(17,KIND=qd)*Z**5/REAL(4480,KIND=qd)
         F_BETA_RRS_DER(2) = REAL(1,KIND=qd)/REAL(24,KIND=qd)-Z**2/REAL(20,KIND=qd) &
                       + REAL(17,KIND=qd)*Z**4/REAL(896,KIND=qd) &
                       - REAL(31,KIND=qd)*Z**6/REAL(6480,KIND=qd)
      ELSE
         TZ = F_THZ_OVER_Z( Z )
         T = TH( Z/REAL(2,KIND=qd) ) 
!T/8 - T^3/8 - T/(4 z^2) + 1/(8 z) - T^2/(8 z)
         F_BETA_RRS_DER(1) = (T-T**3)/REAL(8,KIND=qd)-(REAL(2,KIND=qd)*TZ/Z+T*TZ)/REAL(8,KIND=qd) &
                       + REAL(1,KIND=qd)/(REAL(8,KIND=qd)*Z)
         F_BETA_RRS_DER(2) = REAL(1,KIND=qd)/REAL(16,KIND=qd)       &
                       - T**2/REAL(4,KIND=qd)         &
                       + REAL(3,KIND=qd)*T**4/REAL(16,KIND=qd)  &
                       + TZ/(REAL(2,KIND=qd)*Z**2)      &
                       - REAL(1,KIND=qd)/(REAL(4,KIND=qd)*Z**2) &
                       + TZ**2/REAL(4,KIND=qd)          &
                       - TZ/REAL(8,KIND=qd)             &
                       + T**2*TZ/REAL(8,KIND=qd)      
      ENDIF

      
   END FUNCTION F_BETA_RRS_DER

!***********************************************************************
! DESCRIPTION:
!>calculates this function:
!> \f$
!>     \frac{1}{2} ( \coth\frac{z}{2} -\frac{2}{z} )
!> \f$
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION COTH_REG( Z ) 
      REAL(qd)               :: COTH_REG
      REAL(qd) , INTENT(IN)  :: Z 
      REAL(qd)               :: ARG 

      

      ARG = Z/REAL(2,KIND=qd)
! use a taylor expansion for small arguments
      IF ( ABS( Z ) < REAL(0.3,KIND=qd) ) THEN
         COTH_REG=ARG/REAL(6,KIND=qd)-ARG**3/REAL(90,KIND=qd)+ARG**5/REAL(945,KIND=qd)
      ELSE
         COTH_REG = (REAL(1,KIND=qd)/TH(ARG)-REAL(1,KIND=qd)/ARG)/REAL(2,KIND=qd)
      ENDIF
 
      
   END FUNCTION COTH_REG

!***********************************************************************
! DESCRIPTION:
!>calculates derivative of this function:
!> \f$
!>     \frac{1}{2} ( \coth\frac{z}{2} -\frac{2}{z} )
!> \f$
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION D1_COTH_REG( Z ) 
      REAL(qd)               :: D1_COTH_REG
      REAL(qd) , INTENT(IN)  :: Z 
      REAL(qd)               :: ARG, CT 

      

      ARG = Z/REAL(2,KIND=qd)
! use a taylor expansion for small arguments
      IF ( ABS( Z ) < REAL(0.3,KIND=qd) ) THEN
         D1_COTH_REG=REAL(1,KIND=qd)/REAL(12,KIND=qd)   &
                             -ARG**2/REAL(60,KIND=qd)    &
                             +ARG**4/REAL(378,KIND=qd)   &
                             -ARG**6/REAL(2700,KIND=qd) 
      ELSE
         CT=REAL(1,KIND=qd)/TH(ARG)
         D1_COTH_REG=( REAL(1,KIND=qd)/(ARG*ARG) -CT**2 + REAL(1,KIND=qd) )/REAL(4,KIND=qd)
      ENDIF
 
      
   END FUNCTION D1_COTH_REG

!***********************************************************************
!
! DESCRIPTION:
!>calculates second derivative of
!> \f$
!>    f''(z)= \frac14(\coth^3\frac{z}{2}-\coth\frac{z}{2}-
!>            \frac{8}{z^3})
!> \f$
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION D2_COTH_REG( Z ) 
      REAL(qd)               :: D2_COTH_REG
      REAL(qd) , INTENT(IN)  :: Z 
      REAL(qd)               :: ARG , CT

      

      ARG = Z/REAL(2,KIND=qd)
! use a taylor expansion for small arguments
      IF ( ABS( Z ) < REAL(0.3,KIND=qd) ) THEN
         D2_COTH_REG=-ARG/REAL(60,KIND=qd)      &
                     +ARG**3/REAL(189,KIND=qd)   &
                     -ARG**5/REAL(900,KIND=qd)   &
                     +ARG**7/REAL(20790,KIND=qd) 
      ELSE
         CT=REAL(1,KIND=qd)/TH(ARG)
         D2_COTH_REG=(-REAL(1,KIND=qd)/ARG**3 - CT + CT**3 )/REAL(4,KIND=qd) 
      ENDIF
 
      
   END FUNCTION D2_COTH_REG

!***********************************************************************
!
! DESCRIPTION:
!>calculates taylor expansion of cos for small arguments Z
!> \f$
!>    f(z)  = 1 - z^2/2
!> \f$
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION COS_TAYLOR( Z ) 
      REAL(qd)               :: COS_TAYLOR
      REAL(qd) , INTENT(IN)  :: Z 

      

      COS_TAYLOR = REAL(1,KIND=qd) - Z**2/REAL(2,KIND=qd) 

      
   END FUNCTION COS_TAYLOR

!***********************************************************************
!
! DESCRIPTION:
!>calculates taylor expansion of sin for small arguments Z
!> \f$
!>    f(z)  = z - z^3/6
!> \f$
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION SIN_TAYLOR( Z ) 
      REAL(qd)               :: SIN_TAYLOR
      REAL(qd) , INTENT(IN)  :: Z 

      

      SIN_TAYLOR = Z - Z**3/REAL(6,KIND=qd) 

      
   END FUNCTION SIN_TAYLOR

!****************************************************************************
!
! DESCRIPTION:
!>returns sign of a quad-precision real-valued variable
!
!> @param[in]  X  (quad precision real)
!
!****************************************************************************

   PURE FUNCTION SIGN_QD( X ) 
   INTEGER               :: SIGN_QD 
   REAL(qd) , INTENT(IN)  :: X 
! quad library has no build in functionality for sign function
!#ifdef qd_emulate
      IF ( X > 0 ) THEN
         SIGN_QD = 1
      ELSE
         SIGN_QD =-1
      ENDIF
!#else
!SIGN_QD = SIGN( 1._q, REAL( X, q) )
!#endif
   END FUNCTION SIGN_QD

!****************************************************************************
! DESCRIPTION:
!>calculates the sum of all components of a vector A
!
!> @param[in]  A  (argument)
!****************************************************************************

   FUNCTION SUM_MPR( A  ) 
      USE prec
# 877

      REAL(qd) :: A(:)
      REAL(qd) :: SUM_MPR, S
      INTEGER         :: I, N
   
      
   
      S=0
      N = SIZE( A ) 
      DO I = 1, N
         S = S + A(I)
      ENDDO
      SUM_MPR = S 
   
      
   END FUNCTION SUM_MPR 

END MODULE minimax_functions1D
