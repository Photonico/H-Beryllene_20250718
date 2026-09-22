# 1 "minimax_functions2D.F"
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


# 2 "minimax_functions2D.F" 2 

!***********************************************************************
!
! MODULE: minimax_functions
!
!> @author
!> Merzuk Kaltak, VASP Software GmbH
!
! DESCRIPTION:
!> contains fitting and error functions used in minimax module, mostly
!> basis functions that span the isometric subspaces.
!>
!***********************************************************************

MODULE minimax_functions2D
   USE prec
   USE minimax_functions1D    
# 21


   CONTAINS 

!***********************************************************************
! DESCRIPTION:
!>calculates basis function:
!> \f$
!>     \phi(x,\lambda) = e^{-z \lambda}
!> \f$
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  Z  (argument of error function)
!
!***********************************************************************

   FUNCTION EXPF( Z, L )
      REAL(qd)               :: EXPF
      REAL(qd) , INTENT(IN)  :: Z ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function
      
      
! IF (-Z*L > 100.0_q) THEN; EXPF = 1.0E10_q; ; RETURN; END IF ! Breaks SiC8_ACFDTR (gnu), exp overflow possible.
      EXPF = EXP( -Z*L )
      
   END FUNCTION EXPF

!***********************************************************************
! DESCRIPTION:
!>calculates derivative of basis function w.r.t. first argument
!> \f$
!>     \phi(z,\lambda) = e^{-z \lambda}
!> \f$
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  Z  (argument of error function)
!
!***********************************************************************

   FUNCTION D_EXPF_DZ( Z, L )
      REAL(qd)               :: D_EXPF_DZ
      REAL(qd) , INTENT(IN)  :: Z ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function
      
      
      D_EXPF_DZ = -L*EXPF( Z,L )

      
   END FUNCTION D_EXPF_DZ

!***********************************************************************
! DESCRIPTION:
!>calculates first derivative of basis function w.r.t. second argument
!> \f$
!>     \phi(z,\lambda) = e^{-z \lambda}
!> \f$
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  Z  (argument of error function)
!
!***********************************************************************

   FUNCTION D_EXPF_DL( Z, L )
      REAL(qd)               :: D_EXPF_DL
      REAL(qd) , INTENT(IN)  :: Z ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function
      
      
      D_EXPF_DL = -Z*EXPF( Z,L )

      
   END FUNCTION D_EXPF_DL

!***********************************************************************
! DESCRIPTION:
!>calculates second derivative of basis function w.r.t. first argument
!> \f$
!>     \phi(z,\lambda) = e^{-z \lambda}
!> \f$
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  Z  (argument of error function)
!
!***********************************************************************

   FUNCTION D2_EXPF_DZ2( Z, L )
      REAL(qd)               :: D2_EXPF_DZ2
      REAL(qd) , INTENT(IN)  :: Z ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function
      
      
      D2_EXPF_DZ2 = -L*D_EXPF_DZ( Z,L )

      
   END FUNCTION D2_EXPF_DZ2

!***********************************************************************
! DESCRIPTION:
!>calculates basis function:
!> \f$
!>     \phi(z,\lambda) = \frac{ z }{ z^2 + \lambda^2 }
!> \f$
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  Z  (argument of error function)
!
!***********************************************************************

   FUNCTION UFREQ( Z, L )
      REAL(qd)              :: UFREQ
      REAL(qd), INTENT(IN)  :: Z ! sampling value
      REAL(qd), INTENT(IN)  :: L ! coefficient defining function
      
      
      UFREQ = Z/( Z*Z+L*L )

      
   END FUNCTION UFREQ 

!***********************************************************************
! DESCRIPTION:
!>calculates derivative of basis function w.r.t. first argument
!> \f$
!>     \phi(z,\lambda) = \frac{ z }{ z^2 + \lambda^2 }
!> \f$
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  Z  (argument of error function)
!
!***********************************************************************

   FUNCTION D_UFREQ_DZ( Z, L )
      REAL(qd)               :: D_UFREQ_DZ
      REAL(qd) , INTENT(IN)  :: Z ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function
      
      
      D_UFREQ_DZ = VFREQ(Z,L)**2 - UFREQ(Z,L)**2

      
   END FUNCTION D_UFREQ_DZ

!***********************************************************************
! DESCRIPTION:
!>calculates derivative of basis function w.r.t. second argument
!> \f$
!>     \phi(z,\lambda) = \frac{ z }{ z^2 + \lambda^2 }
!> \f$
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  Z  (argument of error function)
!
!***********************************************************************

   FUNCTION D_UFREQ_DL( Z, L )
      REAL(qd)               :: D_UFREQ_DL
      REAL(qd) , INTENT(IN)  :: Z ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function
      
      
      D_UFREQ_DL = -REAL(2,KIND=qd)*VFREQ(Z,L)*UFREQ(Z,L)

      
   END FUNCTION D_UFREQ_DL

!***********************************************************************
! DESCRIPTION:
!>calculates the second two derivatives of this function:
!> \f$
!>     \frac{ x }{ x^2 + \lambda^2 }
!> \f$
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION D2_UFREQ_DZ2( X, L )
      REAL(qd)               :: D2_UFREQ_DZ2
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function
      
      

      D2_UFREQ_DZ2 = REAL(2,KIND=qd)*UFREQ(X,L)*&
         (UFREQ(X,L)**2 -REAL(3,KIND=qd)*VFREQ(X,L)**2 )

      
   END FUNCTION D2_UFREQ_DZ2

!***********************************************************************
! DESCRIPTION:
!>calculates basis function:
!> \f$
!>     \phi(x,\lambda)=\frac{ \lambda }{ x^2 + \lambda^2 }
!> \f$
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION VFREQ( X, L )
      REAL(qd)               :: VFREQ
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function
      
      

      VFREQ = L/( X*X+L*L )

      
   END FUNCTION VFREQ 

!***********************************************************************
! DESCRIPTION:
!>calculates the first derivative of this function (w.r.t. to x)
!> \f$
!>     \frac{ \lambda }{ x^2 + \lambda^2 }
!> \f$
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION D_VFREQ_DZ( X, L )
      REAL(qd)               :: D_VFREQ_DZ
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function

      
      
      D_VFREQ_DZ = -REAL(2,KIND=qd)*UFREQ(X,L)*VFREQ(X,L)
      
      

   END FUNCTION D_VFREQ_DZ

!***********************************************************************
! DESCRIPTION:
!>calculates the first derivative of this function (w.r.t. to L)
!> \f$
!>     \frac{ \lambda }{ x^2 + \lambda^2 }
!> \f$
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION D_VFREQ_DL( Z, L )
      REAL(qd)               :: D_VFREQ_DL
      REAL(qd) , INTENT(IN)  :: Z ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function

      
      
      D_VFREQ_DL = UFREQ(Z,L)**2 - VFREQ(Z,L)**2
      
      

   END FUNCTION D_VFREQ_DL

!***********************************************************************
! DESCRIPTION:
!>calculates the secondderivative of this function (w.r.t. to x)
!> \f$
!>     \frac{ \lambda }{ x^2 + \lambda^2 }
!> \f$
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION D2_VFREQ_DZ2( X, L )
      REAL(qd)               :: D2_VFREQ_DZ2
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function

      
      
      D2_VFREQ_DZ2 = -REAL(2,KIND=qd)*( D_UFREQ_DZ(X,L)*VFREQ(X,L) &
                              + D_VFREQ_DZ(X,L)*UFREQ(X,L) )
      
      

   END FUNCTION D2_VFREQ_DZ2

!***********************************************************************
! DESCRIPTION:
!>calculates basis function for time domain :
!> \f$
!>     \frac{ \cosh(x/2(1-2*L)) }{ \cosh(x/2) }
!> \f$
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION UTIME( X, L )
      REAL(qd)               :: UTIME
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function
      
      
      UTIME = COSHCOSH( X/REAL(2,KIND=qd), REAL(1,KIND=qd)-REAL(2,KIND=qd)*L )

      
   END FUNCTION UTIME

!***********************************************************************
! DESCRIPTION:
!>calculates odd function for time domain :
!> \f$
!>     \frac{ \sinh(x/2(1-2*L)) }{ \cosh(x/2) }
!> \f$
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION VTIME( X, L )
      REAL(qd)               :: VTIME
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function
      
      
      VTIME = SINHCOSH( X/REAL(2,KIND=qd), REAL(1,KIND=qd)-REAL(2,KIND=qd)*L )

      
   END FUNCTION VTIME

!***********************************************************************
! DESCRIPTION:
!>calculates first derivative of follwing function (w.r.t. to x )
!> \f$
!>     \frac{ \cosh(z/2(1-2*L)) }{ \cosh(z/2) }
!> \f$
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  Z  (argument of error function)
!
!***********************************************************************

   FUNCTION D_UTIME_DZ( Z, L )
      REAL(qd)               :: D_UTIME_DZ
      REAL(qd) , INTENT(IN)  :: Z ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function
      
      
      D_UTIME_DZ = DCOSHCOSHX( Z/REAL(2,KIND=qd), REAL(1,KIND=qd)-REAL(2,KIND=qd)*L )/REAL(2,KIND=qd)

      
   END FUNCTION D_UTIME_DZ

!***********************************************************************
! DESCRIPTION:
!>calculates first derivative of follwing function (w.r.t. to L )
!> \f$
!>     \frac{ \cosh(z/2(1-2*L)) }{ \cosh(z/2) }
!> \f$
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  Z  (argument of error function)
!
!***********************************************************************

   FUNCTION D_UTIME_DL( Z, L )
      REAL(qd)               :: D_UTIME_DL
      REAL(qd) , INTENT(IN)  :: Z ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function
      
      
      D_UTIME_DL = -REAL(2,KIND=qd)*DCOSHCOSHY( Z/REAL(2,KIND=qd), REAL(1,KIND=qd)-REAL(2,KIND=qd)*L )

      
   END FUNCTION D_UTIME_DL

!***********************************************************************
! DESCRIPTION:
!>calculates second derivative of follwing function (w.r.t. to x )
!> \f$
!>     \frac{ \cosh\left(\frac{z}{2}(1-2L)\right) }{ \cosh\frac{z}{2} }
!> \f$
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  Z  (argument of error function)
!
!***********************************************************************

   FUNCTION D2_UTIME_DZ2( Z, L )
      REAL(qd)               :: D2_UTIME_DZ2
      REAL(qd) , INTENT(IN)  :: Z ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function
      
      
      D2_UTIME_DZ2 = DCOSHCOSHX2( Z/REAL(2,KIND=qd), REAL(1,KIND=qd)-REAL(2,KIND=qd)*L )/REAL(4,KIND=qd)

      
   END FUNCTION D2_UTIME_DZ2

!***********************************************************************
! DESCRIPTION:
!>calculates basis function for bosonic Matsubara freuencies:
!> \f$
!>     \frac{ x \tanh(x/2) }{ x^2+L^2 }
!> \f$
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION UTANH( X, L )
      REAL(qd)               :: UTANH
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function
      
      
      IF ( L == 0 ) THEN
         UTANH = F_THZ_OVER_Z( X ) 
      ELSE
! taylor expansion for small arguments
         IF ( ABS( X ) < REAL(1E-8,KIND=qd)  ) THEN
            UTANH = REAL(0.5,KIND=qd)*(X/L)**2                 & 
                          - (REAL(12,KIND=qd)+L**2)/REAL(24,KIND=qd)*(X/L)**4  &
                          + (REAL(120,KIND=qd)+REAL(10,KIND=qd)*L**2+L**4)/REAL(240,KIND=qd)*(X/L)**4
         ELSE
            UTANH = TH( X/REAL(2,KIND=qd) ) * UFREQ( X, L ) 
         ENDIF
      ENDIF
      
   END FUNCTION UTANH 

!***********************************************************************
! DESCRIPTION:
!>calculates first derivative of basis function w.r.t. x :
!> \f$
!>     \frac{ x \tanh(x/2) }{ x^2+L^2 }
!> \f$
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION D_UTANH_DZ( X, L )
      REAL(qd)               :: D_UTANH_DZ
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function
! local
      REAL(qd)               :: T,C,CP
      
    
      T = TH( X/2 )
      C  = UFREQ(X,L)
      CP = D_UFREQ_DZ(X,L)

      D_UTANH_DZ = (REAL(1,KIND=qd)-T**2)*C/REAL(2,KIND=qd) + T*CP

      
   END FUNCTION D_UTANH_DZ 

!***********************************************************************
! DESCRIPTION:
!>calculates first derivative of basis function w.r.t. L :
!> \f$
!>     \frac{ x \tanh(x/2) }{ x^2+L^2 }
!> \f$
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION D_UTANH_DL( X, L )
      REAL(qd)               :: D_UTANH_DL
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function
! local
      
    
      D_UTANH_DL = D_UFREQ_DL(X,L)*TH(X/2)

      
   END FUNCTION D_UTANH_DL 

!***********************************************************************
! DESCRIPTION:
!>calculates second derivative of basis function w.r.t. x :
!> \f$
!>     \frac{ x \tanh(x/2) }{ x^2+L^2 }
!> \f$
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION D2_UTANH_DZ2( X, L )
      REAL(qd)               :: D2_UTANH_DZ2
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function
! local
      REAL(qd)               :: T,C,CP(2)
      
    
      T = TH( X/2 )
      C  = UFREQ(X,L)
      CP(1) = D_UFREQ_DZ(X,L)
      CP(2) = D2_UFREQ_DZ2(X,L)

      D2_UTANH_DZ2 = (REAL(1,KIND=qd)-T**2)*CP(1) + T*CP(2) &
                           - T*(REAL(1,KIND=qd)-T**2)*C/REAL(2,KIND=qd)

      
   END FUNCTION D2_UTANH_DZ2 

!***********************************************************************
! DESCRIPTION:
!>calculates basis function for bosonic Matsubara freuencies:
!> \f$
!>     \frac{ L \tanh x/2}{ x^2+L^2 }
!> \f$
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION VTANH( X, L )
      REAL(qd)               :: VTANH
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function
      
      VTANH = VFREQ( X, L )*TH(X/2) 

   END FUNCTION VTANH 

!***********************************************************************
! DESCRIPTION:
!>calculates frist derivative of this basis function w.r.t. to X
!> \f$
!>     \frac{ L \tanh x/2}{ x^2+L^2 }
!> \f$
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION D_VTANH_DZ( X, L )
      REAL(qd)               :: D_VTANH_DZ
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function
      REAL(qd)               :: T
      
      T = TH(X/2)

      D_VTANH_DZ = D_VFREQ_DZ( X, L )*T + VFREQ(X,L)*(REAL(1,KIND=qd)-T**2)/REAL(2,KIND=qd)

      
   END FUNCTION D_VTANH_DZ

!***********************************************************************
! DESCRIPTION:
!>calculates frist derivative of this basis function w.r.t. to L
!> \f$
!>     \frac{ L \tanh x/2}{ x^2+L^2 }
!> \f$
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION D_VTANH_DL( X, L )
      REAL(qd)               :: D_VTANH_DL
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function
      

      D_VTANH_DL = D_VFREQ_DL( X, L )*TH(X/2)

      
   END FUNCTION D_VTANH_DL

!***********************************************************************
! DESCRIPTION:
!>calculates frist derivative of this basis function w.r.t. to x
!> \f$
!>     \frac{ L \tanh x/2}{ x^2+L^2 }
!> \f$
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION D2_VTANH_DZ2( X, L )
      REAL(qd)               :: D2_VTANH_DZ2
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function
      REAL(qd)               :: T,DT
      
      T = TH(X/2)
      DT = (REAL(1,KIND=qd)-T**2)

      D2_VTANH_DZ2 = T*D2_VFREQ_DZ2(X,L) &
                   + DT*( D_VFREQ_DZ(X,L)-T*VFREQ(X,L)/REAL(2,KIND=qd) )

      
   END FUNCTION D2_VTANH_DZ2

!****************************************************************************
!
! DESCRIPTION:
!>calculates scalar product of two vectors in quad-precision:
!
!> @param[in]  A  (first vector)
!> @param[in]  B  (second vector)
!
!****************************************************************************

FUNCTION DOT_PRODUCT_MPR( A , B ) 
   USE prec
# 653

   REAL(qd) :: A(:), B(:)
   REAL(qd) :: DOT_PRODUCT_MPR, P
   INTEGER :: I, N

   

   P=0
   N = SIZE( A ) 
   IF ( N == SIZE( B ) ) THEN 
# 668

      DOT_PRODUCT_MPR = DOT_PRODUCT( a, b ) 

   ELSE
      WRITE( *, * )&
      ' INTERNAL ERROR IN VASP: DOT_PRODUCT_MPR REPORTS INCONSISTENT SIZES',&
      N, SIZE( B ) 
   ENDIF 

   

ENDFUNCTION DOT_PRODUCT_MPR 

!****************************************************************************
! DESCRIPTION:
!>calculates following function
!> \f$
!>   \frac{\cosh( x y)}{\cosh(x)}
!> \f$
!
!> @param[in]  X  (argument)
!> @param[in]  Y  (argument)
!****************************************************************************

   FUNCTION COSHCOSH( X, Y ) 
      USE prec
# 696

      REAL(qd) :: X, Y
      REAL(qd) :: COSHCOSH
! local
      REAL(qd) :: AX, AXY, D
   
      
    
! This function is a beast!
! For finite temperature quadratures a naive evaluation will
! certainly produce an over or underflow or even a NaN
      AX = ABS( X )
      AXY = ABS( X*Y )

! if both arguments are larger than 1000 :
      IF ( AX > REAL(MAXE,KIND=qd)/REAL(2,KIND=qd) .AND. AXY > REAL(MAXE,KIND=qd)/REAL(2,KIND=qd) ) THEN
! here the difference of arguments is important
         D = X*Y - X 
         COSHCOSH = EXP( SIGN_QD(D)*MIN( ABS(D), REAL(MAXE,KIND=qd)/REAL(2,KIND=qd) ) ) 
      ELSE
! in case denominator is larger than nummerator,
! prevent overflow in nummerator
         IF ( AX - AXY > REAL(MAXE,KIND=qd) ) THEN
            COSHCOSH = EXP( -MIN( ABS( AX - AXY ) , REAL(MAXE,KIND=qd) ) )
         ELSE IF ( AXY - AX > REAL(MAXE,KIND=qd) ) THEN
            COSHCOSH = EXP( MIN( ABS( AX - AXY ) , REAL(MAXE,KIND=qd) ) )
         ELSE
            COSHCOSH = COSH(AXY)/COSH(AX)
         ENDIF
# 732

      ENDIF

      
   END FUNCTION COSHCOSH

!****************************************************************************
! DESCRIPTION:
!>calculates following function
!> \f$
!>   \frac{\sinh( x y)}{\cosh(x)}
!> \f$
!
!> @param[in]  X  (argument)
!> @param[in]  Y  (argument)
!****************************************************************************

   FUNCTION SINHCOSH( X, Y ) 
      USE prec
# 753

      REAL(qd) :: X, Y
      REAL(qd) :: SINHCOSH
! local
      REAL(qd) :: AX, AXY, D
   
      
    
! This function is probably also a beast!
      AX = ABS( X )
      AXY = ABS( X*Y )

! if both arguments are larger than 1000 :
      IF ( AX > REAL(MAXE,KIND=qd)/REAL(2,KIND=qd) .AND. AXY > REAL(MAXE,KIND=qd)/REAL(2,KIND=qd) ) THEN
! here the difference of arguments is important
         D = X*Y - X 
         SINHCOSH = SIGN_QD(X*Y)*EXP( SIGN_QD(D)*MIN( ABS(D), REAL(MAXE,KIND=qd)/REAL(2,KIND=qd) ) ) 
      ELSE
! in case denominator is larger than nummerator,
! prevent overflow in nummerator
         IF ( AX - AXY > REAL(MAXE,KIND=qd) ) THEN
            SINHCOSH = SIGN_QD(X*Y)*EXP( -MIN( ABS( AX - AXY ) , REAL(MAXE,KIND=qd) ) )
         ELSE IF ( AXY - AX > REAL(MAXE,KIND=qd) ) THEN
            SINHCOSH = SIGN_QD(X*Y)*EXP( MIN( ABS( AX - AXY ) , REAL(MAXE,KIND=qd) ) )
         ELSE
            SINHCOSH = SINH(AXY)/COSH(AX)
         ENDIF
# 788

      ENDIF
   
      
   END FUNCTION SINHCOSH

!****************************************************************************
! DESCRIPTION:
!>calculates following derivative
!> \f$
!>   \frac{\partial}{\partial x} \frac{\cosh( x y)}{\cosh(x)}
!>  =y\frac{\sinh( x y )}{\cosh( x )} - \frac{\cosh( x y)}{\cosh(x)} tanh(x)
!> \f$
!
!> @param[in]  X  (argument)
!> @param[in]  Y  (argument)
!****************************************************************************

   FUNCTION DCOSHCOSHX( X, Y ) 
      USE prec
# 810

      REAL(qd) :: X, Y 
      REAL(qd) :: DCOSHCOSHX
      
   
      DCOSHCOSHX = Y * SINHCOSH( X, Y ) - COSHCOSH( X, Y )*TH( X )
   
      
   END FUNCTION DCOSHCOSHX

!****************************************************************************
! DESCRIPTION:
!>calculates following derivative
!> \f$
!>   \frac{\partial^2}{\partial x^2} \frac{\cosh( x y)}{\cosh(x)}
!> \f$
!
!> @param[in]  X  (argument)
!> @param[in]  Y  (argument)
!****************************************************************************

   FUNCTION DCOSHCOSHX2( X, Y ) 
      USE prec
# 835

      REAL(qd) :: X, Y 
      REAL(qd) :: DCOSHCOSHX2
   
      
      DCOSHCOSHX2 = COSHCOSH( X, Y )*( 2*TH( X )**2 + Y**2 - REAL(1,KIND=qd) ) &
                  - 2*Y*SINHCOSH( X, Y )*TH( X )
   
      
   END FUNCTION DCOSHCOSHX2

!****************************************************************************
! DESCRIPTION:
!>calculates following function
!> \f$
!>   \frac{\cosh( x y)^2}{\cosh(x)^2}
!> \f$
!
!> @param[in]  X  (argument)
!> @param[in]  Y  (argument)
!****************************************************************************

   FUNCTION COSH2COSH( X, Y ) 
      REAL(qd) :: X, Y
      REAL(qd) :: COSH2COSH
   
      COSH2COSH = COSHCOSH(X,Y)**2

   END FUNCTION COSH2COSH

!****************************************************************************
! DESCRIPTION:
!>calculates following derivative
!> \f$
!>   \frac{\partial}{\partial x} \frac{\cosh( x y)^2}{\cosh(x)^2}
!> \f$
!
!> @param[in]  X  (argument)
!> @param[in]  Y  (argument)
!****************************************************************************

   FUNCTION DCOSH2COSHX( X, Y ) 
      USE prec
# 880

      REAL(qd) :: X, Y 
      REAL(qd) :: DCOSH2COSHX
   
      DCOSH2COSHX = REAL(2,KIND=qd)*COSH2COSH(X,Y)*( Y * TH( X*Y ) - TH( X ) ) 
   
   END FUNCTION DCOSH2COSHX

!****************************************************************************
! DESCRIPTION:
!>calculates following derivative
!> \f$
!>   \frac{\partial^2}{\partial x^2} \frac{\cosh( x y)^2}{\cosh(x)^2}
!> \f$
!
!> @param[in]  X  (argument)
!> @param[in]  Y  (argument)
!****************************************************************************

   FUNCTION D2COSH2COSHX2( X, Y ) 
      USE prec
# 903

      REAL(qd) :: X, Y 
      REAL(qd) :: D2COSH2COSHX2
   
      D2COSH2COSHX2 = REAL(2,KIND=qd)*COSH2COSH(X,Y)*( &
                      Y**2 * TH( X*Y )**2      &
                    - REAL(4,KIND=qd)*Y*TH( X )*TH( X*Y )  &
                    + REAL(3,KIND=qd)*TH( X )**2        &
                    + Y**2-REAL(1,KIND=qd)  ) 
   
   END FUNCTION D2COSH2COSHX2

!****************************************************************************
! DESCRIPTION:
!>calculates following derivative
!> \f$
!>   \frac{\partial}{\partial y} \frac{\cosh( x y)^2}{\cosh(x)^2}
!> \f$
!
!> @param[in]  X  (argument)
!> @param[in]  Y  (argument)
!****************************************************************************

   FUNCTION DCOSH2COSHY( X, Y ) 
      USE prec
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
# 932

      REAL(qd) :: X, Y 
      REAL(qd) :: DCOSH2COSHY
   
! need to carefully evaluate this function
      DCOSH2COSHY = REAL(2,KIND=qd) * X * COSH2COSH( X, Y ) * TH( X*Y )  
   
   END FUNCTION DCOSH2COSHY

!****************************************************************************
! DESCRIPTION:
!>calculates following derivative
!> \f$
!>   \frac{\partial}{\partial y} \frac{\cosh( x y)}{\cosh(x)}
!> \f$
!
!> @param[in]  X  (argument)
!> @param[in]  Y  (argument)
!****************************************************************************

   FUNCTION DCOSHCOSHY( X, Y ) 
      USE prec
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
# 958

      REAL(qd) :: X, Y 
      REAL(qd) :: DCOSHCOSHY
   
!
      DCOSHCOSHY = X * SINHCOSH( X, Y ) 
   
   END FUNCTION DCOSHCOSHY


!=========================================================================
! Basis definitions, squared functions basically
!=========================================================================
   FUNCTION EXPF2( X, L )
      REAL(qd)               :: EXPF2
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function

      
      
      EXPF2 = EXPF(2*X,L)
      
      

   END FUNCTION EXPF2

   FUNCTION D_EXPF2_DZ( X, L )
      REAL(qd)               :: D_EXPF2_DZ
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function

      
      
      D_EXPF2_DZ = REAL(2,KIND=qd)*D_EXPF_DZ(REAL(2,KIND=qd)*X,L)
      
      

   END FUNCTION D_EXPF2_DZ

   FUNCTION D2_EXPF2_DZ2( X, L )
      REAL(qd)               :: D2_EXPF2_DZ2
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function

      
      
      D2_EXPF2_DZ2 = REAL(4,KIND=qd)*D2_EXPF_DZ2(REAL(2,KIND=qd)*X,L)
      
      

   END FUNCTION D2_EXPF2_DZ2

   FUNCTION D_EXPF2_DL( X, L )
      REAL(qd)               :: D_EXPF2_DL
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function

      
      
      D_EXPF2_DL = D_EXPF_DL(REAL(2,KIND=qd)*X,L)
      
      

   END FUNCTION D_EXPF2_DL

   FUNCTION UFREQ2( X, L )
      REAL(qd)               :: UFREQ2
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function

      
      
      UFREQ2 = UFREQ(X,L)**2
      
      

   END FUNCTION UFREQ2

   FUNCTION D_UFREQ2_DZ( X, L )
      REAL(qd)               :: D_UFREQ2_DZ
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function

      
      
      D_UFREQ2_DZ = REAL(2,KIND=qd)*UFREQ(X,L)*D_UFREQ_DZ(X,L)
      
      

   END FUNCTION D_UFREQ2_DZ

   FUNCTION D_UFREQ2_DL( X, L )
      REAL(qd)               :: D_UFREQ2_DL
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function

      
      
      D_UFREQ2_DL = REAL(2,KIND=qd)*UFREQ(X,L)*D_UFREQ_DL(X,L)
      
      

   END FUNCTION D_UFREQ2_DL

   FUNCTION D2_UFREQ2_DZ2( X, L )
      REAL(qd)               :: D2_UFREQ2_DZ2
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function

      
      
      D2_UFREQ2_DZ2 = REAL(2,KIND=qd)*( D_UFREQ_DZ(X,L)**2 &
                                   + UFREQ(X,L)*D2_UFREQ_DZ2(X,L) )
      
      

   END FUNCTION D2_UFREQ2_DZ2

   FUNCTION VFREQ2( X, L )
      REAL(qd)               :: VFREQ2
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function

      
      
      VFREQ2 = VFREQ(X,L)**2
      
      

   END FUNCTION VFREQ2

   FUNCTION D_VFREQ2_DZ( X, L )
      REAL(qd)               :: D_VFREQ2_DZ
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function

      
      
      D_VFREQ2_DZ = REAL(2,KIND=qd)*VFREQ(X,L)*D_VFREQ_DZ(X,L)
      
      

   END FUNCTION D_VFREQ2_DZ

   FUNCTION D_VFREQ2_DL( X, L )
      REAL(qd)               :: D_VFREQ2_DL
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function

      
      
      D_VFREQ2_DL = REAL(2,KIND=qd)*VFREQ(X,L)*D_VFREQ_DL(X,L)
      
      

   END FUNCTION D_VFREQ2_DL

   FUNCTION D2_VFREQ2_DZ2( X, L )
      REAL(qd)               :: D2_VFREQ2_DZ2
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function

      
      
      D2_VFREQ2_DZ2 = REAL(2,KIND=qd)*( D_VFREQ_DZ(X,L)**2 &
                                   + VFREQ(X,L)*D2_VFREQ_DZ2(X,L) )
      
      

   END FUNCTION D2_VFREQ2_DZ2

! T>0 basis functions
   FUNCTION UTIME2( X, L )
      REAL(qd)               :: UTIME2
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function

      
      
      UTIME2 = UTIME(X,L)**2
      
      

   END FUNCTION UTIME2

   FUNCTION D_UTIME2_DZ( X, L )
      REAL(qd)               :: D_UTIME2_DZ
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function

      
      
      D_UTIME2_DZ = REAL(2,KIND=qd)*UTIME(X,L)*D_UTIME_DZ(X,L)
      
      

   END FUNCTION D_UTIME2_DZ

   FUNCTION D2_UTIME2_DZ2( X, L )
      REAL(qd)               :: D2_UTIME2_DZ2
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function

      
      
      D2_UTIME2_DZ2 = REAL(2,KIND=qd)*( D_UTIME_DZ(X,L)**2 &
                                   + UTIME(X,L)*D2_UTIME_DZ2(X,L) )
      
      

   END FUNCTION D2_UTIME2_DZ2

   FUNCTION D_UTIME2_DL( X, L )
      REAL(qd)               :: D_UTIME2_DL
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function

      
      
      D_UTIME2_DL = REAL(2,KIND=qd)*UTIME(X,L)*D_UTIME_DL(X,L)
      
      

   END FUNCTION D_UTIME2_DL

! T>0 basis functions
   FUNCTION UTANH2( X, L )
      REAL(qd)               :: UTANH2
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function

      
      
      UTANH2 = UTANH(X,L)**2
      
      

   END FUNCTION UTANH2

   FUNCTION D_UTANH2_DZ( X, L )
      REAL(qd)               :: D_UTANH2_DZ
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function

      
      
      D_UTANH2_DZ = REAL(2,KIND=qd)*UTANH(X,L)*D_UTANH_DZ(X,L)
      
      

   END FUNCTION D_UTANH2_DZ

   FUNCTION D2_UTANH2_DZ2( X, L )
      REAL(qd)               :: D2_UTANH2_DZ2
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function

      
      
      D2_UTANH2_DZ2 = REAL(2,KIND=qd)*( D_UTANH_DZ(X,L)**2 &
                                   + UTANH(X,L)*D2_UTANH_DZ2(X,L) )
      
      

   END FUNCTION D2_UTANH2_DZ2

   FUNCTION D_UTANH2_DL( X, L )
      REAL(qd)               :: D_UTANH2_DL
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function

      
      
      D_UTANH2_DL = REAL(2,KIND=qd)*UTANH(X,L)*D_UTANH_DL(X,L)
      
      

   END FUNCTION D_UTANH2_DL

! T>0 basis functions
   FUNCTION VTANH2( X, L )
      REAL(qd)               :: VTANH2
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function

      
      
      VTANH2 = VTANH(X,L)**2
      
      

   END FUNCTION VTANH2

   FUNCTION D_VTANH2_DZ( X, L )
      REAL(qd)               :: D_VTANH2_DZ
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function

      
      
      D_VTANH2_DZ = REAL(2,KIND=qd)*VTANH(X,L)*D_VTANH_DZ(X,L)
      
      

   END FUNCTION D_VTANH2_DZ

   FUNCTION D2_VTANH2_DZ2( X, L )
      REAL(qd)               :: D2_VTANH2_DZ2
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function

      
      
      D2_VTANH2_DZ2 = REAL(2,KIND=qd)*( D_VTANH_DZ(X,L)**2 &
                                   + VTANH(X,L)*D2_VTANH_DZ2(X,L) )
      
      

   END FUNCTION D2_VTANH2_DZ2

   FUNCTION D_VTANH2_DL( X, L )
      REAL(qd)               :: D_VTANH2_DL
      REAL(qd) , INTENT(IN)  :: X ! sampling value
      REAL(qd) , INTENT(IN)  :: L ! coefficient defining function

      
      
      D_VTANH2_DL = REAL(2,KIND=qd)*VTANH(X,L)*D_VTANH_DL(X,L)
      
      

   END FUNCTION D_VTANH2_DL


END MODULE minimax_functions2D
