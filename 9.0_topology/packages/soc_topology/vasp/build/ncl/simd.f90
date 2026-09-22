# 1 "simd.F"
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





# 2 "simd.F" 2 

MODULE SIMD
      TYPE, PUBLIC :: SIMD_REAL64
      REAL*8 :: X(0 : 1 * 1 - 1)
      END TYPE SIMD_REAL64

      TYPE, PUBLIC :: SIMD_MASK_REAL64
      LOGICAL :: X(0 : 1 * 1 - 1)
      END TYPE SIMD_MASK_REAL64

      TYPE, PUBLIC :: SIMD_INT64
      INTEGER*8 :: X(0 : 1 * 1 - 1)
      END TYPE SIMD_INT64

# 45

END MODULE SIMD
