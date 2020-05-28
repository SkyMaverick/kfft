#pragma once

#ifdef KFFT_USE_SYSMATH
    #include <math.h>
#endif
#define KFFT_CONST_PI 3.141592653589793238462643383279502884197169399375105820974944

#if defined(KFFT_CC_GCC)
    #if __GNUC__ < 8

    /* This function emulated because GCC < 8.0.0 don't define it's intrinsics
          See bug https://gcc.gnu.org/bugzilla/show_bug.cgi?id=80582 */

    // clang-format off
        #define _mm256_set_m128d(vh, vl) \
            _mm256_insertf128_pd(_mm256_castpd128_pd256(vl), (vh), 1)
        #define _mm256_setr_m128d(vh, vl) \
            _mm256_insertf128_pd(_mm256_castpd128_pd256(vh), (vl), 1)
    // clang-format on
    #endif
#endif

/*
  Explanation of macros dealing with complex math:
   C_MUL(m,a,b)         : m = a*b
   C_SUB( res, a,b)     : res = a - b
   C_SUBFROM( res , a)  : res -= a
   C_ADDTO( res , a)    : res += a
 * */
// clang-format off
#define C_MUL_AVX1(M, A, B)                                                                                          \
    do {                                                                                                             \
        __m256d AB = _mm256_mul_pd(_mm256_permute_pd(A, 0x6),                                                        \
                                   _mm256_permute_pd(B, 0xC));                                                       \
    M   = _mm256_mul_pd(_mm256_addsub_pd(AB,                                                                         \
                                         _mm256_permute2f128_pd(AB, AB, 0x3)),                                       \
                        _mm256_setr_pd(1.0,1.0,-1.0,1.0));                                                           \
    } while (0)
// clang-format on

// #define S_MOD_AVX(A, B) _mm_sub_pd(A, _mm_mul_pd(_mm_div_pd(A, B), B));
//
// #define C_ADD_AVX(M, A, B) M = _mm_add_pd(A, B)
// #define C_SUB_AVX(M, A, B) M = _mm_sub_pd(A, B)
//
// /* Split optional SSE3 functionality */
//
// #define C_MUL_AVX(M, A, B) \
//     do { \
//         __m128d Tms = _mm_move_sd(B, B); \
//         M = _mm_move_sd(B, B); \
//         M = _mm_unpackhi_pd(M, M); \
//         M = _mm_mul_pd(_mm_mul_pd(M, A), _mm_set_pd(-1, 1)); \
//         M = _mm_add_pd(_mm_shuffle_pd(M, M, 0x1), _mm_mul_pd(_mm_unpacklo_pd(Tms, Tms), A)); \
//     } while (0)
