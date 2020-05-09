#pragma once

#ifdef KFFT_USE_SYSMATH
    #include <math.h>
#endif
#define KFFT_CONST_PI 3.141592653589793238462643383279502884197169399375105820974944

/*
  Explanation of macros dealing with complex math:
   C_MUL(m,a,b)         : m = a*b
   C_SUB( res, a,b)     : res = a - b
   C_SUBFROM( res , a)  : res -= a
   C_ADDTO( res , a)    : res += a
 * */
#define C_ADD_SSE(M, A, B) M = _mm_add_pd(A, B)
#define C_SUB_SSE(M, A, B) M = _mm_sub_pd(A, B)

/* Split optional SSE3 functionality */

#if defined(KFFT_HAVE_SSE3)
    /* C_MULDUP_SSE use for A,B loaded with _mm_loaddup_pd() func */
    #define C_MULDUP_SSE(M, A, B)                                                                  \
        do {                                                                                       \
            M = _mm_mul_pd(M, A);                                                                  \
            M = _mm_shuffle_pd(M, M, 0x1);                                                         \
            M = _mm_addsub_pd(_mm_mul_pd(_mm_move_sd(B, B), A), M);                                \
        } while (0)

    #define C_MUL_SSE(M, A, B)                                                                     \
        do {                                                                                       \
            __m128d Tms = _mm_move_sd(B, B);                                                       \
            M = _mm_move_sd(B, B);                                                                 \
            M = _mm_unpackhi_pd(M, M);                                                             \
            C_MULDUP_SSE(M, A, _mm_unpacklo_pd(Tms, Tms));                                         \
        } while (0)

#else /* KFFT_HAVE_SSE3 */
    #define C_MUL_SSE(M, A, B)                                                                     \
        do {                                                                                       \
            __m128d Tms = _mm_move_sd(B, B);                                                       \
            M = _mm_move_sd(B, B);                                                                 \
            M = _mm_unpackhi_pd(M, M);                                                             \
            M = _mm_mul_pd(_mm_mul_pd(M, A), _mm_set_pd(-1, 1));                                   \
            M = _mm_add_pd(_mm_shuffle_pd(M, M, 0x1), _mm_mul_pd(_mm_unpacklo_pd(Tms, Tms), A));   \
        } while (0)

#endif /* KFFT_HAVE_SSE3 */
