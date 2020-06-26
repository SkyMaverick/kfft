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
#define S_MOD_SSE(A, B) _mm_sub_ps(A, _mm_mul_ps(_mm_div_ps(A, B), B));

#define C_ADD_SSE(M, A, B) M = _mm_add_ps(A, B)
#define C_SUB_SSE(M, A, B) M = _mm_sub_ps(A, B)

/* Split optional SSE3 functionality */
#define CLOAD1122(Ptr)                                                                             \
    _mm_shuffle_ps(_mm_load_ps1((float*)(Ptr)), _mm_load_ps1(((float*)(Ptr)) + 1), 0)
#define CLOAD1212(Ptr)                                                                             \
    _mm_unpackhi_ps(_mm_load_ps1((float*)(Ptr)), _mm_load_ps1(((float*)(Ptr)) + 1))

#if defined(KFFT_HAVE_SSE3)
    /* C_MULDUP_SSE use for A,B loaded with _mm_loaddup_ps() func */
    #define C_MUL_SSE(M, A, B)                                                                     \
        do {                                                                                       \
            __m128 aa, bb, ml;                                                                     \
            aa = _mm_shuffle_ps(A, A, 0x52);                                                       \
            bb = _mm_shuffle_ps(B, B, 0x14);                                                       \
                                                                                                   \
            ml = _mm_mul_ps(aa, bb);                                                               \
            M = _mm_addsub_ps(_mm_movelh_ps(ml, ml), _mm_movehl_ps(ml, ml));                       \
        } while (0)

#else /* KFFT_HAVE_SSE3 */
    #define C_MUL_SSE(M, A, B)                                                                     \
        do {                                                                                       \
            __m128 aa, bb, ml;                                                                     \
            aa = _mm_shuffle_ps(A, A, 0x52);                                                       \
            bb = _mm_shuffle_ps(B, B, 0x14);                                                       \
                                                                                                   \
            ml = _mm_mul_ps(aa, bb);                                                               \
            M = _mm_add_ps(_mm_movelh_ps(ml, ml),                                                  \
                           _mm_mul_ps(_mm_movehl_ps(ml, ml), _mm_setr_ps(-1, 1, -1, 1)));          \
        } while (0)

#endif /* KFFT_HAVE_SSE3 */
