#pragma once

#ifndef USE_SYSMATH
    #define KFFT_CONST_PI 3.141592653589793238462643383279502884197169399375105820974944
#else
    #include <math.h>
    #define KFFT_CONST_PI M_PI
#endif

/*
  Explanation of macros dealing with complex math:

   C_MUL(m,a,b)         : m = a*b
   C_SUB( res, a,b)     : res = a - b
   C_SUBFROM( res , a)  : res -= a
   C_ADDTO( res , a)    : res += a
 * */

#define S_MUL(a, b) ((a) * (b))
#define S_DIV(a, b) ((a) / (b))

#define C_CPY(m, a)                                                                                \
    do {                                                                                           \
        m.r = a.r;                                                                                 \
        m.i = a.i;                                                                                 \
    } while (0)

#define C_MUL(m, a, b)                                                                             \
    do {                                                                                           \
        (m).r = (a).r * (b).r - (a).i * (b).i;                                                     \
        (m).i = (a).r * (b).i + (a).i * (b).r;                                                     \
    } while (0)
#define C_MULBYSCALAR(c, s)                                                                        \
    do {                                                                                           \
        (c).r *= (s);                                                                              \
        (c).i *= (s);                                                                              \
    } while (0)
#define C_DIVBYSCALAR(c, s)                                                                        \
    do {                                                                                           \
        (c).r /= (s);                                                                              \
        (c).i /= (s);                                                                              \
    } while (0)

#define C_ADD(res, a, b)                                                                           \
    do {                                                                                           \
        (res).r = (a).r + (b).r;                                                                   \
        (res).i = (a).i + (b).i;                                                                   \
    } while (0)
#define C_SUB(res, a, b)                                                                           \
    do {                                                                                           \
        (res).r = (a).r - (b).r;                                                                   \
        (res).i = (a).i - (b).i;                                                                   \
    } while (0)
#define C_ADDTO(res, a)                                                                            \
    do {                                                                                           \
        (res).r += (a).r;                                                                          \
        (res).i += (a).i;                                                                          \
    } while (0)

#define C_SUBFROM(res, a)                                                                          \
    do {                                                                                           \
        (res).r -= (a).r;                                                                          \
        (res).i -= (a).i;                                                                          \
    } while (0)

#if defined(KFFT_USE_SIMD)
    #define KFFT_COS(phase) _mm_set1_ps(cos(phase))
    #define KFFT_SIN(phase) _mm_set1_ps(sin(phase))
    #define HALF_OF(x) ((x)*_mm_set1_ps(.5))
#else
    #define KFFT_COS(phase) (kfft_scalar) cos(phase)
    #define KFFT_SIN(phase) (kfft_scalar) sin(phase)
    #define HALF_OF(x) ((x)*.5)
#endif

#define kf_cexp(x, phase)                                                                          \
    do {                                                                                           \
        (x)->r = KFFT_COS(phase);                                                                  \
        (x)->i = KFFT_SIN(phase);                                                                  \
    } while (0)
