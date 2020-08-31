#pragma once

/*!
    \file
    \brief Mathematical primitive functions

    Standart mathematical functions don't implements in standart math library (libm)
    such as complex math, matrix operations, modpow etc.
 */

/// Standart PI constant
#define KFFT_CONST_PI 3.141592653589793238462643383279502884197169399375105820974944

#ifdef KFFT_USE_SYSMATH
    #include <math.h>
#else
    #include "kfft_custom_math.h"
#endif

/*
  Explanation of macros dealing with complex math:
   C_MUL(m,a,b)         : m = a*b
   C_SUB( res, a,b)     : res = a - b
   C_SUBFROM( res , a)  : res -= a
   C_ADDTO( res , a)    : res += a
 * */

#define UNI_SWAP(tmp, A, B)                                                                        \
    do {                                                                                           \
        tmp = A = B = tmp;                                                                         \
    } while (0)

/// Scalar multoplication
#define S_MUL(a, b) ((a) * (b))
/// Scalar division
#define S_DIV(a, b) ((a) / (b))

/// Scalar swap (a<=>b)
#define S_SWAP(a, b)                                                                               \
    do {                                                                                           \
        kfft_scalar tmp = 0;                                                                       \
        UNI_SWAP(tmp, a, b);                                                                       \
    } while (0)

/// Complex swap (A<=>B)
#define C_SWAP(a, b)                                                                               \
    do {                                                                                           \
        S_SWAP(a.r, b.r);                                                                          \
        S_SWAP(a.i, b.i);                                                                          \
    } while (0)

/// Complex copy (A = B)
#define C_CPY(m, a)                                                                                \
    do {                                                                                           \
        (m).r = (a).r;                                                                             \
        (m).i = (a).i;                                                                             \
    } while (0)

/// Complex multyplify (A*B)
#define C_MUL(m, a, b)                                                                             \
    do {                                                                                           \
        (m).r = (a).r * (b).r - (a).i * (b).i;                                                     \
        (m).i = (a).r * (b).i + (a).i * (b).r;                                                     \
    } while (0)

/// Complex mulitiplify with scalar (A * b)
#define C_MULBYSCALAR(c, s)                                                                        \
    do {                                                                                           \
        (c).r *= (s);                                                                              \
        (c).i *= (s);                                                                              \
    } while (0)

/// Complex division by scalar (A / b)
#define C_DIVBYSCALAR(c, s)                                                                        \
    do {                                                                                           \
        (c).r /= (s);                                                                              \
        (c).i /= (s);                                                                              \
    } while (0)

/// Complex summary (A + B)
#define C_ADD(res, a, b)                                                                           \
    do {                                                                                           \
        (res).r = (a).r + (b).r;                                                                   \
        (res).i = (a).i + (b).i;                                                                   \
    } while (0)

/// Complex subtraction (A - B)
#define C_SUB(res, a, b)                                                                           \
    do {                                                                                           \
        (res).r = (a).r - (b).r;                                                                   \
        (res).i = (a).i - (b).i;                                                                   \
    } while (0)

/// Complex summary, result in A (A + B)
#define C_ADDTO(res, a)                                                                            \
    do {                                                                                           \
        (res).r += (a).r;                                                                          \
        (res).i += (a).i;                                                                          \
    } while (0)

/// Complex subtraction, resul in A (A - B)
#define C_SUBFROM(res, a)                                                                          \
    do {                                                                                           \
        (res).r -= (a).r;                                                                          \
        (res).i -= (a).i;                                                                          \
    } while (0)

#ifdef KFFT_USE_SYSMATH
    #define KFFT_SQRT(X) sqrt((X))
#else
    #define KFFT_SQRT(X) kfft_math_sqrt((X))
#endif

#ifdef KFFT_USE_SYSMATH
    #define kf_cexp(x, phase)                                                                      \
        do {                                                                                       \
            (x)->r = (kfft_scalar)cos(phase);                                                      \
            (x)->i = (kfft_scalar)sin(phase);                                                      \
        } while (0)
#else

    #if defined(KFFT_HALF_SCALAR)
static inline unsigned
__kfft_sincos_float(kfft_cpx* X, float num) {
    double co = 0.0, si = 0.0;
    unsigned ret = kfft_sincos_double(&co, &si, (double)num);

    X->r = (float)co;
    X->i = (float)si;

    return ret;
}
        #define kf_cexp(x, phase) __kfft_sincos_float((x), phase)
    #else
        #define kf_cexp(x, phase) kfft_sincos_double(&((x)->r), &((x)->i), phase);
    #endif
#endif
#define HALF_OF(x) ((x)*.5)
// #endif

#define S_MAX(X, Y) ((X) > (Y)) ? (X) : (Y)
#define S_MIN(X, Y) ((X) < (Y)) ? (X) : (Y)
#define S_EQUAL(X, Y) ((X) == (Y)) ? true : false

static inline kfft_scalar
kfft_math_mgnt(const kfft_cpx* A) {
    return KFFT_SQRT(A->r * A->r + A->i * A->i);
}

/// Maximal complex number (by magnitude)
#define C_MAX_ABS(X, Y) (kfft_math_mgnt((X)) > kfft_math_mgnt((Y))) ? (X) : (Y)
/// Minimal complex number (by magnitude)
#define C_MIN_ABS(X, Y) (kfft_math_mgnt((X)) < kfft_math_mgnt((Y))) ? (X) : (Y)
/// Equal complex number (by magnitude)
#define C_EQU_ABS(X, Y) (kfft_math_mgnt((X)) == kfft_math_mgnt((Y))) ? true : false
/// Equal complex number (by R,I)
#define C_EQU(X, Y) (((X.r) == (Y.r)) && ((X.i) == (Y.i))) ? true : false

/* Define as static inline because it's very hot functions */
/*!
    ![Multiplication by modulo N](modpow.svg)
 */
static inline uint32_t
kfft_math_modpow(uint32_t x, uint32_t y, uint32_t m) {
    if (y == 0)
        return 1;
    uint64_t p = kfft_math_modpow(x, y / 2, m) % m;
    p = (p * p) % m;

    return (y % 2 == 0) ? (uint32_t)p : (uint32_t)((x * p) % m);
}

static inline uint32_t
kfft_math_gcd(uint32_t a, uint32_t b) {
    while (b) {
        b ^= a;
        a ^= b;
        b ^= a;

        b = a % b;
    }
    return a;
}
static inline kfft_cpx
kfft_kernel_twiddle(uint32_t i, uint32_t size, bool is_inverse) {
    kfft_cpx ret = {0, 0};

    kfft_scalar phase = -2 * KFFT_CONST_PI * i / size;
    if (is_inverse)
        phase *= -1;

    kf_cexp(&ret, phase);
    return ret;
}

#if defined(KFFT_MEMLESS_MODE)
    #define TWIDDLE(i, P) kfft_kernel_twiddle(i, (P)->nfft, ((P)->flags & KFFT_FLAG_INVERSE))
#else
    #define TWIDDLE(i, P) (P)->twiddles[i]
#endif

#if defined(KFFT_RADER_ALGO)
uint32_t
kfft_math_prmn(uint32_t num);
uint32_t
kfft_math_prmni(uint32_t a, uint32_t m);
#endif /* KFFT_RADER_ALGO */

void
kfft_math_adamar_cpx(kfft_cpx* Fout, kfft_cpx* Fin, uint32_t size);
void
kfft_math_transpose_cpx(const kfft_cpx* Fin, kfft_cpx* Fout, const uint32_t x, const uint32_t y);
void
kfft_math_transpose_scalar(const kfft_scalar* Fin, kfft_scalar* Fout, const uint32_t x,
                           const uint32_t y);
void
kfft_math_transpose_ip_cpx(kfft_cpx* Fin, const uint32_t x, const uint32_t y);
void
kfft_math_transpose_ip_scalar(kfft_scalar* Fin, const uint32_t x, const uint32_t y);
void
kfft_math_magnitude(const kfft_cpx* Fin, kfft_scalar* Fout, uint32_t size);
void
kfft_math_magnitude_ip(kfft_cpx* Fin, uint32_t size);
