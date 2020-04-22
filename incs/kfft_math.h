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

#define S_MUL(a, b) ((a) * (b))
#define S_DIV(a, b) ((a) / (b))

#define C_CPY(m, a)                                                                                \
    do {                                                                                           \
        (m).r = (a).r;                                                                             \
        (m).i = (a).i;                                                                             \
    } while (0)
#define C_SWAP(m, a, b)                                                                            \
    do {                                                                                           \
        C_CPY((m), (a));                                                                           \
        C_CPY((a), (b));                                                                           \
        C_CPY((b), (m));                                                                           \
    } while (0)
#define S_SWAP(m, a, b)                                                                            \
    do {                                                                                           \
        (m) = (a);                                                                                 \
        (a) = (b);                                                                                 \
        (b) = (m);                                                                                 \
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

// #if defined(KFFT_USE_SIMD)
//     #define KFFT_COS(phase) _mm_set1_ps(cos(phase))
//     #define KFFT_SIN(phase) _mm_set1_ps(sin(phase))
//     #define HALF_OF(x) ((x)*_mm_set1_ps(.5))
// #else
    #define KFFT_COS(phase) (kfft_scalar) cos(phase)
    #define KFFT_SIN(phase) (kfft_scalar) sin(phase)
    #define HALF_OF(x) ((x)*.5)
// #endif

#define kf_cexp(x, phase)                                                                          \
    do {                                                                                           \
        (x)->r = KFFT_COS(phase);                                                                  \
        (x)->i = KFFT_SIN(phase);                                                                  \
    } while (0)

/* Define as static inline because it's very hot functions */
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
    if (a == 0)
        return b;
    return kfft_math_gcd(b % a, a);
}

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
