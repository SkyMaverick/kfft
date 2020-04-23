#pragma once

#include "kfft.h"

static inline uint32_t
FUNC_SSE(kfft_math_modpow)(uint32_t x, uint32_t y, uint32_t m) {
    if (y == 0)
        return 1;
    uint64_t p = FUNC_SSE(kfft_math_modpow)(x, y / 2, m) % m;
    p = (p * p) % m;

    return (y % 2 == 0) ? (uint32_t)p : (uint32_t)((x * p) % m);
}

static inline uint32_t
FUNC_SSE(kfft_math_gcd)(uint32_t a, uint32_t b) {
    if (a == 0)
        return b;
    return FUNC_SSE(kfft_math_gcd)(b % a, a);
}

#if defined(KFFT_RADER_ALGO)
uint32_t FUNC_SSE(kfft_math_prmn)(uint32_t num);
uint32_t FUNC_SSE(kfft_math_prmni)(uint32_t a, uint32_t m);
#endif /* KFFT_RADER_ALGO */

void FUNC_SSE(kfft_math_adamar_cpx)(kfft_cpx* Fout, kfft_cpx* Fin, uint32_t size);
void FUNC_SSE(kfft_math_transpose_cpx)(const kfft_cpx* Fin, kfft_cpx* Fout, const uint32_t x,
                                       const uint32_t y);
void FUNC_SSE(kfft_math_transpose_scalar)(const kfft_scalar* Fin, kfft_scalar* Fout,
                                          const uint32_t x, const uint32_t y);
void FUNC_SSE(kfft_math_transpose_ip_cpx)(kfft_cpx* Fin, const uint32_t x, const uint32_t y);
void FUNC_SSE(kfft_math_transpose_ip_scalar)(kfft_scalar* Fin, const uint32_t x, const uint32_t y);
