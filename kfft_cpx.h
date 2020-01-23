#pragma once

#include "kfft_guts.h"

#ifndef KFFT_MEMLESS_MODE
    #define SUPER_TWIDDLE(i, P) P->super_twiddles[i]
#else
static inline kfft_cpx
get_super_twiddle(uint32_t i, kfft_plan_t* P) {
    kfft_cpx ret;

    kfft_scalar phase = -KFFT_CONST_PI * ((kfft_scalar)(i + 1) / P->substate->nfft + .5);
    if (P->substate->flags & KFFT_FLAG_INVERSE)
        phase *= -1;
    kf_cexp(&ret, phase);
    return ret;
}
    #define SUPER_TWIDDLE(i, P) get_super_twiddle(i, P)
#endif

#ifdef KFFT_RADER_ALGO
static inline uint32_t
_kfr_power(uint32_t x, uint32_t y, uint32_t m) {
    if (y == 0)
        return 1;
    uint64_t p = _kfr_power(x, y / 2, m) % m;
    p = (p * p) % m;

    return (y % 2 == 0) ? (uint32_t)p : (uint32_t)((x * p) % m);
}
#endif

KFFT_API kfft_kplan_t*
kfft_config_cpx(const uint32_t nfft, const uint32_t flags, const uint8_t level, kfft_pool_t* A,
                size_t* lenmem);

KFFT_API void
kfft_eval_cpx(kfft_kplan_t* cfg, const kfft_cpx* fin, kfft_cpx* fout);

#define kfft_kfree(X)                                                                              \
    do {                                                                                           \
        free(X);                                                                                   \
        X = NULL;                                                                                  \
    } while (0)
