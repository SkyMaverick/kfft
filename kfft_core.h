#pragma once

#include "kfft_guts.h"

#ifndef KFFT_MEMLESS_MODE
    #define SUPER_TWIDDLE(i, P) P->super_twiddles[i]
#else
static inline kfft_cpx
get_super_twiddle(int i, kfft_plan_t* P) {
    kfft_cpx ret;

    kfft_scalar phase = -KFFT_CONST_PI * ((kfft_scalar)(i + 1) / P->substate->nfft + .5);
    if (P->substate->inverse)
        phase *= -1;
    kf_cexp(&ret, phase);
    return ret;
}
    #define SUPER_TWIDDLE(i, P) get_super_twiddle(i, P)
#endif

unsigned
kfft_prime_root(unsigned num);

unsigned
kfft_primei_root(unsigned a, unsigned m);

kfft_kplan_t*
kfft_kconfig(int nfft, int inverse_fft, int level, void* mem, size_t* lenmem);

void
__kfft(kfft_kplan_t* cfg, const kfft_cpx* fin, kfft_cpx* fout);

#define kfft_kfree(X)                                                                              \
    do {                                                                                           \
        free(X);                                                                                   \
        X = NULL;                                                                                  \
    } while (0)
