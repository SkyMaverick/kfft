#pragma once

#include "kfft_cpx.h"

typedef struct kfft_state {
    kfft_pool_t* mmgr;

    kfft_comp_t* substate;
    kfft_cpx* tmpbuf;
    // TODO memless
    kfft_cpx* super_twiddles;
#ifdef KFFT_USE_SIMD
    void* pad;
#endif
} kfft_real_t;

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

KFFT_API kfft_real_t*
kfft_config_real(const uint32_t nfft, const uint32_t flags, const uintptr_t A, size_t* lenmem);
KFFT_API void
kfft_eval_real(kfft_real_t* cfg, const kfft_scalar* timedata, kfft_cpx* freqdata);
KFFT_API void
kfft_evali_real(kfft_real_t* cfg, const kfft_cpx* freqdata, kfft_scalar* timedata);
KFFT_API uint32_t
kfft_next_fast_size(uint32_t n);
KFFT_API void
kfft_free(kfft_real_t** cfg);

static inline bool
kfft_isnull(kfft_real_t* in) {
    return (in == 0) ? true : false;
}
