#pragma once

#include "kfft_cpx.h"

typedef struct kfft_state {
    kfft_pool_t* mmgr;

    kfft_comp_t* substate;
    kfft_cpx* tmpbuf;
    // TODO memless
    kfft_cpx* super_twiddles;
} kfft_real_t;

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
