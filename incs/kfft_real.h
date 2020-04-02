#pragma once

#include "kfft_cpx.h"

// clang-format off
    #define kfft_trace_real(fmt, ...)                                                                  \
        kfft_trace("[REAL]"" " fmt, __VA_ARGS__)
// clang-format on

typedef struct kfft_state {
    kfft_object_t object;

    kfft_comp_t* substate;
    kfft_cpx* super_twiddles;
} kfft_real_t;

KFFT_API kfft_real_t*
kfft_config_real(const uint32_t nfft, const uint32_t flags, const kfft_pool_t* A, size_t* lenmem);
KFFT_API kfft_return_t
kfft_eval_real(kfft_real_t* cfg, const kfft_scalar* timedata, kfft_cpx* freqdata);
KFFT_API kfft_return_t
kfft_evali_real(kfft_real_t* cfg, const kfft_cpx* freqdata, kfft_scalar* timedata);
KFFT_API uint32_t
kfft_next_fast_size(uint32_t n);
