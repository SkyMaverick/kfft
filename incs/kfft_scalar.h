#pragma once

#include "kfft_cpx.h"

// clang-format off
    #define kfft_trace_scalar(fmt, ...)                                                                  \
        kfft_trace("[SCLR]"" " fmt, __VA_ARGS__)
// clang-format on

typedef struct kfft_state {
    kfft_object_t object;

    uint32_t nfft;
    uint32_t flags;

    kfft_plan_cpx* substate;
    kfft_cpx* super_twiddles;
} kfft_plan_sclr;

KFFT_API kfft_plan_sclr*
kfft_config_scalar(const uint32_t nfft, const uint32_t flags, kfft_pool_t* A, size_t* lenmem);
KFFT_API kfft_return_t
kfft_eval_scalar(kfft_plan_sclr* cfg, const kfft_scalar* timedata, kfft_cpx* freqdata);
KFFT_API kfft_return_t
kfft_eval_scalar_norm(kfft_plan_sclr* cfg, const kfft_scalar* timedata, kfft_scalar* data);
KFFT_API kfft_return_t
kfft_evali_scalar(kfft_plan_sclr* cfg, const kfft_cpx* freqdata, kfft_scalar* timedata);

/*Internal use functions (NOT provide KFFT_API) for optimize extensions */
kfft_return_t
kfft_eval_scalar_internal(kfft_plan_sclr* stu, const kfft_scalar* timedata, kfft_cpx* freqdata,
                          kfft_cpx* tmpbuf);
kfft_return_t
kfft_evali_scalar_internal(kfft_plan_sclr* stu, const kfft_cpx* freqdata, kfft_scalar* timedata,
                           kfft_cpx* tmpbuf);
kfft_return_t
kfft_eval_scalar_norm_internal(kfft_plan_sclr* cfg, const kfft_scalar* timedata, kfft_scalar* data,
                               kfft_cpx* fbuf);

#if defined(KFFT_DYNAPI_ENABLE)
// clang-format off
typedef kfft_plan_sclr*
(*kfft_callback_config_scalar)(const uint32_t nfft, const uint32_t flags,
                               kfft_pool_t* A, size_t* lenmem);
typedef kfft_return_t
(*kfft_callback_eval_scalar)(kfft_plan_sclr* cfg, const kfft_scalar* timedata, kfft_cpx* freqdata);
typedef kfft_return_t
(*kfft_callback_eval_scalar_norm)(kfft_plan_sclr* cfg, const kfft_scalar* timedata, kfft_scalar* data);
typedef kfft_return_t
(*kfft_callback_evali_scalar)(kfft_plan_sclr* cfg, const kfft_cpx* freqdata, kfft_scalar* timedata);
// clang-format on
#endif /* KFFT_DYNAPI_ENABLE */
