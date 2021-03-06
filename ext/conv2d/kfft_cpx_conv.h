#pragma once

#include "kfft.h"

typedef struct {
    kfft_object_t object;

    uint32_t nfft;
    uint32_t x;
    uint32_t y;

    uint32_t flags;

    kfft_plan_c2d* plan_fwd;
    kfft_plan_c2d* plan_inv;
} kfft_plan_c2cnv;

KFFT_API kfft_plan_c2cnv*
kfft_config2_conv_cpx(const uint32_t x, const uint32_t y, const uint32_t flags, kfft_pool_t* A,
                      size_t* lenmem);
KFFT_API kfft_return_t
kfft_eval2_conv_cpx(kfft_plan_c2cnv* plan, const kfft_cpx* fin_A, const kfft_cpx* fin_B,
                    kfft_cpx* fout);

#if defined(KFFT_DYNAPI_ENABLE)
// clang-format off
typedef kfft_plan_c2cnv*
(*kfft_callback_config_conv2_cpx)(const uint32_t x, const uint32_t y, const uint32_t flags, kfft_pool_t* A, size_t* lenmem);
typedef kfft_return_t
(*kfft_callback_eval_conv2_cpx)(kfft_plan_c2cnv* plan, const kfft_cpx* fin_A, const kfft_cpx* fin_B, kfft_cpx* fout);
// clang-format on
#endif /* KFFT_DYNAPI_ENABLE */
