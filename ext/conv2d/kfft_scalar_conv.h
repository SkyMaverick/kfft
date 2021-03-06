#pragma once

#include "kfft.h"

typedef struct {
    kfft_object_t object;

    uint32_t nfft;
    uint32_t x;
    uint32_t y;

    uint32_t flags;

    kfft_plan_s2d* plan_fwd;
    kfft_plan_s2d* plan_inv;
} kfft_plan_s2cnv;

KFFT_API kfft_plan_s2cnv*
kfft_config2_conv_scalar(const uint32_t x, const uint32_t y, const uint32_t flags, kfft_pool_t* A,
                         size_t* lenmem);
KFFT_API kfft_return_t
kfft_eval2_conv_scalar(kfft_plan_s2cnv* plan, const kfft_scalar* fin_A, const kfft_scalar* fin_B,
                       kfft_scalar* fout);

#if defined(KFFT_DYNAPI_ENABLE)
// clang-format off
typedef kfft_plan_s2cnv*
(*kfft_callback_config_conv2_scalar)(const uint32_t x, const uint32_t y, const uint32_t flags,
                                    kfft_pool_t* A, size_t* lenmem);
typedef kfft_return_t
(*kfft_callback_eval_conv2_scalar)(kfft_plan_s2cnv* plan, const kfft_scalar* fin_A,
                                  const kfft_scalar* fin_B, kfft_scalar* fout);
// clang-format on
#endif /* KFFT_DYNAPI_ENABLE */
