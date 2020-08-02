#pragma once

#include "kfft.h"

typedef struct {
    kfft_object_t object;

    uint32_t nfft;
    uint32_t flags;

    kfft_plan_sclr* plan_fwd;
    kfft_plan_sclr* plan_inv;
} kfft_scnv_t;

KFFT_API kfft_scnv_t*
kfft_config_conv_scalar(const uint32_t nfft, const uint32_t flags, kfft_pool_t* A, size_t* lenmem);
KFFT_API kfft_return_t
kfft_eval_conv_scalar(kfft_scnv_t* plan, const kfft_scalar* fin_A, const kfft_scalar* fin_B,
                      kfft_scalar* fout);

#if defined(KFFT_DYNAPI_ENABLE)
// clang-format off
typedef kfft_scnv_t*
(*kfft_callback_config_conv_scalar)(const uint32_t nfft, const uint32_t flags,
                                    kfft_pool_t* A, size_t* lenmem);
typedef kfft_return_t
(*kfft_callback_eval_conv_scalar)(kfft_scnv_t* plan, const kfft_scalar* fin_A,
                                  const kfft_scalar* fin_B, kfft_scalar* fout);
// clang-format on
#endif /* KFFT_DYNAPI_ENABLE */
