#pragma once

#include "kfft.h"

typedef struct {
    kfft_object_t object;

    uint32_t nfft;
    uint32_t flags;

    kfft_plan_cpx* plan_fwd;
    kfft_plan_cpx* plan_inv;
} kfft_ccnv_t;

KFFT_API kfft_ccnv_t*
kfft_config_conv_cpx(const uint32_t nfft, const uint32_t flags, kfft_pool_t* A, size_t* lenmem);
KFFT_API kfft_return_t
kfft_eval_conv_cpx(kfft_ccnv_t* plan, const kfft_cpx* fin_A, const kfft_cpx* fin_B, kfft_cpx* fout);

#if defined(KFFT_DYNAPI_ENABLE)
// clang-format off
typedef kfft_ccnv_t*
(*kfft_callback_config_conv_cpx)(const uint32_t nfft, const uint32_t flags, kfft_pool_t* A, size_t* lenmem);
typedef kfft_return_t
(*kfft_callback_eval_conv_cpx)(kfft_ccnv_t* plan, const kfft_cpx* fin_A, const kfft_cpx* fin_B, kfft_cpx* fout);
// clang-format on
#endif /* KFFT_DYNAPI_ENABLE */
