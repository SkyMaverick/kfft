#pragma once

#include "kfft.h"

typedef struct {
    kfft_object_t object;

    uint32_t nfft;
    uint32_t flags;

    kfft_plan_cpx* plan_fwd;
    kfft_plan_cpx* plan_inv;
} kfft_plan_ccnv;

#define KFFT_CHECK_FLAGS_CCNV(X) (KFFT_CHECK_FLAGS(X) & (~KFFT_FLAG_INVERSE))

KFFT_API kfft_plan_ccnv*
kfft_config_conv_cpx(const uint32_t nfft, const uint32_t flags, kfft_pool_t* A, size_t* lenmem);
KFFT_API kfft_return_t
kfft_eval_conv_cpx(kfft_plan_ccnv* plan, const kfft_cpx* fin_A, const kfft_cpx* fin_B,
                   kfft_cpx* fout);

#if defined(KFFT_DYNAPI_ENABLE)
// clang-format off
typedef kfft_plan_ccnv*
(*kfft_callback_config_conv_cpx)(const uint32_t nfft, const uint32_t flags, kfft_pool_t* A, size_t* lenmem);
typedef kfft_return_t
(*kfft_callback_eval_conv_cpx)(kfft_plan_ccnv* plan, const kfft_cpx* fin_A, const kfft_cpx* fin_B, kfft_cpx* fout);
// clang-format on
#endif /* KFFT_DYNAPI_ENABLE */
