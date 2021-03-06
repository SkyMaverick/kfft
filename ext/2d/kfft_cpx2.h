#pragma once

typedef struct {
    kfft_object_t object;

    uint32_t nfft, x, y;
    uint32_t flags;

    kfft_plan_cpx* plan_x;
    kfft_plan_cpx* plan_y;
} kfft_plan_c2d;

KFFT_API kfft_plan_c2d*
kfft_config2_cpx(const uint32_t x_size, const uint32_t y_size, const uint32_t flags, kfft_pool_t* A,
                 size_t* lenmem);
KFFT_API kfft_return_t
kfft_eval2_cpx(kfft_plan_c2d* cfg, const kfft_cpx* fin, kfft_cpx* fout);
KFFT_API void
kfft_shift2_cpx(kfft_cpx* buf, kfft_cpx* ftmp, const uint32_t sz_x, const uint32_t sz_y,
                const bool is_inverse, kfft_pool_t* mmgr);

#if defined(KFFT_DYNAPI_ENABLE)
// clang-format off
typedef kfft_plan_c2d*
(*kfft_callback_config2_cpx)(const uint32_t x_size, const uint32_t y_size, const uint32_t flags,
                             kfft_pool_t* A, size_t* lenmem);
typedef kfft_return_t
(*kfft_callback_eval2_cpx)(kfft_plan_c2d* cfg, const kfft_cpx* fin, kfft_cpx* fout);
typedef void
(*kfft_callback_shift2_cpx)(kfft_cpx* buf, kfft_cpx* ftmp, const uint32_t sz_x,
                            const uint32_t sz_y, const bool is_inverse, kfft_pool_t* mmgr);
// clang-format on
#endif /* KFFT_DYNAPI_ENABLE */
