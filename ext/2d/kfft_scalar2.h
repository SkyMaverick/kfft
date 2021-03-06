#pragma once

typedef struct {
    kfft_object_t object;

    uint32_t nfft, x, y;
    uint32_t flags;

    kfft_object_t* plan_x;
    kfft_object_t* plan_y;
} kfft_plan_s2d;

KFFT_API kfft_plan_s2d*
kfft_config2_scalar(const uint32_t x_size, const uint32_t y_size, const uint32_t flags,
                    kfft_pool_t* A, size_t* lenmem);
KFFT_API kfft_return_t
kfft_eval2_scalar(kfft_plan_s2d* plan, const kfft_scalar* fin, kfft_cpx* fout);
KFFT_API kfft_return_t
kfft_evali2_scalar(kfft_plan_s2d* plan, const kfft_cpx* fin, kfft_scalar* fout);
KFFT_API void
kfft_shift2_scalar(kfft_scalar* buf, kfft_scalar* ftmp, const uint32_t sz_x, const uint32_t sz_y,
                   const bool is_inverse, kfft_pool_t* mmgr);

#if defined(KFFT_DYNAPI_ENABLE)
// clang-format off
typedef kfft_plan_s2d*
(*kfft_callback_config2_scalar)(const uint32_t x_size, const uint32_t y_size,
                                const uint32_t flags, kfft_pool_t* A, size_t* lenmem);
typedef kfft_return_t
(*kfft_callback_eval2_scalar)(kfft_plan_s2d* plan, const kfft_scalar* fin, kfft_cpx* fout);
typedef kfft_return_t
(*kfft_callback_evali2_scalar)(kfft_plan_s2d* plan, const kfft_cpx* fin, kfft_scalar* fout);
typedef void
(*kfft_callback_shift2_scalar)(kfft_scalar* buf, kfft_scalar* ftmp, const uint32_t sz_x,
                               const uint32_t sz_y, const bool is_inverse, kfft_pool_t* mmgr);
// clang-format on
#endif /* KFFT_DYNAPI_ENABLE */
