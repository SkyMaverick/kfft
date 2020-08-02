#pragma once

typedef struct {
    kfft_object_t object;
    uint32_t nfft;
    uint32_t flags;
    /* Dimensions count. Ex.: x1,y1,z1,x2,y2,z2 ... xn,yn,zn */
    uint32_t dims;
    /* Step. Ex.: x1,y1,z1,a,a,a,x2,y2,z2,b,b,b, ... ,xn,yn,zn */
    uint32_t step;

    kfft_plan_sclr* subst;
} kfft_plan_ssparse;

KFFT_API kfft_plan_ssparse*
kfft_config_sparse_scalar(const uint32_t nfft, const uint32_t flags, const uint32_t dims,
                          uint32_t step, kfft_pool_t* A, size_t* lenmem);
KFFT_API kfft_return_t
kfft_eval_sparse_scalar(kfft_plan_ssparse* cfg, const kfft_scalar* fin, kfft_cpx* fout);
KFFT_API kfft_return_t
kfft_eval_sparse_scalar_norm(kfft_plan_ssparse* plan, const kfft_scalar* fin, kfft_scalar* fout);
KFFT_API kfft_return_t
kfft_evali_sparse_scalar(kfft_plan_ssparse* cfg, const kfft_cpx* fin, kfft_scalar* fout);
KFFT_API void
kfft_shift_sparse_scalar(kfft_scalar* buf, kfft_scalar* ftmp, const uint32_t nfft,
                         const uint32_t dims, uint32_t step, const bool is_inverse,
                         kfft_pool_t* mmgr);

#if defined(KFFT_DYNAPI_ENABLE)
// clang-format off
typedef kfft_plan_ssparse*
(*kfft_callback_config_sparse_scalar)(const uint32_t nfft, const uint32_t flags,
        const uint32_t dims, uint32_t step, kfft_pool_t* A, size_t* lenmem);
typedef kfft_return_t
(*kfft_callback_eval_sparse_scalar)(kfft_plan_ssparse* cfg, const kfft_scalar* fin, kfft_cpx* fout);
typedef kfft_return_t
(*kfft_callback_eval_sparse_scalar_norm)(kfft_plan_ssparse* plan, const kfft_scalar* fin, kfft_scalar* fout);
typedef kfft_return_t
(*kfft_callback_evali_sparse_scalar)(kfft_plan_ssparse* cfg, const kfft_cpx* fin, kfft_scalar* fout);
typedef void (*kfft_callback_shift_sparse_scalar)(kfft_scalar* buf, kfft_scalar* ftmp,
        const uint32_t nfft, const uint32_t dims, uint32_t step, const bool is_inverse, 
        kfft_pool_t* mmgr);
// clang-format on
#endif /* KFFT_DYNAPI_ENABLE */
