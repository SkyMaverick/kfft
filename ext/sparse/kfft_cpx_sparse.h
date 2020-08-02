#pragma once

typedef struct {
    kfft_object_t object;
    uint32_t nfft;
    uint32_t dnfft;

    uint32_t flags;
    /* Dimensions count. Ex.: x1,y1,z1,x2,y2,z2 ... xn,yn,zn */
    uint32_t dims;
    /* Step. Ex.: x1,y1,z1,a,a,a,x2,y2,z2,b,b,b, ... ,xn,yn,zn */
    uint32_t step;

    kfft_plan_cpx* subst;
} kfft_csparse_t;

KFFT_API kfft_csparse_t*
kfft_config_sparse_cpx(const uint32_t nfft, const uint32_t flags, const uint32_t dims,
                       uint32_t step, kfft_pool_t* A, size_t* lenmem);
KFFT_API kfft_return_t
kfft_eval_sparse_cpx(kfft_csparse_t* cfg, const kfft_cpx* fin, kfft_cpx* fout);
KFFT_API void
kfft_shift_sparse_cpx(kfft_cpx* buf, kfft_cpx* ftmp, const uint32_t nfft, const uint32_t dims,
                      uint32_t step, const bool is_inverse, kfft_pool_t* mmgr);

#if defined(KFFT_DYNAPI_ENABLE)
// clang-format off
typedef kfft_csparse_t*
(*kfft_callback_config_sparse_cpx)(const uint32_t nfft, const uint32_t flags, const uint32_t dims,
                                   uint32_t step, kfft_pool_t* A, size_t* lenmem);
typedef kfft_return_t
(*kfft_callback_eval_sparse_cpx)(kfft_csparse_t* cfg, const kfft_cpx* fin, kfft_cpx* fout);
typedef void
(*kfft_callback_shift_sparse_cpx)(kfft_cpx* buf, kfft_cpx* ftmp, const uint32_t nfft, 
                                  const uint32_t dims, uint32_t step, const bool is_inverse,
                                  kfft_pool_t* mmgr);
// clang-format on
#endif /* KFFT_DYNAPI_ENABLE */
