#pragma once

#include "kfft.h"

#ifdef __cplusplus
extern "C" {
#endif

enum {
    KFFT_PLAN_COMPLEX,
    KFFT_PLAN_SCALAR,
    KFFT_PLAN_SCALAR_NORM,
    KFFT_PLAN_COMPLEX_2D,
    KFFT_PLAN_SCALAR_2D,
    KFFT_PLAN_SCALAR_2D_NORM,
    KFFT_PLAN_COMPLEX_SPARSE,
    KFFT_PLAN_SCALAR_SPARSE,
    KFFT_PLAN_SCALAR_SPARSE_NORM,
    KFFT_PLAN_COMPLEX_CONV,
    KFFT_PLAN_SCALAR_CONV,
    KFFT_PLAN_COMPLEX_CONV2D,
    KFFT_PLAN_SCALAR_CONV2D,
};

typedef struct kfft_plan_s {
    uintptr_t state;
    int type;
} kfft_plan_t;

typedef kfft_plan_t* kfft_plan;

#define KTINY_CAST(type, X) (kfft_plan_##type*)((X)->state)
#define KTINY_FLAGS(type, X) (KTINY_CAST(type, (X)))->flags

#define KTINY_CPX(X) (kfft_cpx*)(X)
#define KTINY_SCLR(X) (kfft_scalar*)(X)

typedef struct {
    uint32_t nfft;
} kfft_args_norm;

typedef struct {
    uint32_t x;
    uint32_t y;
} kfft_args_2d;

typedef struct {
    uint32_t nfft;
    uint32_t dims;
    uint32_t step;
} kfft_args_sparse;

static kfft_plan KFFT_UNUSED_FUNC
kfft_tiny_config(unsigned type, uint32_t flags, uintptr_t args) {
    kfft_plan plan = kfft_malloc(sizeof(kfft_plan_t));
    if (plan) {
        switch (type) {

        case KFFT_PLAN_COMPLEX: {
            kfft_args_norm* A = (kfft_args_norm*)args;
            plan->state = (uintptr_t)kfft_config_cpx(A->nfft, flags, NULL, NULL);
            break;
        }
        case KFFT_PLAN_SCALAR:
        case KFFT_PLAN_SCALAR_NORM: {
            kfft_args_norm* A = (kfft_args_norm*)args;
            plan->state = (uintptr_t)kfft_config_scalar(A->nfft, flags, NULL, NULL);
            break;
        }
#if defined(KFFT_2D_ENABLE)
        case KFFT_PLAN_COMPLEX_2D: {
            kfft_args_2d* A = (kfft_args_2d*)args;
            plan->state = (uintptr_t)kfft_config2_cpx(A->x, A->y, flags, NULL, NULL);
            break;
        }
        case KFFT_PLAN_SCALAR_2D:
        case KFFT_PLAN_SCALAR_2D_NORM: {
            kfft_args_2d* A = (kfft_args_2d*)args;
            plan->state = (uintptr_t)kfft_config2_scalar(A->x, A->y, flags, NULL, NULL);
            break;
        }
#endif /* KFFT_2D_ENABLE */

#if defined(KFFT_SPARSE_ENABLE)
        case KFFT_PLAN_COMPLEX_SPARSE: {
            kfft_args_sparse* A = (kfft_args_sparse*)args;
            plan->state =
                (uintptr_t)kfft_config_sparse_cpx(A->nfft, flags, A->dims, A->step, NULL, NULL);
            break;
        }
        case KFFT_PLAN_SCALAR_SPARSE:
        case KFFT_PLAN_SCALAR_SPARSE_NORM: {
            kfft_args_sparse* A = (kfft_args_sparse*)args;
            plan->state =
                (uintptr_t)kfft_config_sparse_scalar(A->nfft, flags, A->dims, A->step, NULL, NULL);
            break;
        }
#endif /* KFFT_SPARSE_ENABLE */

#if defined(KFFT_CONV_ENABLE)
        case KFFT_PLAN_COMPLEX_CONV: {
            kfft_args_norm* A = (kfft_args_norm*)args;
            plan->state = (uintptr_t)kfft_config_conv_cpx(A->nfft, flags, NULL, NULL);
            break;
        }
        case KFFT_PLAN_SCALAR_CONV: {
            kfft_args_norm* A = (kfft_args_norm*)args;
            plan->state = (uintptr_t)kfft_config_conv_scalar(A->nfft, flags, NULL, NULL);
            break;
        }
#endif /* KFFT_CONV_ENABLE */

#if defined(KFFT_CONV2D_ENABLE)
        case KFFT_PLAN_COMPLEX_CONV2D: {
            kfft_args_2d* A = (kfft_args_2d*)args;
            plan->state = (uintptr_t)kfft_config2_conv_cpx(A->x, A->y, flags, NULL, NULL);
            break;
        }
        case KFFT_PLAN_SCALAR_CONV2D: {
            kfft_args_2d* A = (kfft_args_2d*)args;
            plan->state = (uintptr_t)kfft_config2_conv_scalar(A->x, A->y, flags, NULL, NULL);
            break;
        }
#endif /* KFFT_CONV2D_ENABLE */

        default:
            return NULL;
        };

        if (plan->state) {
            plan->type = type;
        } else {
            kfft_free_null((void**)&plan);
        }
    };
    return plan;
}
#define kfft_configuration(T, F, A) kfft_tiny_config((T), (F), (uintptr_t)(A))

static kfft_return_t KFFT_UNUSED_FUNC
kfft_tiny_eval(const kfft_plan plan, const kfft_scalar* fin, kfft_scalar* fout) {
    switch (plan->type) {

    case KFFT_PLAN_COMPLEX:
        return kfft_eval_cpx(KTINY_CAST(cpx, plan), KTINY_CPX(fin), KTINY_CPX(fout));

    case KFFT_PLAN_SCALAR:
        if (KTINY_FLAGS(sclr, plan) & KFFT_FLAG_INVERSE) {
            return kfft_evali_scalar(KTINY_CAST(sclr, plan), KTINY_CPX(fin), fout);
        } else {
            return kfft_eval_scalar(KTINY_CAST(sclr, plan), fin, KTINY_CPX(fout));
        }

    case KFFT_PLAN_SCALAR_NORM:
        return kfft_eval_scalar_norm(KTINY_CAST(sclr, plan), fin, fout);

#if defined(KFFT_2D_ENABLE)
    case KFFT_PLAN_COMPLEX_2D:
        return kfft_eval2_cpx(KTINY_CAST(c2d, plan), KTINY_CPX(fin), KTINY_CPX(fout));

    case KFFT_PLAN_SCALAR_2D:
        if (KTINY_FLAGS(sclr, plan) & KFFT_FLAG_INVERSE) {
            return kfft_evali2_scalar(KTINY_CAST(s2d, plan), KTINY_CPX(fin), fout);
        } else {
            return kfft_eval2_scalar(KTINY_CAST(s2d, plan), fin, KTINY_CPX(fout));
        }

    case KFFT_PLAN_SCALAR_2D_NORM:
        return kfft_eval2_scalar_norm(KTINY_CAST(s2d, plan), fin, fout);
#endif /* KFFT_2D_ENABLE */

#if defined(KFFT_SPARSE_ENABLE)
    case KFFT_PLAN_COMPLEX_SPARSE:
        return kfft_eval_sparse_cpx(KTINY_CAST(csparse, plan), KTINY_CPX(fin), KTINY_CPX(fout));

    case KFFT_PLAN_SCALAR_SPARSE:
        if (KTINY_FLAGS(sclr, plan) & KFFT_FLAG_INVERSE) {
            return kfft_evali_sparse_scalar(KTINY_CAST(ssparse, plan), KTINY_CPX(fin), fout);
        } else {
            return kfft_eval_sparse_scalar(KTINY_CAST(ssparse, plan), fin, KTINY_CPX(fout));
        }

    case KFFT_PLAN_SCALAR_SPARSE_NORM:
        return kfft_eval_sparse_scalar_norm(KTINY_CAST(ssparse, plan), fin, fout);
    };
#endif /* KFFT_SPARSE_ENABLE */
    return KFFT_RET_IMPROPER_PLAN;
}
#define kfft_evalation(P, Fi, Fo) kfft_tiny_eval(P, (kfft_scalar*)(Fi), (kfft_scalar*)(Fo))

static kfft_return_t KFFT_UNUSED_FUNC
kfft_tiny_conv(const kfft_plan plan, const kfft_scalar* fin1, const kfft_scalar* fin2,
               kfft_scalar* fout) {
    switch (plan->type) {
    case KFFT_PLAN_COMPLEX_CONV:
        return kfft_eval_conv_cpx(KTINY_CAST(ccnv, plan), KTINY_CPX(fin1), KTINY_CPX(fin2),
                                  KTINY_CPX(fout));
    case KFFT_PLAN_SCALAR_CONV:
        return kfft_eval_conv_scalar(KTINY_CAST(scnv, plan), fin1, fin2, fout);

    case KFFT_PLAN_COMPLEX_CONV2D:
        return kfft_eval2_conv_cpx(KTINY_CAST(c2cnv, plan), KTINY_CPX(fin1), KTINY_CPX(fin2),
                                   KTINY_CPX(fout));
    case KFFT_PLAN_SCALAR_CONV2D:
        return kfft_eval2_conv_scalar(KTINY_CAST(s2cnv, plan), fin1, fin2, fout);
    };
    return KFFT_RET_IMPROPER_PLAN;
}
#define kfft_convolution(P, Fi, Fo)                                                                \
    kfft_tiny_conv(P, (kfft_scalar*)(Fi1), (kfft_scalar*)(Fi2), (kfft_scalar*)(Fo))

static void KFFT_UNUSED_FUNC
kfft_tiny_release(kfft_plan* plan) {
    kfft_cleanup((void*)((*plan)->state));
    kfft_free(*plan);
}
#define kfft_finalization(P) kfft_tiny_release(&(P))

#ifdef __cplusplus
}
#endif
