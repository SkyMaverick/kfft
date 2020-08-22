#pragma once

/*!
    \file
    \brief Simplificate equivalent KFFT API

    API for quick and easy use KFFT library.
    Used when there is no need for flexible manipulation of memory and plan objects.
 */

#include "kfft.h"

#ifdef __cplusplus
extern "C" {
#endif

/// KFFT plan type definitions
enum ktiny_plan_type {
    KFFT_PLAN_COMPLEX = 0x0000,            ///< 1D complex plan
    KFFT_PLAN_SCALAR = 0x0001,             ///< 1D scalar plan
    KFFT_PLAN_SCALAR_NORM = 0x0002,        ///< 1D scalar normalized-out plan
    KFFT_PLAN_COMPLEX_2D = 0x0003,         ///< 2D complex plan
    KFFT_PLAN_SCALAR_2D = 0x0004,          ///< 2D scalar plan
    KFFT_PLAN_SCALAR_2D_NORM = 0x0005,     ///< 2D scalar normalized-out plan
    KFFT_PLAN_COMPLEX_SPARSE = 0x0006,     ///< 1D Sparse complex plan
    KFFT_PLAN_SCALAR_SPARSE = 0x0007,      ///< 1D Sparse scalar plan
    KFFT_PLAN_SCALAR_SPARSE_NORM = 0x0008, ///< 1D Sparse scalar normalized-out plan
    KFFT_PLAN_COMPLEX_CONV = 0x0009,       ///< 1D complex convolution plan
    KFFT_PLAN_SCALAR_CONV = 0x000A,        ///< 1D scalar convolution plan
    KFFT_PLAN_COMPLEX_CONV2D = 0x000B,     ///< 2D complex convolution plan
    KFFT_PLAN_SCALAR_CONV2D = 0x000C,      ///< 2D scalar convolution plan
};
/*!
    Unified amalgamatied plan structure
 */
typedef struct {
    uintptr_t state; ///< KFFT plan unified pointer
    int type;        ///< KFFT plan type (::ktiny_plan_type)
} kfft_plan_t;

/// Tiny plan pointer type
typedef kfft_plan_t* kfft_plan;

/// Cast tiny plan pointer to KFFT plan types
#define KTINY_CAST(type, X) (kfft_plan_##type*)((X)->state)
/// Get flags in tiny plan
#define KTINY_FLAGS(type, X) (KTINY_CAST(type, (X)))->flags

/// Cast any pointer to ::kfft_cpx pointer
#define KTINY_CPX(X) (kfft_cpx*)(X)
/// Cast any pointer to ::kfft_scalar pointer
#define KTINY_SCLR(X) (kfft_scalar*)(X)
/*! Amalgamited standart config arguments type. @see kfft_config_cpx */
typedef struct {
    uint32_t nfft; ///< sequencce lenght
} kfft_args_norm;

/*! Amalgamited 2D config arguments type. @see kfft_config2_cpx */
typedef struct {
    uint32_t x; ///< sequence X-lenght
    uint32_t y; ///< sequence Y-lenght
} kfft_args_2d;

/*! Amalgamited sparse config arguments type. @see kfft_config_sparse_cpx */
typedef struct {
    uint32_t nfft; ///< sequence lenght
    uint32_t dims; ///< columns count
    uint32_t step; ///< step size
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

/*!
    Amalgamited kfft_confg_* functions.

    \param[in] T - tiny plan type ::ktiny_plan_type
    \param[in] F - config flags ::kfft_eval_mods
    \param[in] A - unify pointer for kfft_args_* structure
    \result - configured untiped plan ::kfft_plan

    \note Flag ::KFFT_FLAG_RENEW ignored
    \warning This function allocates memory. Must to use ::kfft_finalization() to relase this
   memory.
 */
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
        if (KTINY_FLAGS(s2d, plan) & KFFT_FLAG_INVERSE) {
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
        if (KTINY_FLAGS(ssparse, plan) & KFFT_FLAG_INVERSE) {
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

/*!
    Amalgamited kfft_eval_* functions.

    \param[in] - P ::kfft_plan unify pointer
    \param[in] - Fi input sequence pointer
    \param[in] - Fo output sequence pointer
        \note output buffer must be allocation with ::kfft_malloc or ::KFFT_MALLOC functions.
    \result - return status info ::kfft_return_t

    \warning The function does not check the type and correctness of the data buffer.
 */
#define kfft_evalation(P, Fi, Fo) kfft_tiny_eval(P, (kfft_scalar*)(Fi), (kfft_scalar*)(Fo))

static kfft_return_t KFFT_UNUSED_FUNC
kfft_tiny_conv(const kfft_plan plan, const kfft_scalar* fin1, const kfft_scalar* fin2,
               kfft_scalar* fout) {
    switch (plan->type) {
#if defined (KFFT_CONV_ENABLE)
    case KFFT_PLAN_COMPLEX_CONV:
        return kfft_eval_conv_cpx(KTINY_CAST(ccnv, plan), KTINY_CPX(fin1), KTINY_CPX(fin2),
                                  KTINY_CPX(fout));
    case KFFT_PLAN_SCALAR_CONV:
        return kfft_eval_conv_scalar(KTINY_CAST(scnv, plan), fin1, fin2, fout);
#endif /* KFFT_CONV_ENABLE */
#if defined (KFFT_CONV2D_ENABLE)
    case KFFT_PLAN_COMPLEX_CONV2D:
        return kfft_eval2_conv_cpx(KTINY_CAST(c2cnv, plan), KTINY_CPX(fin1), KTINY_CPX(fin2),
                                   KTINY_CPX(fout));
    case KFFT_PLAN_SCALAR_CONV2D:
        return kfft_eval2_conv_scalar(KTINY_CAST(s2cnv, plan), fin1, fin2, fout);
#endif /* KFFT_CONV2D_ENABLE */
    };
    return KFFT_RET_IMPROPER_PLAN;
}

/*!
    Amalgamited convolution functions.

    \param[in] - P ::kfft_plan unify pointer
    \param[in] - Fi1 input sequence pointer
    \param[in] - Fi2 input sequence pointer
    \param[in] - Fo output sequence pointer
        \note output buffer must be allocation with ::kfft_malloc or ::KFFT_MALLOC functions.

    \result - return status info ::kfft_return_t

    \warning The function does not check the type and correctness of the data buffer.
 */
#define kfft_convolution(P, Fi, Fo)                                                                \
    kfft_tiny_conv(P, (kfft_scalar*)(Fi1), (kfft_scalar*)(Fi2), (kfft_scalar*)(Fo))

static void KFFT_UNUSED_FUNC
kfft_tiny_release(kfft_plan* plan) {
    kfft_cleanup((void*)((*plan)->state));
    kfft_free(*plan);
}

/*!
    Cleanup and free ::kfft_plan. NULL input pointer.

    \param[in] - P ::kfft_plan
    \result - None
 */
#define kfft_finalization(P) kfft_tiny_release(&(P))

#ifdef __cplusplus
}
#endif
