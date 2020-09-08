/* The guts header contains all the multiplication and addition macros that are defined for
 fixed or floating point complex numbers.  It also delares the kf_ internal functions.
 */
#include "kfft.h"

#include "kfft_math.h"
#include "kfft_trace.h"

static inline kfft_cpx
kfft_sclr_twiddle(uint32_t i, const kfft_plan_sclr* P) {
    kfft_cpx ret;

    kfft_scalar phase = -KFFT_CONST_PI * ((kfft_scalar)(i + 1) / P->nfft + .5);
    if (P->flags & KFFT_FLAG_INVERSE)
        phase *= -1;
    kf_cexp(&ret, phase);
    return ret;
}

#if defined(KFFT_MEMLESS_MODE)
    #define SUPER_TWIDDLE(i, P) kfft_sclr_twiddle(i, P)
#else
    #define SUPER_TWIDDLE(i, P) P->super_twiddles[i]
#endif /* KFFT_MEMLESS_MODE */

#ifdef KFFT_TRACE

static void
kfft_trace_plan(kfft_plan_sclr* P) {
    KFFT_OMP(omp critical(trace_log)) {
        kfft_trace_scalar("%s: %p", "Create KFFT scalar plan", (void*)P);
        kfft_trace_raw("\n\t %s - %u", "nfft", P->nfft);
        kfft_trace_raw("\n\t %s - %u : ", "flags", P->flags);
        kfft_trace_raw("\n\t %s - %p", "Uses complex plan", (void*)(P->basis));
        kfft_trace_raw("\n\t %s - %p\n", "scalar twiddles", (void*)(P->super_twiddles));
    }
}

#endif /* KFFT_TRACE */

/*******************************************************************************
              Non-symmetric signal (use complex buffer transform)
 ******************************************************************************/

static inline size_t
calculate_trivial(const uint32_t nfft, const uint32_t flags) {
    size_t ret = sizeof(kfft_plan_sclr);
#if !defined(KFFT_MEMLESS_MODE)
    ret += sizeof(kfft_cpx) * nfft;
#endif
    size_t subsize = 0;

    kfft_config_cpx(nfft, flags, NULL, &subsize);
    return ret + subsize;
}

KFFT_API kfft_plan_sclr*
config_trivial(const uint32_t nfft, const uint32_t flags, kfft_pool_t* A, size_t* lenmem) {
    kfft_plan_sclr* P = NULL;
    size_t memneeded = calculate_trivial(nfft, flags);

    KFFT_ALGO_PLAN_PREPARE(P, flags, kfft_plan_sclr, memneeded, A, lenmem);

    if (__likely__(P)) {
        P->basis = kfft_config_cpx(nfft, KFFT_CHECK_FLAGS(flags), P->object.mmgr, NULL);
        if (__unlikely__(P->basis == NULL)) {
            KFFT_ALGO_PLAN_TERMINATE(P, A);
            return NULL;
        }

        P->nfft = nfft;
        P->flags = flags;

//#if !defined(KFFT_MEMLESS_MODE)
//        if (__likely__(nfft > 1)) {
//            P->super_twiddles = kfft_pool_alloc(P->object.mmgr, sizeof(kfft_cpx) * nfft);
//            if (__unlikely__(P->super_twiddles == NULL)) {
//                KFFT_ALGO_PLAN_TERMINATE(P, A);
//                return NULL;
//            }
//        }
//        for (uint32_t i = 0; i < nfft; ++i) {
//            P->super_twiddles[i] = kfft_sclr_twiddle(i, P);
//        }
//#endif /* not KFFT_MEMLESS_MODE */
//
#ifdef KFFT_TRACE
        kfft_trace_plan(P);
#endif
    }
    return P;
}

static kfft_return_t
eval_trivial(kfft_plan_sclr* plan, const kfft_scalar* fin, kfft_cpx* fout, kfft_cpx* ftmp) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    kfft_cpx* fbuf =
        (ftmp) ? ftmp : KFFT_TMP_ALLOC(sizeof(kfft_cpx) * plan->nfft, KFFT_PLAN_ALIGN(plan));
    if (fbuf) {
        for (uint32_t i = 0; i < plan->nfft; i++)
            fbuf[i].r = fin[i]; //< COPY to complex buffer

        ret = kfft_eval_cpx(plan->basis, fbuf, fout);

        if (ftmp == NULL)
            KFFT_TMP_FREE(fbuf, KFFT_PLAN_ALIGN(plan));
    } else {
        ret = KFFT_RET_BUFFER_FAIL;
    }
    return ret;
}

static kfft_return_t
evali_trivial(kfft_plan_sclr* plan, const kfft_cpx* fin, kfft_scalar* fout, kfft_cpx* ftmp) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    kfft_cpx* fbuf =
        (ftmp) ? ftmp : KFFT_TMP_ALLOC(sizeof(kfft_cpx) * plan->nfft, KFFT_PLAN_ALIGN(plan));
    if (fbuf) {
        ret = kfft_eval_cpx(plan->basis, fin, fbuf);

        for (uint32_t i = 0; i < plan->nfft; i++)
            fout[i] = fbuf[i].r; //< COPY to scalar buffer

        if (ftmp == NULL)
            KFFT_TMP_FREE(fbuf, KFFT_PLAN_ALIGN(plan));
    } else {
        ret = KFFT_RET_BUFFER_FAIL;
    }
    return ret;
}

/*******************************************************************************
                Symmetric signal (use Nayquist frequences N/2)
 ******************************************************************************/

static kfft_return_t
eval_nayquist(kfft_plan_sclr* plan, const kfft_scalar* fin, kfft_cpx* fout, kfft_cpx* ftmp) {
    kfft_return_t ret = KFFT_RET_SUCCESS;
    return ret;
}

static kfft_return_t
evali_nayquist(kfft_plan_sclr* plan, const kfft_cpx* fin, kfft_scalar* fout, kfft_cpx* ftmp) {
    kfft_return_t ret = KFFT_RET_SUCCESS;
    return ret;
}

static kfft_plan_sclr*
config_nayquist(const uint32_t nfft, const uint32_t flags, kfft_pool_t* A, size_t* lenmem) {
    return NULL;
}

/*******************************************************************************
                            API functionality
 ******************************************************************************/

kfft_return_t
kfft_eval_scalar_internal(kfft_plan_sclr* plan, const kfft_scalar* fin, kfft_cpx* fout,
                          kfft_cpx* ftmp) {

    if (__unlikely__(plan->basis->flags & KFFT_FLAG_INVERSE))
        return KFFT_RET_IMPROPER_PLAN;

    return (plan->nfft % 2) ? eval_trivial(plan, fin, fout, ftmp)
                            : eval_nayquist(plan, fin, fout, ftmp);
}

kfft_return_t
kfft_evali_scalar_internal(kfft_plan_sclr* plan, const kfft_cpx* fin, kfft_scalar* fout,
                           kfft_cpx* ftmp) {

    if (__unlikely__(!(plan->basis->flags & KFFT_FLAG_INVERSE)))
        return KFFT_RET_IMPROPER_PLAN;

    return (plan->nfft % 2) ? evali_trivial(plan, fin, fout, ftmp)
                            : evali_nayquist(plan, fin, fout, ftmp);
}

KFFT_API kfft_plan_sclr*
kfft_config_scalar(const uint32_t nfft, const uint32_t flags, kfft_pool_t* A, size_t* lenmem) {
    return (nfft % 2) ? config_trivial(nfft, flags, A, lenmem)
                      : config_nayquist(nfft, flags, A, lenmem);
}
KFFT_API kfft_return_t
kfft_eval_scalar(kfft_plan_sclr* plan, const kfft_scalar* fin, kfft_cpx* fout) {
    return kfft_eval_scalar_internal(plan, fin, fout, NULL);
}
KFFT_API kfft_return_t
kfft_evali_scalar(kfft_plan_sclr* plan, const kfft_cpx* fin, kfft_scalar* fout) {
    return kfft_evali_scalar_internal(plan, fin, fout, NULL);
}
// KFFT_API kfft_plan_sclr*
// kfft_config_scalar(const uint32_t nfft, const uint32_t flags, kfft_pool_t* A, size_t* lenmem) {
//     kfft_plan_sclr* P = NULL;
//     size_t memneeded = kfft_calculate(nfft, flags);
//
//     KFFT_ALGO_PLAN_PREPARE(P, flags, kfft_plan_sclr, memneeded, A, lenmem);
//
//     if (__likely__(P)) {
//         P->basis = kfft_config_cpx(HALF_NFFT(nfft), KFFT_CHECK_FLAGS(flags), P->object.mmgr,
//         NULL); if (__unlikely__(P->basis == NULL)) {
//             KFFT_ALGO_PLAN_TERMINATE(P, A);
//             return NULL;
//         }
//
//         P->nfft = nfft;
//         P->flags = flags;
//
// #if !defined(KFFT_MEMLESS_MODE)
//         if (__likely__(nfft > 1)) {
//             P->super_twiddles = kfft_pool_alloc(P->object.mmgr, sizeof(kfft_cpx) *
//             (HALF_NFFT(nfft))); if (__unlikely__(P->super_twiddles == NULL)) {
//                 KFFT_ALGO_PLAN_TERMINATE(P, A);
//                 return NULL;
//             }
//         }
//         for (uint32_t i = 0; i < HALF_NFFT(nfft); ++i) {
//             P->super_twiddles[i] = kfft_sclr_twiddle(i, P);
//         }
// #endif /* not KFFT_MEMLESS_MODE */
//
// #ifdef KFFT_TRACE
//         kfft_trace_plan(P);
// #endif
//     }
//     return P;
// }
//
// static inline kfft_return_t
// eval_forward_internal(const kfft_plan_sclr* plan, const kfft_cpx* fin, kfft_cpx* fout) {
//     kfft_return_t ret = KFFT_RET_SUCCESS;
//     kfft_cpx fpnk, fpk, f1k, f2k, tw;
//
//     uint32_t k, ncfft = plan->basis->nfft;
//
//     C_CPY(fout[0], fin[0]);
//
//     for (k = 1; k <= ncfft; ++k) {
//         fpk = fin[k];
//         fpnk.r = fin[ncfft - k].r;
//         fpnk.i = -fin[ncfft - k].i;
//
//         C_ADD(f1k, fpk, fpnk);
//         C_SUB(f2k, fpk, fpnk);
//         C_MUL(tw, f2k, SUPER_TWIDDLE(k - 1, plan) /* P->super_twiddles[k - 1] */);
//
//         fout[k].r = HALF_OF(f1k.r + tw.r);
//         fout[k].i = HALF_OF(f1k.i + tw.i);
//         fout[ncfft - k].r = HALF_OF(f1k.r - tw.r);
//         fout[ncfft - k].i = HALF_OF(tw.i - f1k.i);
//     }
//
//     return ret;
// }
//
// static kfft_return_t
// eval_func(kfft_plan_sclr* plan, kfft_cpx* fout, kfft_cpx* ftmp) {
//     kfft_return_t ret = KFFT_RET_SUCCESS;
//
//     ret = kfft_eval_cpx(plan->basis, fout, ftmp);
//     if (__likely__(ret == KFFT_RET_SUCCESS)) {
//         ret = eval_forward_internal(plan, ftmp, fout);
//     }
//
//     return ret;
// }
//
// kfft_return_t
// kfft_eval_scalar_internal(kfft_plan_sclr* plan, const kfft_scalar* fin, kfft_cpx* fout,
//                           kfft_cpx* ftmp) {
//     kfft_return_t ret = KFFT_RET_SUCCESS;
//     /* input buffer timedata is stored row-wise */
//     kfft_plan_sclr* P = (kfft_plan_sclr*)plan;
//
//     if (__unlikely__(P->basis->flags & KFFT_FLAG_INVERSE)) {
//         return KFFT_RET_IMPROPER_PLAN;
//     }
//
//     uint32_t ncfft = plan->basis->nfft;
//
//     KFFT_ZEROMEM(fout, sizeof(kfft_cpx) * ncfft);
//     memcpy((void*)fout, (void*)fin, sizeof(kfft_scalar) * plan->nfft);
//
//     if (__unlikely__(ftmp == NULL)) {
//         kfft_cpx* tbuf = KFFT_TMP_ALLOC(sizeof(kfft_cpx) * ncfft, KFFT_PLAN_ALIGN(plan));
//         if (__likely__(tbuf)) {
//
//             ret = eval_func(plan, fout, tbuf);
//             KFFT_TMP_FREE(tbuf, KFFT_PLAN_ALIGN(plan));
//         } else {
//             ret = KFFT_RET_BUFFER_FAIL;
//         }
//     } else {
//         ret = eval_func(plan, fout, ftmp);
//     }
//
//     return ret;
// }
//
// KFFT_API kfft_return_t
// kfft_eval_scalar(kfft_plan_sclr* plan, const kfft_scalar* fin, kfft_cpx* fout) {
//     return kfft_eval_scalar_internal(plan, fin, fout, NULL);
// }
//
// static inline kfft_return_t
// eval_inverse_internal(const kfft_plan_sclr* plan, const kfft_cpx* fin, kfft_cpx* fout) {
//     kfft_return_t ret = KFFT_RET_SUCCESS;
// //    kfft_cpx fk, fnkc, fek, fok, tmp;
// //
// //    uint32_t k, ncfft = plan->basis->nfft;
// //
// //    C_CPY(fout[0], fin[0]);
// //    C_MULBYSCALAR(fout[0], 2);
// //
// //    for (k = 1; k <= ncfft; ++k) {
// //        fk = fin[k];
// //        fnkc.r = fin[ncfft - k].r;
// //        fnkc.i = -fin[ncfft - k].i;
// //
// //        C_ADD(fek, fk, fnkc);
// //        C_SUB(tmp, fk, fnkc);
// //        C_MUL(fok, tmp, SUPER_TWIDDLE(k - 1, plan) /* P->super_twiddles[k - 1] */);
// //        C_ADD(fout[k], fek, fok);
// //        C_SUB(fout[ncfft - k], fek, fok);
// //        fout[ncfft - k].i *= -1;
// //    }
//     return ret;
// }
//
// static kfft_return_t
// evali_func(kfft_plan_sclr* plan, const kfft_cpx* fin, kfft_scalar* fout, kfft_cpx* ftmp) {
//     kfft_return_t ret = KFFT_RET_SUCCESS;
// //    uint32_t ncfft = plan->nfft;
// //
// //    KFFT_TMP_ZEROMEM(ftmp, sizeof(kfft_cpx) * ncfft);
// //
// //    ret = eval_inverse_internal(plan, fin, ftmp);
// //    if (__likely__(ret == KFFT_RET_SUCCESS)) {
// //        kfft_cpx* fbuf = ftmp + plan->nfft;
// //        ret = kfft_eval_cpx(plan->basis, ftmp, fbuf);
// //        if (__likely__(ret == KFFT_RET_SUCCESS)) {
// //            for (uint32_t i = 0; i < ncfft; i++) {
// //                fout[i] = fbuf[i].r / 2;
// //            }
// //        }
// //    }
//     return ret;
// }
//
// kfft_return_t
// kfft_evali_scalar_internal(kfft_plan_sclr* plan, const kfft_cpx* fin, kfft_scalar* fout,
//                            kfft_cpx* ftmp) {
//     kfft_return_t ret = KFFT_RET_SUCCESS;
// //    kfft_plan_sclr* P = (kfft_plan_sclr*)plan;
// //
// //    if (__unlikely__(!(P->basis->flags & KFFT_FLAG_INVERSE))) {
// //        return KFFT_RET_IMPROPER_PLAN;
// //    }
// //
// //    uint32_t ncfft = P->nfft;
// //
// //    if (__unlikely__(ftmp == NULL)) {
// //        kfft_cpx* tbuf = KFFT_TMP_ALLOC(2 * sizeof(kfft_cpx) * ncfft, KFFT_PLAN_ALIGN(plan));
// //        if (__likely__(tbuf)) {
// //            ret = evali_func(plan, fin, fout, tbuf);
// //            KFFT_TMP_FREE(tbuf, KFFT_PLAN_ALIGN(plan));
// //        } else {
// //            ret = KFFT_RET_BUFFER_FAIL;
// //        }
// //    } else {
// //        ret = evali_func(plan, fin, fout, ftmp);
// //    }
//     return ret;
// }
//
// KFFT_API kfft_return_t
// kfft_evali_scalar(kfft_plan_sclr* plan, const kfft_cpx* fin, kfft_scalar* fout) {
//     /* input buffer timedata is stored row-wise */
//     return kfft_evali_scalar_internal(plan, fin, fout, NULL);
// }
//
// #undef HALF_NFFT
