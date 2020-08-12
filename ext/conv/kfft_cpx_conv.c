#include "kfft.h"
#include "kfft_trace.h"

// clang-format off
#define kfft_trace_ccnv(fmt, ...)                                                           \
    kfft_trace("[CNV_CPX]"" " fmt, __VA_ARGS__)
// clang-format on

#if defined(KFFT_TRACE)
static void
kfft_trace_plan(kfft_plan_ccnv* P) {
    kfft_trace_ccnv("%s: %p", "Create KFFT convolution complex plan", (void*)P);
    kfft_trace("\n\t %s - %u", "Total lenght", P->nfft);
    kfft_trace("\n\t %s - %p", "Plan for one lenght", (void*)P->plan_fwd);
    kfft_trace("\n\t %s - %p\n", "Plan inverse for one lenght", (void*)P->plan_inv);
}
#endif /*KFFT_TRACE */
static inline kfft_return_t
kfft_init(kfft_plan_ccnv* st) {
    KFFT_OMP(omp parallel sections shared(st)) {
        KFFT_OMP(omp section) {
            st->plan_fwd =
                kfft_config_cpx(st->nfft, KFFT_CHECK_FLAGS(st->flags), KFFT_PLAN_MMGR(st), NULL);
        }
        KFFT_OMP(omp section) {
            st->plan_inv =
                kfft_config_cpx(st->nfft, KFFT_CHECK_FLAGS(st->flags | KFFT_FLAG_INVERSE),
                                KFFT_PLAN_MMGR(st), NULL);
        }
    }
    return ((st->plan_fwd) && (st->plan_inv)) ? KFFT_RET_SUCCESS : KFFT_RET_ALLOC_FAIL;
}

static inline size_t
kfft_calculate(const uint32_t nfft, const uint32_t flags) {
    size_t r1,r2,ret = sizeof(kfft_plan_ccnv);
    r1 = r2 = 0;

    KFFT_OMP(omp parallel sections ) {
        KFFT_OMP(omp section) {
            kfft_config_cpx(nfft, KFFT_CHECK_FLAGS(flags), NULL, &r1);
        }
        KFFT_OMP(omp section) {
            kfft_config_cpx(nfft, KFFT_CHECK_FLAGS(flags | KFFT_FLAG_INVERSE), NULL, &r2);
        }
    }
    ret += r1 + r2;

    return ret;
}

KFFT_API kfft_plan_ccnv*
kfft_config_conv_cpx(const uint32_t nfft, const uint32_t flags, kfft_pool_t* A, size_t* lenmem) {
    kfft_plan_ccnv* st = NULL;
    size_t memneeded = kfft_calculate(nfft, flags);
    KFFT_ALGO_PLAN_PREPARE(st, flags, kfft_plan_ccnv, memneeded, A, lenmem);
    if (st) {
        st->nfft = nfft;
        st->flags = flags;
        if (kfft_init(st) != KFFT_RET_SUCCESS) {
            KFFT_ALGO_PLAN_TERMINATE(st, A);
            return NULL;
        }
#ifdef KFFT_TRACE
        kfft_trace_plan(st);
#endif
    }
    return st;
}

KFFT_API kfft_return_t
kfft_eval_conv_cpx(kfft_plan_ccnv* plan, const kfft_cpx* fin_A, const kfft_cpx* fin_B,
                   kfft_cpx* fout) {

    kfft_return_t ret, retA, retB;
    ret = retA = retB = KFFT_RET_BUFFER_FAIL;

    kfft_cpx* bufA = KFFT_TMP_ALLOC(sizeof(kfft_cpx) * plan->nfft, KFFT_PLAN_ALIGN(plan));
    if (bufA) {
        kfft_cpx* bufB = (fout == fin_B)
                             ? KFFT_TMP_ALLOC(sizeof(kfft_cpx) * plan->nfft, KFFT_PLAN_ALIGN(plan))
                             : fout;
        if (bufB) {
            KFFT_OMP(omp parallel sections shared(plan)) {
                KFFT_OMP(omp section) { retA = kfft_eval_cpx(plan->plan_fwd, fin_A, bufA); }
                KFFT_OMP(omp section) { retB = kfft_eval_cpx(plan->plan_fwd, fin_B, bufB); }
            }
            if ((retA == KFFT_RET_SUCCESS) && (retB == KFFT_RET_SUCCESS)) {
                VEXFUNC(plan, kfft_math_adamar_cpx, bufB, bufA, plan->nfft);
                ret = kfft_eval_cpx(plan->plan_inv, bufB, fout);
            } else {
                ret = (retA != KFFT_RET_SUCCESS) ? retA : retB;
            }
            if (fout == fin_B)
                KFFT_TMP_FREE(bufB, KFFT_PLAN_ALIGN(plan));
        }
        KFFT_TMP_FREE(bufA, KFFT_PLAN_ALIGN(plan));
    }
    return ret;
}
