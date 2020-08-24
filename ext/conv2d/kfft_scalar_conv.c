#include "kfft.h"
#include "kfft_trace.h"

// clang-format off
#define kfft_trace_scnv2(fmt, ...)                                                           \
    kfft_trace("[SCNV_2D]"" " fmt, __VA_ARGS__)
// clang-format on

#if defined(KFFT_TRACE)
static void
kfft_trace_plan(kfft_plan_s2cnv* P) {
    KFFT_OMP(omp critical(trace_log)) {
        kfft_trace_scnv2("%s: %p", "Create KFFT convolution scalar plan", (void*)P);
        kfft_trace_raw("\n\t %s - %u", "Total lenght", P->nfft);
        kfft_trace_raw("\n\t %s - %u", "Size X", P->x);
        kfft_trace_raw("\n\t %s - %u", "Size Y", P->y);
        kfft_trace_raw("\n\t %s - %p", "Plan for one lenght", (void*)P->plan_fwd);
        kfft_trace_raw("\n\t %s - %p\n", "Plan inverse for one lenght", (void*)P->plan_inv);
    }
}
#endif /*KFFT_TRACE */
static inline kfft_return_t
kfft_init(kfft_plan_s2cnv* plan) {
    KFFT_OMP(omp parallel sections shared(plan)) {
        KFFT_OMP(omp section) {
            plan->plan_fwd = kfft_config2_scalar(plan->x, plan->y, KFFT_CHECK_FLAGS(plan->flags),
                                                 KFFT_PLAN_MMGR(plan), NULL);
        }
        KFFT_OMP(omp section) {
            plan->plan_inv = kfft_config2_scalar(plan->x, plan->y,
                                                 KFFT_CHECK_FLAGS(plan->flags | KFFT_FLAG_INVERSE),
                                                 KFFT_PLAN_MMGR(plan), NULL);
        }
    }
    return ((plan->plan_fwd) && (plan->plan_inv)) ? KFFT_RET_SUCCESS : KFFT_RET_ALLOC_FAIL;
}

static inline size_t
kfft_calculate(const uint32_t x, const uint32_t y, const uint32_t flags) {
    size_t r1, r2, ret = sizeof(kfft_plan_s2cnv);
    r1 = r2 = 0;

    KFFT_OMP(omp parallel sections) {
        KFFT_OMP(omp section) { kfft_config2_scalar(x, y, KFFT_CHECK_FLAGS(flags), NULL, &r1); }
        KFFT_OMP(omp section) {
            kfft_config2_scalar(x, y, KFFT_CHECK_FLAGS(flags | KFFT_FLAG_INVERSE), NULL, &r2);
        }
    }
    ret += r1 + r2;
    return ret;
}

KFFT_API kfft_plan_s2cnv*
kfft_config2_conv_scalar(const uint32_t x, const uint32_t y, const uint32_t flags, kfft_pool_t* A,
                         size_t* lenmem) {
    kfft_plan_s2cnv* plan = NULL;
    size_t memneeded = kfft_calculate(x, y, flags);
    KFFT_ALGO_PLAN_PREPARE(plan, flags, kfft_plan_s2cnv, memneeded, A, lenmem);
    if (plan) {

        plan->nfft = x * y;
        plan->x = x;
        plan->y = y;
        plan->flags = flags;

        if (kfft_init(plan) != KFFT_RET_SUCCESS) {
            KFFT_ALGO_PLAN_TERMINATE(plan, A);
            return NULL;
        }
#ifdef KFFT_TRACE
        kfft_trace_plan(plan);
#endif
    }
    return plan;
}

KFFT_API kfft_return_t
kfft_eval2_conv_scalar(kfft_plan_s2cnv* plan, const kfft_scalar* fin_A, const kfft_scalar* fin_B,
                       kfft_scalar* fout) {

    kfft_return_t ret, retA, retB;
    ret = retA = retB = KFFT_RET_BUFFER_FAIL;

    kfft_cpx* bufA = KFFT_TMP_ALLOC(sizeof(kfft_cpx) * plan->nfft, KFFT_PLAN_ALIGN(plan));
    if (bufA) {
        kfft_cpx* bufB = KFFT_TMP_ALLOC(sizeof(kfft_cpx) * plan->nfft, KFFT_PLAN_ALIGN(plan));
        if (bufB) {
            KFFT_OMP(omp parallel sections shared(plan)) {
                KFFT_OMP(omp section) { retA = kfft_eval2_scalar(plan->plan_fwd, fin_A, bufA); }
                KFFT_OMP(omp section) { retB = kfft_eval2_scalar(plan->plan_fwd, fin_B, bufB); }
            }
            if ((retA == KFFT_RET_SUCCESS) && (retB == KFFT_RET_SUCCESS)) {
                VEXFUNC(plan, kfft_math_adamar_cpx, bufB, bufA, plan->nfft);
                ret = kfft_evali2_scalar(plan->plan_inv, bufB, fout);
            } else {
                ret = (retA != KFFT_RET_SUCCESS) ? retA : retB;
            }
            KFFT_TMP_FREE(bufB, KFFT_PLAN_ALIGN(plan));
        }
        KFFT_TMP_FREE(bufA, KFFT_PLAN_ALIGN(plan));
    }
    return ret;
}
