#include "kfft.h"
#include "kfft_trace.h"

// clang-format off
#define kfft_trace_ccnv2(fmt, ...)                                                           \
    kfft_trace("[CCNV_2D]"" " fmt, __VA_ARGS__)
// clang-format on

#if defined(KFFT_TRACE)
static void
kfft_trace_plan(kfft_plan_c2cnv* P) {
    kfft_trace_ccnv2("%s: %p", "Create KFFT 2D convolution complex plan", (void*)P);
    kfft_trace("\n\t %s - %u", "Total lenght", P->nfft);
    kfft_trace("\n\t %s - %u", "Size X", P->x);
    kfft_trace("\n\t %s - %u", "Size Y", P->y);
    kfft_trace("\n\t %s - %p", "Plan for one lenght", (void*)P->plan_fwd);
    kfft_trace("\n\t %s - %p\n", "Plan inverse for one lenght", (void*)P->plan_inv);
}
#endif /*KFFT_TRACE */

static inline kfft_return_t
kfft_init(kfft_plan_c2cnv* plan) {
    KFFT_OMP(omp parallel sections shared(plan)) {
        KFFT_OMP(omp section) {
            plan->plan_fwd = kfft_config2_cpx(plan->x, plan->y, KFFT_CHECK_FLAGS(plan->flags),
                                              KFFT_PLAN_MMGR(plan), NULL);
        }
        KFFT_OMP(omp section) {
            plan->plan_inv = kfft_config2_cpx(plan->x, plan->y,
                                              KFFT_CHECK_FLAGS(plan->flags | KFFT_FLAG_INVERSE),
                                              KFFT_PLAN_MMGR(plan), NULL);
        }
    }
    return ((plan->plan_fwd) && (plan->plan_inv)) ? KFFT_RET_SUCCESS : KFFT_RET_ALLOC_FAIL;
}

static inline size_t
kfft_calculate(const uint32_t x, const uint32_t y, const uint32_t flags) {
    size_t ret = sizeof(kfft_plan_c2cnv);
    size_t delta = 0;
    KFFT_OMP(omp parallel sections shared(ret) private(delta)) {
        KFFT_OMP(omp section) {
            kfft_config2_cpx(x, y, KFFT_CHECK_FLAGS(flags), NULL, &delta);
            ret += delta;
        }
        KFFT_OMP(omp section) {
            kfft_config2_cpx(x, y, KFFT_CHECK_FLAGS(flags | KFFT_FLAG_INVERSE), NULL, &delta);
            ret += delta;
        }
    }

    return ret;
}

KFFT_API kfft_plan_c2cnv*
kfft_config2_conv_cpx(const uint32_t x, const uint32_t y, const uint32_t flags, kfft_pool_t* A,
                      size_t* lenmem) {
    kfft_plan_c2cnv* plan = NULL;
    size_t memneeded = kfft_calculate(x, y, flags);
    KFFT_ALGO_PLAN_PREPARE(plan, flags, kfft_plan_c2cnv, memneeded, A, lenmem);
    if (plan) {

        plan->nfft = x * y;
        plan->x = x * y;
        plan->y = x * y;
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
kfft_eval2_conv_cpx(kfft_plan_c2cnv* plan, const kfft_cpx* fin_A, const kfft_cpx* fin_B,
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
                KFFT_OMP(omp section) { retA = kfft_eval2_cpx(plan->plan_fwd, fin_A, bufA); }
                KFFT_OMP(omp section) { retB = kfft_eval2_cpx(plan->plan_fwd, fin_B, bufB); }
            }
            if ((retA == KFFT_RET_SUCCESS) && (retB == KFFT_RET_SUCCESS)) {
                VEXFUNC(plan, kfft_math_adamar_cpx, bufB, bufA, plan->nfft);
                ret = kfft_eval2_cpx(plan->plan_inv, bufB, fout);
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
