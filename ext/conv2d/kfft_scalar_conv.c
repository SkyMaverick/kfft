#include "kfft.h"
#include "kfft_trace.h"

// clang-format off
#define kfft_trace_scnv2(fmt, ...)                                                           \
    kfft_trace("[SCNV_2D]"" " fmt, __VA_ARGS__)
// clang-format on

#if defined(KFFT_TRACE)
static void
kfft_trace_plan(kfft_scnv2_t* P) {
    kfft_trace_scnv2("%s: %p", "Create KFFT convolution scalar plan", (void*)P);
    kfft_trace("\n\t %s - %u", "Total lenght", P->nfft);
    kfft_trace("\n\t %s - %u", "Size X", P->x);
    kfft_trace("\n\t %s - %u", "Size Y", P->y);
    kfft_trace("\n\t %s - %p", "Plan for one lenght", (void*)P->plan_fwd);
    kfft_trace("\n\t %s - %p\n", "Plan inverse for one lenght", (void*)P->plan_inv);
}
#endif /*KFFT_TRACE */
static inline kfft_return_t
kfft_init(kfft_scnv2_t* st) {
    KFFT_OMP(omp parallel sections shared(st)) {
        KFFT_OMP(omp section) {
            st->plan_fwd = kfft_config2_scalar(st->x, st->y, KFFT_CHECK_FLAGS(st->flags),
                                               KFFT_PLAN_MMGR(st), NULL);
        }
        KFFT_OMP(omp section) {
            st->plan_inv =
                kfft_config2_scalar(st->x, st->y, KFFT_CHECK_FLAGS(st->flags | KFFT_FLAG_INVERSE),
                                    KFFT_PLAN_MMGR(st), NULL);
        }
    }
    return ((st->plan_fwd) && (st->plan_inv)) ? KFFT_RET_SUCCESS : KFFT_RET_ALLOC_FAIL;
}

static inline size_t
kfft_calculate(const uint32_t x, const uint32_t y, const uint32_t flags) {
    size_t ret = sizeof(kfft_scnv2_t);
    size_t delta = 0;
    KFFT_OMP(omp parallel sections shared(ret) private(delta)) {
        KFFT_OMP(omp section) {
            kfft_config2_scalar(x, y, KFFT_CHECK_FLAGS(flags), NULL, &delta);
            ret += delta;
        }
        KFFT_OMP(omp section) {
            kfft_config2_scalar(x, y, KFFT_CHECK_FLAGS(flags | KFFT_FLAG_INVERSE), NULL, &delta);
            ret += delta;
        }
    }

    return ret;
}

KFFT_API kfft_scnv2_t*
kfft_config2_conv_scalar(const uint32_t x, const uint32_t y, const uint32_t flags, kfft_pool_t* A,
                         size_t* lenmem) {
    kfft_scnv2_t* st = NULL;
    size_t memneeded = kfft_calculate(x, y, flags);
    KFFT_ALGO_PLAN_PREPARE(st, flags, kfft_scnv2_t, memneeded, A, lenmem);
    if (st) {

        st->nfft = x * y;
        st->x = x;
        st->y = y;
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
kfft_eval2_conv_scalar(kfft_scnv2_t* plan, const kfft_scalar* fin_A, const kfft_scalar* fin_B,
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