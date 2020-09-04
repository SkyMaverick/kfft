#include "kfft.h"
#include "kfft_trace.h"

// clang-format off
#define kfft_trace_2d(fmt, ...)                                                           \
    kfft_trace("[SCLR_2D]"" " fmt, __VA_ARGS__)
// clang-format on

#define CPXSP(X) ((X)->basis)

#if defined(KFFT_TRACE)
static void
kfft_trace_plan(kfft_plan_s2d* P) {
    KFFT_OMP(omp critical(trace_log)) {
        kfft_trace_2d("%s: %p", "Create KFFT complex plan", (void*)P);
        kfft_trace_raw("\n\t %s - %u", "Total lenght", P->nfft);
        kfft_trace_raw("\n\t %s - %u", "Size X", P->x);
        kfft_trace_raw("\n\t %s - %u", "Size Y", P->y);
        kfft_trace_raw("\n\t %s - %p", "Plan for X dim", (void*)P->plan_x);
        kfft_trace_raw("\n\t %s - %p\n", "Plan for Y dim", (void*)P->plan_y);
    }
}
#endif /*KFFT_TRACE */

static inline kfft_return_t
kfft_init(kfft_plan_s2d* plan) {
    if (__likely__(plan->x != plan->y)) {
        KFFT_OMP(omp parallel sections shared(plan)) {
            KFFT_OMP(omp section) {
                plan->plan_x = kfft_config_scalar(plan->x, KFFT_CHECK_FLAGS(plan->flags),
                                                  plan->object.mmgr, NULL);
            }
            KFFT_OMP(omp section) {
                plan->plan_y = kfft_config_scalar(plan->y, KFFT_CHECK_FLAGS(plan->flags),
                                                  plan->object.mmgr, NULL);
            }
        }
    } else {
        plan->plan_y = plan->plan_x =
            kfft_config_scalar(plan->x, KFFT_CHECK_FLAGS(plan->flags), plan->object.mmgr, NULL);
    }
    return (plan->plan_y && plan->plan_x) ? KFFT_RET_SUCCESS : KFFT_RET_ALLOC_FAIL;
}

static inline size_t
kfft_calculate(const uint32_t szx, const uint32_t szy, const uint32_t flags) {
    size_t r1, r2, ret = sizeof(kfft_plan_s2d);
    r1 = r2 = 0;

    if (__likely__(szy > 1)) {
        if (__unlikely__(szx == szy)) {
            kfft_config_scalar(szx, KFFT_CHECK_FLAGS(flags), NULL, &r1);
            ret += r1;
        } else {
            KFFT_OMP(omp parallel sections) {
                KFFT_OMP(omp section) {
                    kfft_config_scalar(szx, KFFT_CHECK_FLAGS(flags), NULL, &r1);
                }
                KFFT_OMP(omp section) {
                    kfft_config_scalar(szy, KFFT_CHECK_FLAGS(flags), NULL, &r2);
                }
            }
            ret += r1 + r2;
        }
    } else {
        KFFT_UNUSED_VAR(r2);
        kfft_config_cpx(szx, KFFT_CHECK_FLAGS(flags), NULL, &r1);
        ret += r1;
    }
    return ret;
}

KFFT_API kfft_plan_s2d*
kfft_config2_scalar(const uint32_t x_size, const uint32_t y_size, const uint32_t flags,
                    kfft_pool_t* A, size_t* lenmem) {
    kfft_plan_s2d* plan = NULL;
    size_t memneeded = kfft_calculate(x_size, y_size, flags);

    KFFT_ALGO_PLAN_PREPARE(plan, flags, kfft_plan_s2d, memneeded, A, lenmem);
    if (__likely__(plan)) {
        plan->nfft = x_size * y_size;
        plan->x = x_size;
        plan->y = y_size;
        plan->flags = flags;

        if (__unlikely__(kfft_init(plan) != KFFT_RET_SUCCESS)) {
            KFFT_ALGO_PLAN_TERMINATE(plan, A);
            return NULL;
        }
#ifdef KFFT_TRACE
        kfft_trace_plan(plan);
#endif
    }
    return plan;
}

static inline kfft_return_t
kfft_2transform(kfft_plan_s2d* plan, const kfft_scalar* fin, kfft_cpx* fout) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    kfft_cpx *fbuf, *ftmp;
    size_t memneed = plan->nfft + plan->plan_x->nfft;
    ftmp = KFFT_TMP_ALLOC(memneed * sizeof(kfft_cpx), KFFT_PLAN_ALIGN(plan));
    if (__likely__(ftmp)) {
        fbuf = ftmp + plan->nfft;
        kfft_trace_2d("%s: %p\n", "X-axes transform with plan", (void*)(plan->plan_x));
        for (uint32_t i = 0; i < plan->y; i++) {
            uint64_t bp = plan->x * i;
            ret = kfft_eval_scalar_internal(plan->plan_x, &(fin[bp]), &(ftmp[bp]), fbuf);
        }

        kfft_trace_2d("%s: %p\n", "Transposition matrix plan", (void*)plan);
        kfft_math_transpose_cpx(ftmp, fout, plan->x, plan->y);

        kfft_trace_2d("%s: %p\n", "Y-axes transform with plan", (void*)(plan->plan_y));

        KFFT_OMP(omp parallel for schedule(static))
        for (uint32_t i = 0; i < plan->x; i++) {
            uint64_t bp = plan->y * i;
            ret = kfft_eval_cpx(CPXSP(plan->plan_y), &(fout[bp]), &(ftmp[bp]));
        }
        kfft_trace_2d("%s: %p\n", "Transposition matrix plan", (void*)plan);
        kfft_math_transpose_cpx(ftmp, fout, plan->y, plan->x);

        KFFT_TMP_FREE(ftmp, KFFT_PLAN_ALIGN(plan));
    } else {
        ret = KFFT_RET_BUFFER_FAIL;
    }
    return ret;
}

KFFT_API kfft_return_t
kfft_eval2_scalar(kfft_plan_s2d* plan, const kfft_scalar* fin, kfft_cpx* fout) {
    return (plan->flags & KFFT_FLAG_INVERSE) ? KFFT_RET_IMPROPER_PLAN
                                             : kfft_2transform(plan, fin, fout);
}

#if defined(KFFT_MEMLESS_MODE)
static inline kfft_return_t
kfft_2transform_inverse_memless(kfft_plan_s2d* plan, const kfft_cpx* fin, kfft_scalar* fout) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    kfft_cpx *fbuf, *ftmp;

    /* Memory area: */
    size_t memneed = (plan->nfft + 2 * plan->plan_y->nfft) * sizeof(kfft_cpx);

    ftmp = KFFT_TMP_ALLOC(memneed, KFFT_PLAN_ALIGN(plan));

    if (__likely__(ftmp)) {
        fbuf = ftmp + plan->nfft;

        kfft_trace_2d("%s: %p\n", "X-axes transform with plan", (void*)(plan->plan_x));

        KFFT_OMP( omp parallel for schedule(static))
        for (uint32_t i = 0; i < plan->y; i++) {
            uint64_t bp = plan->x * i;
            ret = kfft_eval_cpx(CPXSP(plan->plan_x), &(fin[bp]), &(ftmp[bp]));
        }

        kfft_trace_2d("%s: %p\n", "Transposition matrix plan (in-place)", (void*)plan);
        kfft_math_transpose_ip_cpx(ftmp, plan->x, plan->y);

        kfft_trace_2d("%s: %p\n", "Y-axes transform with plan", (void*)(plan->plan_y));

        for (uint32_t i = 0; i < plan->x; i++) {
            uint64_t bp = plan->y * i;
            ret = kfft_evali_scalar_internal(plan->plan_y, &(ftmp[bp]), (&(fout[bp])), fbuf);
        }
        kfft_trace_2d("%s: %p\n", "Transposition matrix plan (in-place)", (void*)plan);
        kfft_math_transpose_ip_scalar(fout, plan->y, plan->x);
        KFFT_TMP_FREE(ftmp, KFFT_PLAN_ALIGN(plan));
    } else {
        ret = KFFT_RET_BUFFER_FAIL;
    }
    return ret;
}

#else  /* KFFT_MEMLESS_MODE) */

static inline kfft_return_t
kfft_2transform_inverse_normal(kfft_plan_s2d* plan, const kfft_cpx* fin, kfft_scalar* fout) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    kfft_cpx *ftps, *fbuf, *ftmp;
    size_t memneed = 2 * (plan->nfft + plan->plan_y->nfft) * sizeof(kfft_cpx);
    ftmp = KFFT_TMP_ALLOC(memneed, KFFT_PLAN_ALIGN(plan));

    if (__likely__(ftmp)) {
        ftps = ftmp + plan->nfft;
        fbuf = ftps + plan->nfft;

        kfft_trace_2d("%s: %p\n", "X-axes transform with plan", (void*)(plan->plan_x));

        KFFT_OMP( omp parallel for schedule(static))
        for (uint32_t i = 0; i < plan->y; i++) {
            uint64_t bp = plan->x * i;
            ret = kfft_eval_cpx(CPXSP(plan->plan_x), &(fin[bp]), &(ftmp[bp]));
        }

        kfft_trace_2d("%s: %p\n", "Transposition matrix plan", (void*)plan);
        kfft_math_transpose_cpx(ftmp, ftps, plan->x, plan->y);

        kfft_trace_2d("%s: %p\n", "Y-axes transform with plan", (void*)(plan->plan_y));

        for (uint32_t i = 0; i < plan->x; i++) {
            uint64_t bp = plan->y * i;
            ret = kfft_evali_scalar_internal(plan->plan_y, &(ftps[bp]),
                                             (&(((kfft_scalar*)ftmp)[bp])), fbuf);
        }
        kfft_trace_2d("%s: %p\n", "Transposition matrix plan", (void*)plan);
        kfft_math_transpose_scalar((kfft_scalar*)ftmp, fout, plan->y, plan->x);
        KFFT_TMP_FREE(ftmp, KFFT_PLAN_ALIGN(plan));
    } else {
        ret = KFFT_RET_BUFFER_FAIL;
    }
    return ret;
}
#endif /* KFFT_MEMLESS_MODE) */

static inline kfft_return_t
kfft_2transform_inverse(kfft_plan_s2d* plan, const kfft_cpx* fin, kfft_scalar* fout) {
#if defined(KFFT_MEMLESS_MODE)
    return kfft_2transform_inverse_memless(plan, fin, fout);
#else
    return kfft_2transform_inverse_normal(plan, fin, fout);
#endif /* KFFT_MEMLESS_MODE */
}

KFFT_API kfft_return_t
kfft_evali2_scalar(kfft_plan_s2d* plan, const kfft_cpx* fin, kfft_scalar* fout) {
    return (plan->flags & KFFT_FLAG_INVERSE) ? kfft_2transform_inverse(plan, fin, fout)
                                             : KFFT_RET_IMPROPER_PLAN;
}

void
shift_internal(kfft_scalar* buf, kfft_scalar* ftmp, const uint32_t sz_x, const uint32_t sz_y,
               const bool is_inverse, kfft_pool_t* mmgr) {
    kfft_trace_2d("%s\n", "X-axes shift transform");

    KFFT_OMP(omp parallel for schedule(static))
    for (uint32_t i = 0; i < sz_y; i++) {
        uint64_t bp = sz_x * i;
        kfft_shift_scalar(&(buf[bp]), sz_x, is_inverse, mmgr);
    }
    if (__likely__(ftmp != NULL)) {
        kfft_trace_2d("%s\n", "Transposition matrix");
        kfft_math_transpose_scalar(buf, ftmp, sz_x, sz_y);

        kfft_trace_2d("%s\n", "Y-axes shift transform");

        KFFT_OMP( omp parallel for schedule(static))
        for (uint32_t i = 0; i < sz_x; i++) {
            uint64_t bp = sz_y * i;
            kfft_shift_scalar(&(ftmp[bp]), sz_y, is_inverse, mmgr);
        }
        kfft_trace_2d("%s\n", "Transposition matrix");
        kfft_math_transpose_scalar(ftmp, buf, sz_y, sz_x);
    } else { /* ftmp != NULL */
        kfft_trace_2d("%s\n", "Transposition matrix (in-place)");
        kfft_math_transpose_ip_scalar(buf, sz_x, sz_y);

        kfft_trace_2d("%s\n", "Y-axes shift transform");

        KFFT_OMP( omp parallel for schedule(static))
        for (uint32_t i = 0; i < sz_x; i++) {
            uint64_t bp = sz_y * i;
            kfft_shift_scalar(&(buf[bp]), sz_y, is_inverse, mmgr);
        }
        kfft_trace_2d("%s\n", "Transposition matrix (in-place)");
        kfft_math_transpose_ip_scalar(buf, sz_y, sz_x);
    } /* ftmp != NULL */
}

KFFT_API void
kfft_shift2_scalar(kfft_scalar* buf, kfft_scalar* ftmp, const uint32_t sz_x, const uint32_t sz_y,
                   const bool is_inverse, kfft_pool_t* mmgr) {
#if !defined(KFFT_MEMLESS_MODE)
    if (__unlikely__(ftmp == NULL)) {
        kfft_scalar* tbuf = KFFT_TMP_ALLOC(sizeof(kfft_scalar) * sz_x * sz_y, mmgr->align);
        if (__likely__(tbuf)) {
            shift_internal(buf, tbuf, sz_x, sz_y, is_inverse, mmgr);
            KFFT_TMP_FREE(tbuf, mmgr->align);
        }
    } else
#endif /* KFFT_MEMLESS_MODE */
        shift_internal(buf, ftmp, sz_x, sz_y, is_inverse, mmgr);
}
#undef kfft_trace_2d
