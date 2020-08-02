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
    kfft_trace_2d("%s: %p", "Create KFFT complex plan", (void*)P);
    kfft_trace("\n\t %s - %u", "Total lenght", P->nfft);
    kfft_trace("\n\t %s - %u", "Size X", P->x);
    kfft_trace("\n\t %s - %u", "Size Y", P->y);

    kfft_trace("\n\t %s - %p", "Plan for X dim", (void*)P->plan_x);
    kfft_trace("\n\t %s - %p\n", "Plan for Y dim", (void*)P->plan_y);
}
#endif /*KFFT_TRACE */

static inline kfft_return_t
kfft_init(kfft_plan_s2d* st) {
    if (st->x != st->y) {
        KFFT_OMP(omp parallel sections shared(st)) {
            KFFT_OMP(omp section) {
                st->plan_x =
                    kfft_config_scalar(st->x, KFFT_CHECK_FLAGS(st->flags), st->object.mmgr, NULL);
            }
            KFFT_OMP(omp section) {
                st->plan_y =
                    kfft_config_scalar(st->y, KFFT_CHECK_FLAGS(st->flags), st->object.mmgr, NULL);
            }
        }
    } else {
        st->plan_y = st->plan_x =
            kfft_config_scalar(st->x, KFFT_CHECK_FLAGS(st->flags), st->object.mmgr, NULL);
    }
    return (st->plan_y && st->plan_x) ? KFFT_RET_SUCCESS : KFFT_RET_ALLOC_FAIL;
}

static inline size_t
kfft_calculate(const uint32_t szx, const uint32_t szy, const uint32_t flags) {
    size_t ret = sizeof(kfft_plan_s2d);
    size_t delta = 0;

    if ((szy != szx) || (szy > 1)) {
        KFFT_OMP(omp parallel sections shared(ret) private(delta)) {
            KFFT_OMP(omp section) {
                kfft_config_scalar(szx, KFFT_CHECK_FLAGS(flags), NULL, &delta);
                ret += delta;
            }
            KFFT_OMP(omp section) {
                kfft_config_scalar(szy, KFFT_CHECK_FLAGS(flags), NULL, &delta);
                ret += delta;
            }
        }
    } else {
        kfft_config_cpx(szx, KFFT_CHECK_FLAGS(flags), NULL, &delta);
        ret += delta;
    }
    return ret;
}

KFFT_API kfft_plan_s2d*
kfft_config2_scalar(const uint32_t x_size, const uint32_t y_size, const uint32_t flags,
                    kfft_pool_t* A, size_t* lenmem) {
    kfft_plan_s2d* st = NULL;
    size_t memneeded = kfft_calculate(x_size, y_size, flags);

    KFFT_ALGO_PLAN_PREPARE(st, flags, kfft_plan_s2d, memneeded, A, lenmem);
    if (st) {
        st->nfft = x_size * y_size;
        st->x = x_size;
        st->y = y_size;
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

static inline kfft_return_t
kfft_2transform(kfft_plan_s2d* st, const kfft_scalar* fin, kfft_cpx* fout) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    kfft_cpx *fbuf, *ftmp;
    ftmp = KFFT_TMP_ALLOC(2 * st->nfft * sizeof(kfft_cpx), KFFT_PLAN_ALIGN(st));
    if (ftmp) {
        fbuf = ftmp + st->nfft;
        kfft_trace_2d("%s: %p\n", "X-axes transform with plan", (void*)(st->plan_x));
        for (uint32_t i = 0; i < st->y; i++) {
            uint64_t bp = st->x * i;
            ret = kfft_eval_scalar_internal(st->plan_x, &(fin[bp]), &(ftmp[bp]), fbuf);
        }

        kfft_trace_2d("%s: %p\n", "Transposition matrix plan", (void*)st);
        kfft_math_transpose_cpx(ftmp, fout, st->x, st->y);

        kfft_trace_2d("%s: %p\n", "Y-axes transform with plan", (void*)(st->plan_y));

        KFFT_OMP(omp parallel for schedule(static))
        for (uint32_t i = 0; i < st->x; i++) {
            uint64_t bp = st->y * i;
            ret = kfft_eval_cpx(CPXSP(st->plan_y), &(fout[bp]), &(ftmp[bp]));
        }
        kfft_trace_2d("%s: %p\n", "Transposition matrix plan", (void*)st);
        kfft_math_transpose_cpx(ftmp, fout, st->y, st->x);

        KFFT_TMP_FREE(ftmp, KFFT_PLAN_ALIGN(st));
    } else {
        ret = KFFT_RET_BUFFER_FAIL;
    }
    return ret;
}

KFFT_API kfft_return_t
kfft_eval2_scalar(kfft_plan_s2d* cfg, const kfft_scalar* fin, kfft_cpx* fout) {
    return (cfg->flags & KFFT_FLAG_INVERSE) ? KFFT_RET_IMPROPER_PLAN
                                            : kfft_2transform(cfg, fin, fout);
}

#if defined(KFFT_MEMLESS_MODE)
static inline kfft_return_t
kfft_2transform_inverse_memless(kfft_plan_s2d* st, const kfft_cpx* fin, kfft_scalar* fout) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    kfft_cpx *fbuf, *ftmp;
    ftmp = KFFT_TMP_ALLOC(2 * st->nfft * sizeof(kfft_cpx), KFFT_PLAN_ALIGN(st));

    if (ftmp) {
        fbuf = ftmp + st->nfft;

        kfft_trace_2d("%s: %p\n", "X-axes transform with plan", (void*)(st->plan_x));

        KFFT_OMP( omp parallel for schedule(static))
        for (uint32_t i = 0; i < st->y; i++) {
            uint64_t bp = st->x * i;
            ret = kfft_eval_cpx(CPXSP(st->plan_x), &(fin[bp]), &(ftmp[bp]));
        }

        kfft_trace_2d("%s: %p\n", "Transposition matrix plan (in-place)", (void*)st);
        kfft_math_transpose_ip_cpx(ftmp, st->x, st->y);

        kfft_trace_2d("%s: %p\n", "Y-axes transform with plan", (void*)(st->plan_y));

        for (uint32_t i = 0; i < st->x; i++) {
            uint64_t bp = st->y * i;
            ret = kfft_evali_scalar_internal(st->plan_y, &(ftmp[bp]), (&(fout[bp])), fbuf);
        }
        kfft_trace_2d("%s: %p\n", "Transposition matrix plan (in-place)", (void*)st);
        kfft_math_transpose_ip_scalar(fout, st->y, st->x);
        KFFT_TMP_FREE(ftmp, KFFT_PLAN_ALIGN(st));
    } else {
        ret = KFFT_RET_BUFFER_FAIL;
    }
    return ret;
}

#else  /* KFFT_MEMLESS_MODE) */

static inline kfft_return_t
kfft_2transform_inverse_normal(kfft_plan_s2d* st, const kfft_cpx* fin, kfft_scalar* fout) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    kfft_cpx *ftps, *fbuf, *ftmp;
    ftmp = KFFT_TMP_ALLOC(3 * st->nfft * sizeof(kfft_cpx), KFFT_PLAN_ALIGN(st));

    if (ftmp) {
        ftps = ftmp + st->nfft;
        fbuf = ftps + st->nfft;

        kfft_trace_2d("%s: %p\n", "X-axes transform with plan", (void*)(st->plan_x));

        KFFT_OMP( omp parallel for schedule(static))
        for (uint32_t i = 0; i < st->y; i++) {
            uint64_t bp = st->x * i;
            ret = kfft_eval_cpx(CPXSP(st->plan_x), &(fin[bp]), &(ftmp[bp]));
        }

        kfft_trace_2d("%s: %p\n", "Transposition matrix plan", (void*)st);
        kfft_math_transpose_cpx(ftmp, ftps, st->x, st->y);

        kfft_trace_2d("%s: %p\n", "Y-axes transform with plan", (void*)(st->plan_y));

        for (uint32_t i = 0; i < st->x; i++) {
            uint64_t bp = st->y * i;
            ret = kfft_evali_scalar_internal(st->plan_y, &(ftps[bp]), (&(((kfft_scalar*)ftmp)[bp])),
                                             fbuf);
        }
        kfft_trace_2d("%s: %p\n", "Transposition matrix plan", (void*)st);
        kfft_math_transpose_scalar((kfft_scalar*)ftmp, fout, st->y, st->x);
        KFFT_TMP_FREE(ftmp, KFFT_PLAN_ALIGN(st));
    } else {
        ret = KFFT_RET_BUFFER_FAIL;
    }
    return ret;
}
#endif /* KFFT_MEMLESS_MODE) */

static inline kfft_return_t
kfft_2transform_inverse(kfft_plan_s2d* st, const kfft_cpx* fin, kfft_scalar* fout) {
#if defined(KFFT_MEMLESS_MODE)
    return kfft_2transform_inverse_memless(st, fin, fout);
#else
    return kfft_2transform_inverse_normal(st, fin, fout);
#endif /* KFFT_MEMLESS_MODE */
}

KFFT_API kfft_return_t
kfft_evali2_scalar(kfft_plan_s2d* cfg, const kfft_cpx* fin, kfft_scalar* fout) {
    return (cfg->flags & KFFT_FLAG_INVERSE) ? kfft_2transform_inverse(cfg, fin, fout)
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
    if (ftmp != NULL) {
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
    if (ftmp == NULL) {
        kfft_scalar* tbuf = KFFT_TMP_ALLOC(sizeof(kfft_scalar) * sz_x * sz_y, mmgr->align);
        if (tbuf) {
            shift_internal(buf, tbuf, sz_x, sz_y, is_inverse, mmgr);
            KFFT_TMP_FREE(tbuf, mmgr->align);
        }
    } else
#endif /* KFFT_MEMLESS_MODE */
        shift_internal(buf, ftmp, sz_x, sz_y, is_inverse, mmgr);
}

static inline kfft_return_t
kfft_2transform_norm(kfft_plan_s2d* st, const kfft_scalar* fin, kfft_scalar* fout) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    kfft_cpx *fbuf, *ftmp;
    ftmp = KFFT_TMP_ALLOC(2 * st->nfft * sizeof(kfft_cpx), KFFT_PLAN_ALIGN(st));
    if (ftmp) {
        fbuf = ftmp + st->nfft;
        kfft_trace_2d("%s: %p\n", "X-axes transform with plan", (void*)(st->plan_x));
        for (uint32_t i = 0; i < st->y; i++) {
            uint64_t bp = st->x * i;
            ret = kfft_eval_scalar_internal(st->plan_x, &(fin[bp]), &(ftmp[bp]), fbuf);
        }

        kfft_trace_2d("%s: %p\n", "Transposition matrix plan", (void*)st);
        kfft_math_transpose_cpx(ftmp, fbuf, st->x, st->y);

        kfft_trace_2d("%s: %p\n", "Y-axes transform with plan", (void*)(st->plan_y));

        KFFT_OMP(omp parallel for schedule(static))
        for (uint32_t i = 0; i < st->x; i++) {
            uint64_t bp = st->y * i;
            ret = kfft_eval_cpx(CPXSP(st->plan_y), &(fbuf[bp]), &(ftmp[bp]));
        }
        kfft_trace_2d("%s: %p\n", "Transposition matrix plan", (void*)st);
        kfft_math_transpose_cpx(ftmp, fbuf, st->y, st->x);

        kfft_math_magnitude(fbuf, fout, st->nfft);

        KFFT_TMP_FREE(ftmp, KFFT_PLAN_ALIGN(st));
    } else {
        ret = KFFT_RET_BUFFER_FAIL;
    }
    return ret;
}

KFFT_API kfft_return_t
kfft_eval2_scalar_norm(kfft_plan_s2d* cfg, const kfft_scalar* fin, kfft_scalar* fout) {
    return (cfg->flags & KFFT_FLAG_INVERSE) ? KFFT_RET_IMPROPER_PLAN
                                            : kfft_2transform_norm(cfg, fin, fout);
}

#undef kfft_trace_2d
