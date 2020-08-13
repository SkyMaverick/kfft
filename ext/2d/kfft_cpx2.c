#include "kfft.h"
#include "kfft_trace.h"

// clang-format off
#define kfft_trace_2d(fmt, ...)                                                           \
    kfft_trace("[CPX_2D]"" " fmt, __VA_ARGS__)
// clang-format on

#if defined(KFFT_TRACE)
static void
kfft_trace_plan(kfft_plan_c2d* P) {
    kfft_trace_2d("%s: %p", "Create KFFT complex plan", (void*)P);
    kfft_trace("\n\t %s - %u", "Total lenght", P->nfft);
    kfft_trace("\n\t %s - %u", "Size X", P->x);
    kfft_trace("\n\t %s - %u", "Size Y", P->y);

    kfft_trace("\n\t %s - %p", "Plan for X dim", (void*)P->plan_x);
    kfft_trace("\n\t %s - %p\n", "Plan for Y dim", (void*)P->plan_y);
}
#endif /*KFFT_TRACE */

static inline kfft_return_t
kfft_init(kfft_plan_c2d* st) {
    if (st->x != st->y) {
        KFFT_OMP(omp parallel sections shared(st)) {
            KFFT_OMP(omp section) {
                st->plan_x =
                    kfft_config_cpx(st->x, KFFT_CHECK_FLAGS(st->flags), st->object.mmgr, NULL);
            }
            KFFT_OMP(omp section) {
                st->plan_y =
                    kfft_config_cpx(st->y, KFFT_CHECK_FLAGS(st->flags), st->object.mmgr, NULL);
            }
        }
    } else {
        st->plan_y = st->plan_x =
            kfft_config_cpx(st->x, KFFT_CHECK_FLAGS(st->flags), st->object.mmgr, NULL);
    }
    return (st->plan_y && st->plan_x) ? KFFT_RET_SUCCESS : KFFT_RET_ALLOC_FAIL;
}

static inline size_t
kfft_calculate(const uint32_t szx, const uint32_t szy, const uint32_t flags) {
    size_t ret = sizeof(kfft_plan_c2d);
    size_t r1,r2;
    r1 = r2 = 0;

    if ((szy != szx) || (szy > 1)) {
        KFFT_OMP(omp parallel sections) {
            KFFT_OMP(omp section) {
                kfft_config_cpx(szx, KFFT_CHECK_FLAGS(flags), NULL, &r1);
            }
            KFFT_OMP(omp section) {
                kfft_config_cpx(szy, KFFT_CHECK_FLAGS(flags), NULL, &r2);
            }
        }
        ret += r1 + r2;
    } else {
        KFFT_UNUSED_VAR(r2);
        kfft_config_cpx(szx, KFFT_CHECK_FLAGS(flags), NULL, &r1);
        ret += r1;
    }
    return ret;
}

KFFT_API kfft_plan_c2d*
kfft_config2_cpx(const uint32_t x_size, const uint32_t y_size, const uint32_t flags, kfft_pool_t* A,
                 size_t* lenmem) {
    kfft_plan_c2d* st = NULL;
    size_t memneeded = kfft_calculate(x_size, y_size, flags);

    KFFT_ALGO_PLAN_PREPARE(st, flags, kfft_plan_c2d, memneeded, A, lenmem);
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
kfft_2transform_normal(kfft_plan_c2d* st, const kfft_cpx* fin, kfft_cpx* fout) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    kfft_cpx* ftmp = KFFT_TMP_ALLOC(st->nfft * sizeof(kfft_cpx), KFFT_PLAN_ALIGN(st));
    if (ftmp) {
        kfft_trace_2d("%s: %p\n", "X-axes transform with plan", (void*)(st->plan_x));

    KFFT_OMP( omp parallel for schedule(static))
    for (uint32_t i = 0; i < st->y; i++) {
        uint64_t bp = st->x * i;
        ret = kfft_eval_cpx(st->plan_x, &(fin[bp]), &(ftmp[bp]));
    }

    kfft_trace_2d("%s: %p\n", "Transposition matrix plan", (void*)st);
    kfft_math_transpose_cpx(ftmp, fout, st->x, st->y);

    kfft_trace_2d("%s: %p\n", "Y-axes transform with plan", (void*)(st->plan_y));

    KFFT_OMP( omp parallel for schedule(static))
    for (uint32_t i = 0; i < st->x; i++) {
        uint64_t bp = st->y * i;
        ret = kfft_eval_cpx(st->plan_y, &(fout[bp]), &(ftmp[bp]));
    }
    kfft_trace_2d("%s: %p\n", "Transposition matrix plan", (void*)st);
    kfft_math_transpose_cpx(ftmp, fout, st->y, st->x);

    KFFT_TMP_FREE(ftmp, KFFT_PLAN_ALIGN(st));
    } else {
        ret = KFFT_RET_BUFFER_FAIL;
    }
    return ret;
}

#if defined(KFFT_MEMLESS_MODE)
static inline kfft_return_t
kfft_2transform_memless(kfft_plan_c2d* st, kfft_cpx* fin) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    kfft_trace_2d("%s: %p\n", "X-axes transform with plan", (void*)(st->plan_x));

    KFFT_OMP( omp parallel for schedule(static))
    for (uint32_t i = 0; i < st->y; i++) {
        uint64_t bp = st->x * i;
        ret = kfft_eval_cpx(st->plan_x, &(fin[bp]), &(fin[bp]));
    }

    kfft_trace_2d("%s: %p\n", "Transposition matrix plan (in-place)", (void*)st);
    kfft_math_transpose_ip_cpx(fin, st->x, st->y);

    kfft_trace_2d("%s: %p\n", "Y-axes transform with plan", (void*)(st->plan_y));

    KFFT_OMP( omp parallel for schedule(static))
    for (uint32_t i = 0; i < st->x; i++) {
        uint64_t bp = st->y * i;
        ret = kfft_eval_cpx(st->plan_y, &(fin[bp]), &(fin[bp]));
    }
    kfft_trace_2d("%s: %p\n", "Transposition matrix plan (in-place)", (void*)st);
    kfft_math_transpose_ip_cpx(fin, st->y, st->x);

    return ret;
}
#endif /* KFFT_MEMLESS_MODE */

KFFT_API kfft_return_t
kfft_eval2_cpx(kfft_plan_c2d* cfg, const kfft_cpx* fin, kfft_cpx* fout) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    size_t memneeded = cfg->nfft * sizeof(kfft_cpx);
#if defined(KFFT_MEMLESS_MODE)
    if (cfg->flags & KFFT_FLAG_GENERIC_ONLY) {
        if (fin != fout)
            memcpy(fout, fin, memneeded);
        ret = kfft_2transform_memless(cfg, fout);
    } else {
#endif /* KFFT_MEMLESS_MODE */
        if (fin == fout) {
            kfft_cpx* Fbuf = KFFT_TMP_ALLOC(memneeded, KFFT_PLAN_ALIGN(cfg));
            if (Fbuf) {
                ret = kfft_2transform_normal(cfg, fin, Fbuf);
                if (ret == KFFT_RET_SUCCESS) {
                    memcpy(fout, Fbuf, memneeded);
                }
                KFFT_TMP_FREE(Fbuf, KFFT_PLAN_ALIGN(cfg));
            } else {
                ret = KFFT_RET_BUFFER_FAIL;
            } /* Fbuf */
        } else {
            ret = kfft_2transform_normal(cfg, fin, fout);
        } /* fin == fout */
#if defined(KFFT_MEMLESS_MODE)
    }  /* KFFT_FLAG_GENERIC_ONLY */
#endif /* memless mode */
    return ret;
}

static void
shift_internal(kfft_cpx* buf, kfft_cpx* ftmp, const uint32_t sz_x, const uint32_t sz_y,
               const bool is_inverse, kfft_pool_t* mmgr) {
    kfft_trace_2d("%s\n", "X-axes shift transform");

    KFFT_OMP( omp parallel for schedule(static))
    for (uint32_t i = 0; i < sz_y; i++) {
        uint64_t bp = sz_x * i;
        kfft_shift_cpx(&(buf[bp]), sz_x, is_inverse, mmgr);
    }
    if (ftmp != NULL) {
        kfft_trace_2d("%s\n", "Transposition matrix");
        kfft_math_transpose_cpx(buf, ftmp, sz_x, sz_y);

        kfft_trace_2d("%s\n", "Y-axes shift transform");

        KFFT_OMP( omp parallel for schedule(static))
        for (uint32_t i = 0; i < sz_x; i++) {
            uint64_t bp = sz_y * i;
            kfft_shift_cpx(&(ftmp[bp]), sz_y, is_inverse, mmgr);
        }
        kfft_trace_2d("%s\n", "Transposition matrix");
        kfft_math_transpose_cpx(ftmp, buf, sz_y, sz_x);
    } else { /* ftmp != NULL */
        kfft_trace_2d("%s\n", "Transposition matrix (in-place)");
        kfft_math_transpose_ip_cpx(buf, sz_x, sz_y);

        kfft_trace_2d("%s\n", "Y-axes shift transform");

        KFFT_OMP( omp parallel for schedule(static))
        for (uint32_t i = 0; i < sz_x; i++) {
            uint64_t bp = sz_y * i;
            kfft_shift_cpx(&(buf[bp]), sz_y, is_inverse, mmgr);
        }
        kfft_trace_2d("%s\n", "Transposition matrix (in-place)");
        kfft_math_transpose_ip_cpx(buf, sz_y, sz_x);
    } /* ftmp != NULL */
}

KFFT_API void
kfft_shift2_cpx(kfft_cpx* buf, kfft_cpx* ftmp, const uint32_t sz_x, const uint32_t sz_y,
                const bool is_inverse, kfft_pool_t* mmgr) {
#if !defined(KFFT_MEMLESS_MODE)
    if (ftmp == NULL) {
        kfft_cpx* tbuf = KFFT_TMP_ALLOC(sizeof(kfft_cpx) * sz_x * sz_y, mmgr->align);
        if (tbuf) {
            shift_internal(buf, tbuf, sz_x, sz_y, is_inverse, mmgr);
            KFFT_TMP_FREE(tbuf, mmgr->align);
        }
    } else
#endif /* KFFT_MEMLESS_MODE */
        shift_internal(buf, ftmp, sz_x, sz_y, is_inverse, mmgr);
}

#undef kfft_trace_2d
