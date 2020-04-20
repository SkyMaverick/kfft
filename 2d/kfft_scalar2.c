#include "kfft.h"
#include "kfft_trace.h"

#if defined(KFFT_USE_OPENMP)
    #include <omp.h>
#endif

// clang-format off
#define kfft_trace_2d(fmt, ...)                                                           \
    kfft_trace("[SCLR_2D]"" " fmt, __VA_ARGS__)
// clang-format on

#define CPXSP(X) ((X)->substate)

#if defined(KFFT_TRACE)
static void
kfft_trace_plan(kfft_sclr2_t* P) {
    kfft_trace_2d("%s: %p", "Create KFFT complex plan", (void*)P);
    kfft_trace("\n\t %s - %u", "Total lenght", P->nfft);
    kfft_trace("\n\t %s - %u", "Size X", P->x);
    kfft_trace("\n\t %s - %u", "Size Y", P->y);

    kfft_trace("\n\t %s - %p", "Plan for X dim", (void*)P->plan_x);
    kfft_trace("\n\t %s - %p\n", "Plan for Y dim", (void*)P->plan_y);
}
#endif /*KFFT_TRACE */

static inline kfft_return_t
kfft_init(kfft_sclr2_t* st) {

    st->plan_x = kfft_config_scalar(st->x, st->flags, st->object.mmgr, NULL);
    if (st->plan_x) {
        if (st->x != st->y) {
            st->plan_y = kfft_config_scalar(st->y, st->flags, st->object.mmgr, NULL);
        } else {
            st->plan_y = st->plan_x;
        }
    }

    return (st->plan_y && st->plan_x) ? KFFT_RET_SUCCESS : KFFT_RET_ALLOC_FAIL;
}

static inline size_t
kfft_calculate(const uint32_t szx, const uint32_t szy, const uint32_t flags) {
    size_t ret = sizeof(kfft_sclr2_t);
    size_t delta = 0;

    kfft_config_scalar(szx, flags, NULL, &delta);
    ret += delta;

    if ((szy != szx) || (szy > 1)) {
        kfft_config_scalar(szy, flags, NULL, &delta);
        ret += delta;
    }

    return ret;
}

KFFT_API kfft_sclr2_t*
kfft_config2_scalar(const uint32_t x_size, const uint32_t y_size, const uint32_t flags,
                    kfft_pool_t* A, size_t* lenmem) {

    kfft_sclr2_t* st = NULL;
    size_t memneeded = kfft_calculate(x_size, y_size, flags);

    kfft_pool_t* mmgr = NULL;
    bool flag_create = false;

    if (lenmem == NULL) {
        if (A == NULL) {
            mmgr = kfft_allocator_create(memneeded);
            flag_create = true;

            kfft_trace_2d("%s: %p\n", "Create new allocator and plan", (void*)mmgr);
        } else {
            mmgr = A;
            kfft_trace_2d("%s: %p\n", "Use allocator and create plan", (void*)mmgr);
        }
        if (mmgr)
            st = kfft_internal_alloc(mmgr, sizeof(kfft_sclr2_t));
    } else {
        if (A && *lenmem >= memneeded) {
            mmgr = A;

            if (flags & KFFT_FLAG_RENEW)
                kfft_allocator_clear(mmgr);

            st = kfft_internal_alloc(mmgr, sizeof(kfft_sclr2_t));
            kfft_trace_2d("%s: %p\n", "Reuse allocator and create plan", (void*)mmgr);
        }
        *lenmem = memneeded;
    }

    if (!st) {
    bailout:
        if (mmgr && (flag_create == true))
            kfft_allocator_free(mmgr);
        return NULL;
    }

    st->object.mmgr = mmgr;

    st->nfft = x_size * y_size;
    st->x = x_size;
    st->y = y_size;
    st->flags = flags;

    if (kfft_init(st) != KFFT_RET_SUCCESS)
        goto bailout;

#ifdef KFFT_TRACE
    kfft_trace_plan(st);
#endif

    return st;
}

static inline kfft_return_t
kfft_2transform(kfft_sclr2_t* st, const kfft_scalar* fin, kfft_cpx* fout) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    kfft_cpx *fbuf, *ftmp;
    ftmp = KFFT_TMP_ALLOC(2 * st->nfft * sizeof(kfft_cpx));
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
#pragma omp for schedule(static)
        for (uint32_t i = 0; i < st->x; i++) {
            uint64_t bp = st->y * i;
            ret = kfft_eval_cpx(CPXSP(st->plan_y), &(fout[bp]), &(ftmp[bp]));
        }
        kfft_trace_2d("%s: %p\n", "Transposition matrix plan", (void*)st);
        kfft_math_transpose_cpx(ftmp, fout, st->y, st->x);

        KFFT_TMP_FREE(ftmp);
    } else {
        ret = KFFT_RET_BUFFER_FAIL;
    }
    return ret;
}

KFFT_API kfft_return_t
kfft_eval2_scalar(kfft_sclr2_t* cfg, const kfft_scalar* fin, kfft_cpx* fout) {
    return (cfg->flags & KFFT_FLAG_INVERSE) ? KFFT_RET_IMPROPER_PLAN
                                            : kfft_2transform(cfg, fin, fout);
}

static inline kfft_return_t
kfft_2transform_inverse(kfft_sclr2_t* st, const kfft_cpx* fin, kfft_scalar* fout) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    kfft_cpx *ftps, *fbuf, *ftmp;
    ftmp = KFFT_TMP_ALLOC(3 * st->nfft * sizeof(kfft_cpx));

    if (ftmp) {
        ftps = ftmp + st->nfft;
        fbuf = ftps + st->nfft;

        kfft_trace_2d("%s: %p\n", "X-axes transform with plan", (void*)(st->plan_x));
#pragma omp for schedule(static)
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
        KFFT_TMP_FREE(ftmp);
    } else {
        ret = KFFT_RET_BUFFER_FAIL;
    }
    return ret;
}

KFFT_API kfft_return_t
kfft_evali2_scalar(kfft_sclr2_t* cfg, const kfft_cpx* fin, kfft_scalar* fout) {
    return (cfg->flags & KFFT_FLAG_INVERSE) ? kfft_2transform_inverse(cfg, fin, fout)
                                            : KFFT_RET_IMPROPER_PLAN;
}

void
shift_internal(kfft_scalar* buf, kfft_scalar* ftmp, const uint32_t sz_x, const uint32_t sz_y,
               const bool is_inverse) {
    kfft_trace_2d("%s\n", "X-axes shift transform");
#pragma omp for schedule(static)
    for (uint32_t i = 0; i < sz_y; i++) {
        uint64_t bp = sz_x * i;
        kfft_shift_scalar(&(buf[bp]), sz_x, is_inverse);
    }

    kfft_trace_2d("%s\n", "Transposition matrix");
    kfft_math_transpose_scalar(buf, ftmp, sz_x, sz_y);

    kfft_trace_2d("%s\n", "Y-axes shift transform");
#pragma omp for schedule(static)
    for (uint32_t i = 0; i < sz_x; i++) {
        uint64_t bp = sz_y * i;
        kfft_shift_scalar(&(ftmp[bp]), sz_y, is_inverse);
    }
    kfft_trace_2d("%s\n", "Transposition matrix");
    kfft_math_transpose_scalar(ftmp, buf, sz_y, sz_x);
}

KFFT_API void
kfft_shift2_scalar(kfft_scalar* buf, kfft_scalar* ftmp, const uint32_t sz_x, const uint32_t sz_y,
                   const bool is_inverse) {
    if (ftmp == NULL) {
        kfft_scalar* tbuf = KFFT_TMP_ALLOC(sizeof(kfft_cpx) * sz_x * sz_y);
        if (tbuf) {
            shift_internal(buf, tbuf, sz_x, sz_y, is_inverse);
            KFFT_TMP_FREE(tbuf);
        }
    } else {
        shift_internal(buf, ftmp, sz_x, sz_y, is_inverse);
    }
}
#undef kfft_trace_2d