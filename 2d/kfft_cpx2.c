#include "kfft.h"
#include "kfft_trace.h"

// clang-format off
#define kfft_trace_2d(fmt, ...)                                                           \
    kfft_trace("[CPX_2D]"" " fmt, __VA_ARGS__)
// clang-format on

static void
kfft_trace_plan(kfft_comp2_t* P) {
    kfft_trace_2d("%s: %p", "Create KFFT complex plan", (void*)P);
    kfft_trace("\n\t %s - %u", "Total lenght", P->nfft);
    kfft_trace("\n\t %s - %u", "Size X", P->x);
    kfft_trace("\n\t %s - %u", "Size Y", P->y);

    kfft_trace("\n\t %s - %p", "Plan for X dim", (void*)P->plan_x);
    kfft_trace("\n\t %s - %p\n", "Plan for Y dim", (void*)P->plan_y);
}

static inline kfft_return_t
kfft_init(kfft_comp2_t* st) {

    st->plan_x = kfft_config_cpx(st->x, st->flags, 0, st->object.mmgr, NULL);
    if (st->plan_x) {
        if (st->x != st->y) {
            st->plan_y = kfft_config_cpx(st->y, st->flags, 0, st->object.mmgr, NULL);
        } else {
            st->plan_y = st->plan_x;
        }
    }

    return (st->plan_y && st->plan_x) ? KFFT_RET_SUCCESS : KFFT_RET_ALLOC_FAIL;
}

static inline size_t
kfft_calculate(const uint32_t szx, const uint32_t szy, const uint32_t flags) {
    size_t ret = sizeof(kfft_comp2_t);
    size_t delta = 0;

    kfft_config_cpx(szx, flags, 0, NULL, &delta);
    ret += delta;

    if ((szy != szx) || (szy > 1)) {
        kfft_config_cpx(szy, flags, 0, NULL, &delta);
        ret += delta;
    }

    return ret;
}

KFFT_API kfft_comp2_t*
kfft_config2_cpx(const uint32_t x_size, const uint32_t y_size, const uint32_t flags, kfft_pool_t* A,
                 size_t* lenmem) {

    kfft_comp2_t* st = NULL;
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
            st = kfft_internal_alloc(mmgr, sizeof(kfft_comp2_t));
    } else {
        if (A && *lenmem >= memneeded) {
            mmgr = A;

            if (flags & KFFT_FLAG_RENEW)
                kfft_allocator_clear(mmgr);

            st = kfft_internal_alloc(mmgr, sizeof(kfft_comp2_t));
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

static kfft_return_t
kfft_2transform(kfft_comp2_t* st, const kfft_cpx* fin, kfft_cpx* ftmp, kfft_cpx* fout) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    kfft_trace_2d("%s: %p\n", "X-axes transform with plan", (void*)(st->plan_x));
    for (uint32_t i = 0; i < st->y; i++) {
        uint64_t bp = st->x * i;
        ret = kfft_eval_cpx(st->plan_x, &(fin[bp]), &(ftmp[bp]));
        if (ret != KFFT_RET_SUCCESS)
            goto bailout;
    }

    kfft_trace_2d("%s: %p\n", "Transposition matrix plan", (void*)st);
    kfft_math_transpose_cpx(ftmp, fout, st->x, st->y);

    kfft_trace_2d("%s: %p\n", "Y-axes transform with plan", (void*)(st->plan_y));
    for (uint32_t i = 0; i < st->x; i++) {
        uint64_t bp = st->y * i;
        ret = kfft_eval_cpx(st->plan_y, &(fout[bp]), &(ftmp[bp]));
        if (ret != KFFT_RET_SUCCESS)
            goto bailout;
    }
    kfft_trace_2d("%s: %p\n", "Transposition matrix plan", (void*)st);
    kfft_math_transpose_cpx(ftmp, fout, st->y, st->x);

bailout:
    return ret;
}

KFFT_API kfft_return_t
kfft_eval2_cpx(kfft_comp2_t* cfg, const kfft_cpx* fin, kfft_cpx* fout) {
    kfft_return_t ret = KFFT_RET_SUCCESS;
    size_t memneeded = cfg->nfft * sizeof(kfft_cpx);

    kfft_cpx* Ft = KFFT_TMP_ALLOC(memneeded);
    if (Ft) {
        if (fin == fout) {
            kfft_cpx* Fbuf = KFFT_TMP_ALLOC(memneeded);
            if (Fbuf) {
                ret = kfft_2transform(cfg, fin, Ft, Fbuf);
                if (ret == KFFT_RET_SUCCESS) {
                    memcpy(fout, Fbuf, memneeded);
                }
                KFFT_TMP_FREE(Fbuf);
            } else {
                ret = KFFT_RET_BUFFER_FAIL;
            } /* Fbuf */
        } else {
            ret = kfft_2transform(cfg, fin, Ft, fout);
        } /* fin == fout */
        KFFT_TMP_FREE(Ft);
    } else {
        ret = KFFT_RET_BUFFER_FAIL;
    } /* Ft */

    return ret;
}

#undef kfft_trace_2d
