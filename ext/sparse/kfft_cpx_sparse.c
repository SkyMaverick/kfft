#include "kfft.h"
#include "kfft_trace.h"

#if defined(KFFT_USE_OPENMP)
    #include <omp.h>
#endif

// clang-format off
#define kfft_trace_spr(fmt, ...)                                                           \
    kfft_trace("[CPX_SPR]"" " fmt, __VA_ARGS__)
// clang-format on

#if defined(KFFT_TRACE)
static void
kfft_trace_plan(kfft_csparse_t* P) {
    kfft_trace_spr("%s: %p", "Create KFFT complex plan", (void*)P);
    kfft_trace("\n\t %s - %u", "Total lenght", P->nfft);
    kfft_trace("\n\t %s - %u", "Dims", P->dims);
    kfft_trace("\n\t %s - %u", "Steps", P->step);

    kfft_trace("\n\t %s - %p", "Plan for one lenght", (void*)P->subst);
}
#endif /*KFFT_TRACE */

static inline kfft_return_t
kfft_init(kfft_csparse_t* st) {
    st->subst = kfft_config_cpx(st->nfft, KFFT_CHECK_FLAGS(st->flags), st->object.mmgr, NULL);
    if (st->subst) {
        st->nfft *= st->dims + st->step - 1;
        return KFFT_RET_SUCCESS;
    }
    return KFFT_RET_ALLOC_FAIL;
}

static inline size_t
kfft_calculate(const uint32_t nfft, const uint32_t dims, uint32_t step, const uint32_t flags) {
    size_t ret;
    kfft_config_cpx(nfft, KFFT_CHECK_FLAGS(flags), NULL, &ret);

    ret += sizeof(kfft_csparse_t);

    return ret;
}

KFFT_API kfft_csparse_t*
kfft_config_sparse_cpx(const uint32_t nfft, const uint32_t flags, const uint32_t dims,
                       uint32_t step, kfft_pool_t* A, size_t* lenmem) {
    kfft_csparse_t* st = NULL;
    size_t memneeded = kfft_calculate(nfft, dims, step, flags);

    kfft_pool_t* mmgr = NULL;
    bool flag_create = false;

    if (lenmem == NULL) {
        if (A == NULL) {
            mmgr = kfft_allocator_create(memneeded);
            flag_create = true;

            kfft_trace_spr("%s: %p\n", "Create new allocator and plan", (void*)mmgr);
        } else {
            mmgr = A;
            kfft_trace_spr("%s: %p\n", "Use allocator and create plan", (void*)mmgr);
        }
        if (mmgr)
            st = kfft_internal_alloc(mmgr, sizeof(kfft_csparse_t));
    } else {
        if (A && *lenmem >= memneeded) {
            mmgr = A;

            if (flags & KFFT_FLAG_RENEW)
                kfft_allocator_clear(mmgr);

            st = kfft_internal_alloc(mmgr, sizeof(kfft_csparse_t));
            kfft_trace_spr("%s: %p\n", "Reuse allocator and create plan", (void*)mmgr);
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

    st->nfft = nfft; // !!! Temporary
    st->flags = flags;
    st->dims = dims;
    st->step = step;

    if (kfft_init(st) != KFFT_RET_SUCCESS)
        goto bailout;

#ifdef KFFT_TRACE
    kfft_trace_plan(st);
#endif
    return st;
}

KFFT_API kfft_return_t
kfft_eval_sparse_cpx(kfft_csparse_t* cfg, const kfft_cpx* fin, kfft_cpx* fout) {
    return KFFT_RET_SUCCESS;
}
