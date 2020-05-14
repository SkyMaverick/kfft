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

    kfft_trace("\n\t %s - %p\n", "Plan for one lenght", (void*)P->subst);
}
#endif /*KFFT_TRACE */

static inline kfft_return_t
kfft_init(kfft_csparse_t* st) {
    st->subst = kfft_config_cpx(st->nfft, KFFT_CHECK_FLAGS(st->flags), st->object.mmgr, NULL);
    return (st->subst) ? KFFT_RET_SUCCESS : KFFT_RET_ALLOC_FAIL;
}

static inline size_t
kfft_calculate(const uint32_t nfft, const uint32_t flags) {
    size_t ret;
    kfft_config_cpx(nfft, KFFT_CHECK_FLAGS(flags), NULL, &ret);
    ret += sizeof(kfft_csparse_t);

    return ret;
}

KFFT_API kfft_csparse_t*
kfft_config_sparse_cpx(const uint32_t nfft, const uint32_t flags, const uint32_t dims,
                       uint32_t step, kfft_pool_t* A, size_t* lenmem) {
    kfft_csparse_t* st = NULL;

    size_t dim_nfft = (nfft + step) / (dims + step);
    size_t memneeded = kfft_calculate(dim_nfft, flags);

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

    st->nfft = dim_nfft;
    st->dims = dims;
    st->step = step;

    st->flags = flags;

    if (kfft_init(st) != KFFT_RET_SUCCESS)
        goto bailout;

#ifdef KFFT_TRACE
    kfft_trace_plan(st);
#endif
    return st;
}

#if defined(KFFT_MEMLESS_MODE)
static kfft_return_t
kfft_process_memless(kfft_csparse_t* plan, const kfft_cpx* fin, kfft_cpx* fout) {
    kfft_return_t ret = KFFT_RET_SUCCESS;
    uint32_t memneeded = plan->subst->nfft * sizeof(kfft_cpx);

    kfft_cpx* fbuf = KFFT_TMP_ALLOC(memneeded, KFFT_PLAN_ALIGN(plan));
    if (fbuf) {
        for (uint32_t n = 0; n < plan->dims; n++) {
            // forward scramble input buffer
            for (uint32_t i = 0; i < plan->subst->nfft; i++) {
                C_CPY(fbuf[i], fin[i * (plan->dims + plan->step) + n]);
            }
            ret = kfft_eval_cpx(plan->subst, fbuf, fbuf);
            if (ret == KFFT_RET_SUCCESS) {
                // backward put values in output buffer
                for (uint32_t i = 0; i < plan->subst->nfft; i++) {
                    C_CPY(fout[i * (plan->dims + plan->step) + n], fbuf[i]);
                }
            }
        }
        KFFT_TMP_FREE(fbuf, KFFT_PLAN_ALIGN(plan));
    } else {
        ret = KFFT_RET_BUFFER_FAIL;
    }
    return ret;
}
#else /* KFFT_MEMLESS_MODE */

static kfft_return_t
kfft_process(kfft_csparse_t* plan, const kfft_cpx* fin, kfft_cpx* fout) {

    kfft_return_t ret = KFFT_RET_SUCCESS;
    uint32_t memneeded = plan->subst->nfft * sizeof(kfft_cpx);

    #if !defined(KFFT_OS_WINDOWS)
        #pragma omp parallel for schedule(static)
    #endif
    for (uint32_t n = 0; n < plan->dims; n++) {
        kfft_cpx* fbuf = KFFT_TMP_ALLOC(memneeded, KFFT_PLAN_ALIGN(plan));
        if (fbuf) {
            // forward scramble input buffer
            for (uint32_t i = 0; i < plan->subst->nfft; i++) {
                C_CPY(fbuf[i], fin[i * (plan->dims + plan->step) + n]);
            }
            ret = kfft_eval_cpx(plan->subst, fbuf, fbuf);
            if (ret == KFFT_RET_SUCCESS) {
                // backward put values in output buffer
                for (uint32_t i = 0; i < plan->subst->nfft; i++) {
                    C_CPY(fout[i * (plan->dims + plan->step) + n], fbuf[i]);
                }
            }
            KFFT_TMP_FREE(fbuf, KFFT_PLAN_ALIGN(plan));
        } else {
            ret = KFFT_RET_BUFFER_FAIL;
        }
    }
    return ret;
}
#endif /* KFFT_MEMLESS_MODE */

KFFT_API kfft_return_t
kfft_eval_sparse_cpx(kfft_csparse_t* cfg, const kfft_cpx* fin, kfft_cpx* fout) {

    if ((cfg->dims < 2) && (cfg->step))
        return kfft_eval_cpx(cfg->subst, fin, fout);

#if defined(KFFT_MEMLESS_MODE)
    return kfft_process_memless(cfg, fin, fout);
#else
    return kfft_process(cfg, fin, fout);
#endif /* KFFT_MEMLESS_MODE */
}
