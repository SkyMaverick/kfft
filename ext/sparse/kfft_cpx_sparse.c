#include "kfft.h"
#include "kfft_trace.h"

// clang-format off
#define kfft_trace_spr(fmt, ...)                                                           \
    kfft_trace("[CPX_SPR]"" " fmt, __VA_ARGS__)
// clang-format on

#if defined(KFFT_TRACE)
static void
kfft_trace_plan(kfft_plan_csparse* P) {
    KFFT_OMP(omp critical(trace_log)) {
        kfft_trace_spr("%s: %p", "Create KFFT complex plan", (void*)P);
        kfft_trace_raw("\n\t %s - %u", "Total lenght", P->nfft);
        kfft_trace_raw("\n\t %s - %u", "Dim lenght", P->dnfft);
        kfft_trace_raw("\n\t %s - %u", "Dims", P->dims);
        kfft_trace_raw("\n\t %s - %u", "Steps", P->step);
        kfft_trace_raw("\n\t %s - %p\n", "Plan for one lenght", (void*)P->subst);
    }
}
#endif /*KFFT_TRACE */

static inline kfft_return_t
kfft_init(kfft_plan_csparse* st) {
    st->subst = kfft_config_cpx(st->dnfft, KFFT_CHECK_FLAGS(st->flags), st->object.mmgr, NULL);
    return (st->subst) ? KFFT_RET_SUCCESS : KFFT_RET_ALLOC_FAIL;
}

static inline size_t
kfft_calculate(const uint32_t nfft, const uint32_t flags) {
    size_t ret;
    kfft_config_cpx(nfft, KFFT_CHECK_FLAGS(flags), NULL, &ret);
    ret += sizeof(kfft_plan_csparse);

    return ret;
}

KFFT_API kfft_plan_csparse*
kfft_config_sparse_cpx(const uint32_t nfft, const uint32_t flags, const uint32_t dims,
                       uint32_t step, kfft_pool_t* A, size_t* lenmem) {
    kfft_plan_csparse* st = NULL;

    size_t dim_nfft = (nfft + step) / (dims + step);
    size_t memneeded = kfft_calculate(dim_nfft, flags);

    KFFT_ALGO_PLAN_PREPARE(st, flags, kfft_plan_csparse, memneeded, A, lenmem);
    if (__likely__(st)) {
        st->nfft = nfft;
        st->dnfft = dim_nfft;
        st->dims = dims;
        st->step = step;

        st->flags = flags;

        if (__unlikely__(kfft_init(st) != KFFT_RET_SUCCESS)) {
            KFFT_ALGO_PLAN_TERMINATE(st, A);
            return NULL;
        }
#ifdef KFFT_TRACE
        kfft_trace_plan(st);
#endif
    }
    return st;
}

#if defined(KFFT_MEMLESS_MODE)
static kfft_return_t
kfft_process_memless(kfft_plan_csparse* plan, const kfft_cpx* fin, kfft_cpx* fout) {
    kfft_return_t ret = KFFT_RET_SUCCESS;
    uint32_t memneeded = plan->subst->nfft * sizeof(kfft_cpx);

    kfft_cpx* fbuf = KFFT_TMP_ALLOC(memneeded, KFFT_PLAN_ALIGN(plan));
    if (__likely__(fbuf)) {
        for (uint32_t n = 0; n < plan->dims; n++) {
            // forward scramble input buffer
            for (uint32_t i = 0; i < plan->subst->nfft; i++) {
                C_CPY(fbuf[i], fin[i * (plan->dims + plan->step) + n]);
            }
            ret = kfft_eval_cpx(plan->subst, fbuf, fbuf);
            if (__likely__(ret == KFFT_RET_SUCCESS)) {
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
#else  /* KFFT_MEMLESS_MODE */

static kfft_return_t
kfft_process(kfft_plan_csparse* plan, const kfft_cpx* fin, kfft_cpx* fout) {
    kfft_return_t ret = KFFT_RET_SUCCESS;
    uint32_t memneeded = plan->subst->nfft * sizeof(kfft_cpx);

    KFFT_OMP(omp parallel for schedule(static))
    for (uint32_t n = 0; n < plan->dims; n++) {
        kfft_cpx* fbuf = KFFT_TMP_ALLOC(memneeded, KFFT_PLAN_ALIGN(plan));
        if (__likely__(fbuf)) {
            // forward scramble input buffer
            for (uint32_t i = 0; i < plan->subst->nfft; i++) {
                C_CPY(fbuf[i], fin[i * (plan->dims + plan->step) + n]);
            }
            ret = kfft_eval_cpx(plan->subst, fbuf, fbuf);
            if (__likely__(ret == KFFT_RET_SUCCESS)) {
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
kfft_eval_sparse_cpx(kfft_plan_csparse* cfg, const kfft_cpx* fin, kfft_cpx* fout) {
    if (__unlikely__((cfg->dims < 2) && (!(cfg->step))))
        return kfft_eval_cpx(cfg->subst, fin, fout);

#if defined(KFFT_MEMLESS_MODE)
    return kfft_process_memless(cfg, fin, fout);
#else
    return kfft_process(cfg, fin, fout);
#endif /* KFFT_MEMLESS_MODE */
}

static inline void
shift_internal(kfft_cpx* buf, kfft_cpx* ftmp, const uint32_t nfft, const uint32_t dims,
               uint32_t step, const bool is_inverse, kfft_pool_t* mmgr) {
    for (uint32_t n = 0; n < dims; n++) {
        for (uint32_t i = 0; i < nfft; i++)
            C_CPY(ftmp[i], buf[i * (dims + step) + n]);

        kfft_shift_cpx(ftmp, nfft, is_inverse, mmgr);

        for (uint32_t i = 0; i < nfft; i++)
            C_CPY(buf[i * (dims + step) + n], ftmp[i]);
    }
}

KFFT_API void
kfft_shift_sparse_cpx(kfft_cpx* buf, kfft_cpx* ftmp, const uint32_t nfft, const uint32_t dims,
                      uint32_t step, const bool is_inverse, kfft_pool_t* mmgr) {
    size_t dim_nfft = (nfft + step) / (dims + step);
#if !defined(KFFT_MEMLESS_MODE)
    if (__unlikely__(ftmp == NULL)) {
        kfft_cpx* tbuf = KFFT_TMP_ALLOC(sizeof(kfft_cpx) * dim_nfft, mmgr->align);
        if (__likely__(tbuf)) {
            shift_internal(buf, tbuf, dim_nfft, dims, step, is_inverse, mmgr);
            KFFT_TMP_FREE(tbuf, mmgr->align);
        }
    } else
#endif /* KFFT_MEMLESS_MODE */
        shift_internal(buf, ftmp, dim_nfft, dims, step, is_inverse, mmgr);
}
