#include "kfft.h"
#include "kfft_trace.h"

// clang-format off
#define kfft_trace_spr(fmt, ...)                                                           \
    kfft_trace("[SCR_SPR]"" " fmt, __VA_ARGS__)
// clang-format on

#if defined(KFFT_TRACE)
static void
kfft_trace_plan(kfft_ssparse_t* P) {
    kfft_trace_spr("%s: %p", "Create KFFT scalar plan", (void*)P);
    kfft_trace("\n\t %s - %u", "Total lenght", P->nfft);
    kfft_trace("\n\t %s - %u", "Dims", P->dims);
    kfft_trace("\n\t %s - %u", "Steps", P->step);

    kfft_trace("\n\t %s - %p\n", "Plan for one lenght", (void*)P->subst);
}
#endif /*KFFT_TRACE */

static inline kfft_return_t
kfft_init(kfft_ssparse_t* st) {
    st->subst = kfft_config_scalar(st->nfft, KFFT_CHECK_FLAGS(st->flags), st->object.mmgr, NULL);
    return (st->subst) ? KFFT_RET_SUCCESS : KFFT_RET_ALLOC_FAIL;
}

static inline size_t
kfft_calculate(const uint32_t nfft, const uint32_t flags) {
    size_t ret;
    kfft_config_scalar(nfft, KFFT_CHECK_FLAGS(flags), NULL, &ret);
    ret += sizeof(kfft_ssparse_t);

    return ret;
}

KFFT_API kfft_ssparse_t*
kfft_config_sparse_scalar(const uint32_t nfft, const uint32_t flags, const uint32_t dims,
                          uint32_t step, kfft_pool_t* A, size_t* lenmem) {
    kfft_ssparse_t* st = NULL;

    size_t dim_nfft = (nfft + step) / (dims + step);
    size_t memneeded = kfft_calculate(dim_nfft, flags);

    KFFT_ALGO_PLAN_PREPARE(st, flags, kfft_ssparse_t, memneeded, A, lenmem);
    if (st) {
        st->nfft = dim_nfft;
        st->dims = dims;
        st->step = step;

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

#if defined(KFFT_MEMLESS_MODE)
static inline kfft_return_t
kfft_eval_process_memless(kfft_ssparse_t* plan, const kfft_scalar* fin, kfft_cpx* fout) {
    kfft_return_t ret = KFFT_RET_SUCCESS;
    uint32_t memneeded = plan->nfft * (sizeof(kfft_scalar) + sizeof(kfft_cpx));

    kfft_cpx* fbuf_cpx = KFFT_TMP_ALLOC(memneeded, KFFT_PLAN_ALIGN(plan));
    if (fbuf_cpx) {
        for (uint32_t n = 0; n < plan->dims; n++) {
            kfft_scalar* fbuf_scr = (kfft_scalar*)((kfft_cpx*)fbuf_cpx + plan->nfft);
            // forward scramble input buffer
            for (uint32_t i = 0; i < plan->nfft; i++) {
                fbuf_scr[i] = fin[i * (plan->dims + plan->step) + n];
            }
            ret = kfft_eval_scalar(plan->subst, fbuf_scr, fbuf_cpx);
            if (ret == KFFT_RET_SUCCESS) {
                // backward put values in output buffer
                for (uint32_t i = 0; i < plan->nfft; i++) {
                    C_CPY(fout[i * (plan->dims + plan->step) + n], fbuf_cpx[i]);
                }
            }
            KFFT_TMP_FREE(fbuf_cpx, KFFT_PLAN_ALIGN(plan));
        }
    } else {
        ret = KFFT_RET_BUFFER_FAIL;
    }

    return ret;
}
#endif /* KFFT_MEMLESS_MODE */

static inline kfft_return_t
kfft_eval_process(kfft_ssparse_t* plan, const kfft_scalar* fin, kfft_cpx* fout) {
    kfft_return_t ret = KFFT_RET_SUCCESS;
    uint32_t memneeded = plan->nfft * (sizeof(kfft_scalar) + sizeof(kfft_cpx));

    KFFT_OMP( omp parallel for schedule(static))
    for (uint32_t n = 0; n < plan->dims; n++) {
        kfft_cpx* fbuf_cpx = KFFT_TMP_ALLOC(memneeded, KFFT_PLAN_ALIGN(plan));
        if (fbuf_cpx) {
            kfft_scalar* fbuf_scr = (kfft_scalar*)((kfft_cpx*)fbuf_cpx + plan->nfft);
            // forward scramble input buffer
            for (uint32_t i = 0; i < plan->nfft; i++) {
                fbuf_scr[i] = fin[i * (plan->dims + plan->step) + n];
            }
            ret = kfft_eval_scalar(plan->subst, fbuf_scr, fbuf_cpx);
            if (ret == KFFT_RET_SUCCESS) {
                // backward put values in output buffer
                for (uint32_t i = 0; i < plan->nfft; i++) {
                    C_CPY(fout[i * (plan->dims + plan->step) + n], fbuf_cpx[i]);
                }
            }
            KFFT_TMP_FREE(fbuf_cpx, KFFT_PLAN_ALIGN(plan));
        } else {
            ret = KFFT_RET_BUFFER_FAIL;
        }
    }
    return ret;
}

KFFT_API kfft_return_t
kfft_eval_sparse_scalar(kfft_ssparse_t* plan, const kfft_scalar* fin, kfft_cpx* fout) {
    if ((plan->dims < 2) && (plan->step))
        return kfft_eval_scalar(plan->subst, fin, fout);

#if defined(KFFT_MEMLESS_MODE)
    return kfft_eval_process_memless(plan, fin, fout);
#else
    return kfft_eval_process(plan, fin, fout);
#endif /* KFFT_MEMLESS_MODE */
}

#if defined(KFFT_MEMLESS_MODE)
KFFT_API kfft_return_t
kfft_evali_process_memless(kfft_ssparse_t* plan, const kfft_cpx* fin, kfft_scalar* fout) {
    kfft_return_t ret = KFFT_RET_SUCCESS;
    uint32_t memneeded = plan->nfft * (sizeof(kfft_scalar) + sizeof(kfft_cpx));

    kfft_cpx* fbuf_cpx = KFFT_TMP_ALLOC(memneeded, KFFT_PLAN_ALIGN(plan));
    if (fbuf_cpx) {
        for (uint32_t n = 0; n < plan->dims; n++) {
            kfft_scalar* fbuf_scr = (kfft_scalar*)((kfft_cpx*)fbuf_cpx + plan->nfft);
            // forward scramble input buffer
            for (uint32_t i = 0; i < plan->nfft; i++) {
                C_CPY(fbuf_cpx[i], fin[i * (plan->dims + plan->step) + n]);
            }
            ret = kfft_evali_scalar(plan->subst, fbuf_cpx, fbuf_scr);
            if (ret == KFFT_RET_SUCCESS) {
                // backward put values in output buffer
                for (uint32_t i = 0; i < plan->nfft; i++) {
                    fout[i * (plan->dims + plan->step) + n] = fbuf_scr[i];
                }
            }
            KFFT_TMP_FREE(fbuf_cpx, KFFT_PLAN_ALIGN(plan));
        }
    } else {
        ret = KFFT_RET_BUFFER_FAIL;
    }
    return ret;
}
#endif /* KFFT_MEMLESS_MODE */

KFFT_API kfft_return_t
kfft_evali_process(kfft_ssparse_t* plan, const kfft_cpx* fin, kfft_scalar* fout) {
    kfft_return_t ret = KFFT_RET_SUCCESS;
    uint32_t memneeded = plan->nfft * (sizeof(kfft_scalar) + sizeof(kfft_cpx));

    KFFT_OMP( omp parallel for schedule(static))
    for (uint32_t n = 0; n < plan->dims; n++) {
        kfft_cpx* fbuf_cpx = KFFT_TMP_ALLOC(memneeded, KFFT_PLAN_ALIGN(plan));
        if (fbuf_cpx) {
            kfft_scalar* fbuf_scr = (kfft_scalar*)((kfft_cpx*)fbuf_cpx + plan->nfft);
            // forward scramble input buffer
            for (uint32_t i = 0; i < plan->nfft; i++) {
                C_CPY(fbuf_cpx[i], fin[i * (plan->dims + plan->step) + n]);
            }
            ret = kfft_evali_scalar(plan->subst, fbuf_cpx, fbuf_scr);
            if (ret == KFFT_RET_SUCCESS) {
                // backward put values in output buffer
                for (uint32_t i = 0; i < plan->nfft; i++) {
                    fout[i * (plan->dims + plan->step) + n] = fbuf_scr[i];
                }
            }
            KFFT_TMP_FREE(fbuf_cpx, KFFT_PLAN_ALIGN(plan));
        } else {
            ret = KFFT_RET_BUFFER_FAIL;
        }
    }
    return ret;
}

KFFT_API kfft_return_t
kfft_evali_sparse_scalar(kfft_ssparse_t* plan, const kfft_cpx* fin, kfft_scalar* fout) {
    if ((plan->dims < 2) && (plan->step))
        return kfft_evali_scalar(plan->subst, fin, fout);

#if defined(KFFT_MEMLESS_MODE)
    return kfft_evali_process_memless(plan, fin, fout);
#else
    return kfft_evali_process(plan, fin, fout);
#endif /* KFFT_MEMLESS_MODE */
}

static inline void
shift_internal(kfft_scalar* buf, kfft_scalar* ftmp, const uint32_t nfft, const uint32_t dims,
               uint32_t step, const bool is_inverse, kfft_pool_t* mmgr) {
    for (uint32_t n = 0; n < dims; n++) {
        for (uint32_t i = 0; i < nfft; i++)
            ftmp[i] = buf[i * (dims + step) + n];

        kfft_shift_scalar(ftmp, nfft, is_inverse, mmgr);

        for (uint32_t i = 0; i < nfft; i++)
            buf[i * (dims + step) + n] = ftmp[i];
    }
}

KFFT_API void
kfft_shift_sparse_scalar(kfft_scalar* buf, kfft_scalar* ftmp, const uint32_t nfft,
                         const uint32_t dims, uint32_t step, const bool is_inverse,
                         kfft_pool_t* mmgr) {
    size_t dim_nfft = (nfft + step) / (dims + step);
#if !defined(KFFT_MEMLESS_MODE)
    if (ftmp == NULL) {
        kfft_scalar* tbuf = KFFT_TMP_ALLOC(sizeof(kfft_scalar) * dim_nfft, mmgr->align);
        if (tbuf) {
            shift_internal(buf, tbuf, dim_nfft, dims, step, is_inverse, mmgr);
            KFFT_TMP_FREE(tbuf, mmgr->align);
        }
    } else
#endif /* KFFT_MEMLESS_MODE */
        shift_internal(buf, ftmp, dim_nfft, dims, step, is_inverse, mmgr);
    return;
}
