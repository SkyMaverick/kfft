/* The guts header contains all the multiplication and addition macros that are defined for
 fixed or floating point complex numbers.  It also delares the kf_ internal functions.
 */
#include "kfft.h"

#include "kfft_math.h"
#include "kfft_trace.h"

static inline kfft_cpx
kfft_sclr_twiddle(uint32_t i, const kfft_plan_sclr* P) {
    kfft_cpx ret;

    kfft_scalar phase = -KFFT_CONST_PI * ((kfft_scalar)(i + 1) / P->basis->nfft + .5);
    if (P->flags & KFFT_FLAG_INVERSE)
        phase *= -1;
    kf_cexp(&ret, phase);
    return ret;
}

#if defined(KFFT_MEMLESS_MODE)
    #define SUPER_TWIDDLE(i, P) kfft_sclr_twiddle(i, P)
#else
    #define SUPER_TWIDDLE(i, P) P->super_twiddles[i]
#endif /* KFFT_MEMLESS_MODE */

#ifdef KFFT_TRACE

static void
kfft_trace_plan(kfft_plan_sclr* P) {
    KFFT_OMP(omp critical(trace_log)) {
        kfft_trace_scalar("%s: %p", "Create KFFT scalar plan", (void*)P);
        kfft_trace_raw("\n\t %s - %u", "nfft", P->nfft);
        kfft_trace_raw("\n\t %s - %u : ", "flags", P->flags);
        kfft_trace_raw("\n\t %s - %p", "Uses complex plan", (void*)(P->basis));
        kfft_trace_raw("\n\t %s - %p\n", "Scalar twiddles", (void*)(P->super_twiddles));
    }
}

#endif /* KFFT_TRACE */

/*******************************************************************************
              Non-symmetric signal (use complex buffer transform)
 ******************************************************************************/

static inline size_t
calculate_trivial(const uint32_t nfft, const uint32_t flags) {
    size_t ret = 0;
    kfft_config_cpx(nfft, flags, NULL, &ret);
    if (ret > 0)
        ret += sizeof(kfft_plan_sclr);
    return ret;
}

KFFT_API kfft_plan_sclr*
config_trivial(const uint32_t nfft, const uint32_t flags, kfft_pool_t* A, size_t* lenmem) {
    kfft_plan_sclr* P = NULL;
    size_t memneeded = calculate_trivial(nfft, flags);

    KFFT_ALGO_PLAN_PREPARE(P, flags, kfft_plan_sclr, memneeded, A, lenmem);

    if (__likely__(P)) {
        P->basis = kfft_config_cpx(nfft, KFFT_CHECK_FLAGS(flags), A, NULL);
        if (__unlikely__(P->basis == NULL)) {
            KFFT_ALGO_PLAN_TERMINATE(P, A);
            return NULL;
        }

        P->nfft = nfft;
        P->flags = flags;

#ifdef KFFT_TRACE
        kfft_trace_plan(P);
#endif
    }
    return P;
}

static kfft_return_t
eval_trivial(kfft_plan_sclr* plan, const kfft_scalar* fin, kfft_cpx* fout, kfft_cpx* ftmp) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    kfft_cpx* fbuf =
        (ftmp) ? ftmp : KFFT_TMP_ALLOC(sizeof(kfft_cpx) * plan->nfft, KFFT_PLAN_ALIGN(plan));
    if (__likely__(fbuf)) {
        KFFT_TMP_ZEROMEM(fbuf, sizeof(kfft_cpx) * plan->nfft);

        for (uint32_t i = 0; i < plan->nfft; i++)
            fbuf[i].r = fin[i]; //< COPY to complex buffer

        ret = kfft_eval_cpx(plan->basis, fbuf, fout);

        if (__unlikely__(ftmp == NULL))
            KFFT_TMP_FREE(fbuf, KFFT_PLAN_ALIGN(plan));
    } else {
        ret = KFFT_RET_BUFFER_FAIL;
    }
    return ret;
}

static kfft_return_t
evali_trivial(kfft_plan_sclr* plan, const kfft_cpx* fin, kfft_scalar* fout, kfft_cpx* ftmp) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    kfft_cpx* fbuf =
        (ftmp) ? ftmp : KFFT_TMP_ALLOC(sizeof(kfft_cpx) * plan->nfft, KFFT_PLAN_ALIGN(plan));
    if (__likely__(fbuf)) {
        ret = kfft_eval_cpx(plan->basis, fin, fbuf);

        for (uint32_t i = 0; i < plan->nfft; i++)
            fout[i] = fbuf[i].r; //< COPY to scalar buffer

        if (__unlikely__(ftmp == NULL))
            KFFT_TMP_FREE(fbuf, KFFT_PLAN_ALIGN(plan));
    } else {
        ret = KFFT_RET_BUFFER_FAIL;
    }
    return ret;
}

/*******************************************************************************
                Symmetric signal (use nyquist frequences N/2)
 ******************************************************************************/
static inline size_t
calculate_nyquist(const uint32_t nfft, const uint32_t flags) {
    size_t ret = 0;
    kfft_config_cpx(nfft / 2, flags, NULL, &ret);
    if (ret > 0) {
#if !defined(KFFT_MEMLESS_MODE)
        ret += nfft * sizeof(kfft_cpx) / 2;
#endif /* not KFFT_MEMLESS_MODE */
        ret += sizeof(kfft_plan_sclr);
    }
    return ret;
}

static kfft_plan_sclr*
config_nyquist(const uint32_t nfft, const uint32_t flags, kfft_pool_t* A, size_t* lenmem) {
    kfft_plan_sclr* P = NULL;
    size_t memneeded = calculate_nyquist(nfft, flags);

    KFFT_ALGO_PLAN_PREPARE(P, flags, kfft_plan_sclr, memneeded, A, lenmem);

    if (__likely__(P)) {

        P->nfft = nfft;
        P->flags = flags;

        uint32_t half_nfft = nfft / 2;

        P->basis = kfft_config_cpx(half_nfft, KFFT_CHECK_FLAGS(flags), P->object.mmgr, NULL);
        if (__unlikely__(P->basis == NULL)) {
            KFFT_ALGO_PLAN_TERMINATE(P, A);
            return NULL;
        }

#if !defined(KFFT_MEMLESS_MODE)
        //        if (__likely__(half_nfft > 1)) {
        P->super_twiddles = kfft_pool_alloc(P->object.mmgr, sizeof(kfft_cpx) * (half_nfft));
        if (__unlikely__(P->super_twiddles == NULL)) {
            KFFT_ALGO_PLAN_TERMINATE(P, A);
            return NULL;
        }
        //        }
        for (uint32_t i = 0; i < half_nfft; ++i) {
            P->super_twiddles[i] = kfft_sclr_twiddle(i, P);
        }
#endif /* not KFFT_MEMLESS_MODE */

#ifdef KFFT_TRACE
        kfft_trace_plan(P);
#endif
    }
    return P;
}

static inline kfft_return_t
rebuild_forward(const kfft_plan_sclr* plan, const kfft_cpx* fin, kfft_cpx* fout) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    uint32_t k, ncfft, nfft;

    nfft = plan->nfft;
    ncfft = plan->basis->nfft;

    kfft_cpx fpnk, fpk, f1k, f2k, tw, tdc;

    tdc.r = fin[0].r;
    tdc.i = fin[0].i;

    fout[0].r = tdc.r + tdc.i;
    fout[ncfft].r = tdc.r - tdc.i;
    fout[ncfft].i = fout[0].i = 0;

    for (k = 1; k <= ncfft / 2; k++) {
        fpk = fin[k];
        fpnk.r = fin[ncfft - k].r;
        fpnk.i = -fin[ncfft - k].i;

        C_ADD(f1k, fpk, fpnk);
        C_SUB(f2k, fpk, fpnk);
        C_MUL(tw, f2k, SUPER_TWIDDLE(k - 1, plan));

        fout[k].r = HALF_OF(f1k.r + tw.r);
        fout[k].i = HALF_OF(f1k.i + tw.i);
        fout[ncfft - k].r = HALF_OF(f1k.r - tw.r);
        fout[ncfft - k].i = HALF_OF(tw.i - f1k.i);
    }

    if (__unlikely__(plan->flags & KFFT_FLAG_EXPAND_SCALAR)) {
        kfft_trace_scalar("Expand buffer with plan: %p\n", (void*)plan);
        for (k = 1; k < nfft; k++) {
            fout[nfft - k].r = fout[k].r;
            fout[nfft - k].i = -fout[k].i;
        }
    }

    return ret;
}

static inline kfft_return_t
rebuild_inverse(const kfft_plan_sclr* plan, const kfft_cpx* fin, kfft_cpx* fout) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    uint32_t k, ncfft = plan->basis->nfft;

    fout[0].r = fin[0].r + fin[ncfft].r;
    fout[0].i = fin[0].r - fin[ncfft].r;

    for (k = 1; k <= ncfft / 2; ++k) {
        kfft_cpx fk, fnkc, fek, fok, tmp;
        fk = fin[k];
        fnkc.r = fin[ncfft - k].r;
        fnkc.i = -fin[ncfft - k].i;

        C_ADD(fek, fk, fnkc);
        C_SUB(tmp, fk, fnkc);
        C_MUL(fok, tmp, SUPER_TWIDDLE(k - 1, plan));
        C_ADD(fout[k], fek, fok);
        C_SUB(fout[ncfft - k], fek, fok);

        fout[ncfft - k].i *= -1;
    }
    return ret;
}

static kfft_return_t
eval_nyquist(kfft_plan_sclr* plan, const kfft_scalar* fin, kfft_cpx* fout, kfft_cpx* ftmp) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    size_t seq_size = plan->basis->nfft + 1;
    kfft_cpx* fbuf =
        (ftmp) ? ftmp : KFFT_TMP_ALLOC(sizeof(kfft_cpx) * seq_size, KFFT_PLAN_ALIGN(plan));
    if (__likely__(fbuf)) {
        KFFT_TMP_ZEROMEM(fbuf, seq_size);

        ret = kfft_eval_cpx(plan->basis, (kfft_cpx*)fin, fbuf);

        if (ret == KFFT_RET_SUCCESS)
            ret = rebuild_forward(plan, fbuf, fout);

        if (__unlikely__(ftmp == NULL))
            KFFT_TMP_FREE(fbuf, KFFT_PLAN_ALIGN(plan));
    } else {
        ret = KFFT_RET_BUFFER_FAIL;
    }

    return ret;
}

static kfft_return_t
evali_nyquist(kfft_plan_sclr* plan, const kfft_cpx* fin, kfft_scalar* fout, kfft_cpx* ftmp) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

    size_t seq_size = plan->basis->nfft + 1;
    kfft_cpx* fbuf =
        (ftmp) ? ftmp : KFFT_TMP_ALLOC(sizeof(kfft_cpx) * seq_size, KFFT_PLAN_ALIGN(plan));
    if (__likely__(fbuf)) {
        KFFT_TMP_ZEROMEM(fbuf, seq_size);

        ret = rebuild_inverse(plan, fin, fbuf);
        if (ret == KFFT_RET_SUCCESS)
            kfft_eval_cpx(plan->basis, fbuf, (kfft_cpx*)fout);

        if (__unlikely__(ftmp == NULL))
            KFFT_TMP_FREE(fbuf, KFFT_PLAN_ALIGN(plan));
    } else {
        ret = KFFT_RET_BUFFER_FAIL;
    }

    return ret;
}

/*******************************************************************************
                            API functionality
 ******************************************************************************/

kfft_return_t
kfft_eval_scalar_internal(kfft_plan_sclr* plan, const kfft_scalar* fin, kfft_cpx* fout,
                          kfft_cpx* ftmp) {

    if (__unlikely__(plan->basis->flags & KFFT_FLAG_INVERSE))
        return KFFT_RET_IMPROPER_PLAN;

    return (plan->nfft % 2) ? eval_trivial(plan, fin, fout, ftmp)
                            : eval_nyquist(plan, fin, fout, ftmp);
}

kfft_return_t
kfft_evali_scalar_internal(kfft_plan_sclr* plan, const kfft_cpx* fin, kfft_scalar* fout,
                           kfft_cpx* ftmp) {

    if (__unlikely__(!(plan->basis->flags & KFFT_FLAG_INVERSE)))
        return KFFT_RET_IMPROPER_PLAN;

    return (plan->nfft % 2) ? evali_trivial(plan, fin, fout, ftmp)
                            : evali_nyquist(plan, fin, fout, ftmp);
}

KFFT_API kfft_plan_sclr*
kfft_config_scalar(const uint32_t nfft, const uint32_t flags, kfft_pool_t* A, size_t* lenmem) {
    return (nfft % 2) ? config_trivial(nfft, flags, A, lenmem)
                      : config_nyquist(nfft, flags, A, lenmem);
}
KFFT_API kfft_return_t
kfft_eval_scalar(kfft_plan_sclr* plan, const kfft_scalar* fin, kfft_cpx* fout) {
    return kfft_eval_scalar_internal(plan, fin, fout, NULL);
}
KFFT_API kfft_return_t
kfft_evali_scalar(kfft_plan_sclr* plan, const kfft_cpx* fin, kfft_scalar* fout) {
    return kfft_evali_scalar_internal(plan, fin, fout, NULL);
}
