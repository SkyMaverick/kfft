/* The guts header contains all the multiplication and addition macros that are defined for
 fixed or floating point complex numbers.  It also delares the kf_ internal functions.
 */
#include "kfft.h"

#include "kfft_math.h"
#include "kfft_trace.h"

static inline kfft_cpx
kfft_plan_sclrwiddle(uint32_t i, const kfft_plan_sclr* P) {
    kfft_cpx ret;

    kfft_scalar phase = -KFFT_CONST_PI * ((kfft_scalar)(i + 1) / P->nfft + .5);
    if (P->flags & KFFT_FLAG_INVERSE)
        phase *= -1;
    kf_cexp(&ret, phase);
    return ret;
}

#if defined(KFFT_MEMLESS_MODE)
    #define SUPER_TWIDDLE(i, P) kfft_plan_sclrwiddle(i, P)
#else
    #define SUPER_TWIDDLE(i, P) P->super_twiddles[i]
#endif /* KFFT_MEMLESS_MODE */

static inline size_t
kfft_calculate(const uint32_t nfft, const uint32_t flags) {
    size_t ret = sizeof(kfft_plan_sclr);
#if !defined(KFFT_MEMLESS_MODE)
    ret += sizeof(kfft_cpx) * (nfft / 2);
#endif
    size_t subsize = 0;

    kfft_config_cpx(nfft, flags, NULL, &subsize);
    return ret + subsize;
}

#ifdef KFFT_TRACE

static void
kfft_trace_plan(kfft_plan_sclr* P) {
    kfft_trace_scalar("%s: %p", "Create KFFT scalar plan", (void*)P);
    kfft_trace("\n\t %s - %u", "nfft", P->nfft);
    kfft_trace("\n\t %s - %u : ", "flags", P->flags);
    kfft_trace("\n\t %s - %p", "Uses complex plan", (void*)(P->substate));
    kfft_trace("\n\t %s - %p\n", "scalar twiddles", (void*)(P->super_twiddles));
}

#endif /* KFFT_TRACE */
/* ********************************************************************************
      TODO  Functionality
******************************************************************************** */

KFFT_API kfft_plan_sclr*
kfft_config_scalar(const uint32_t nfft, const uint32_t flags, kfft_pool_t* A, size_t* lenmem) {
    kfft_plan_sclr* st = NULL;
    size_t memneeded = kfft_calculate(nfft, flags);

    KFFT_ALGO_PLAN_PREPARE(st, flags, kfft_plan_sclr, memneeded, A, lenmem);

    if (st) {
        st->substate = kfft_config_cpx(nfft, KFFT_CHECK_FLAGS(flags), st->object.mmgr, NULL);
        if (st->substate == NULL) {
            KFFT_ALGO_PLAN_TERMINATE(st, A);
            return NULL;
        }

        st->nfft = st->substate->nfft;
        st->flags = st->substate->flags;

#if !defined(KFFT_MEMLESS_MODE)
        if (nfft > 1) {
            st->super_twiddles = kfft_pool_alloc(st->object.mmgr, sizeof(kfft_cpx) * (nfft / 2));
            if (st->super_twiddles == NULL) {
                KFFT_ALGO_PLAN_TERMINATE(st, A);
                return NULL;
            }
        }
        for (uint32_t i = 0; i < nfft / 2; ++i) {
            st->super_twiddles[i] = kfft_plan_sclrwiddle(i, st);
        }
#endif /* not KFFT_MEMLESS_MODE */

#ifdef KFFT_TRACE
        kfft_trace_plan(st);
#endif
    }
    return st;
}

static inline kfft_return_t
eval_forward_internal(const kfft_plan_sclr* st, const kfft_cpx* Fin, kfft_cpx* Fout) {
    kfft_return_t ret = KFFT_RET_SUCCESS;
    kfft_cpx fpnk, fpk, f1k, f2k, tw;

    uint32_t k, ncfft = st->nfft;

    C_CPY(Fout[0], Fin[0]);

    for (k = 1; k <= ncfft / 2; ++k) {
        fpk = Fin[k];
        fpnk.r = Fin[ncfft - k].r;
        fpnk.i = -Fin[ncfft - k].i;

        C_ADD(f1k, fpk, fpnk);
        C_SUB(f2k, fpk, fpnk);
        C_MUL(tw, f2k, SUPER_TWIDDLE(k - 1, st) /* st->super_twiddles[k - 1] */);

        Fout[k].r = HALF_OF(f1k.r + tw.r);
        Fout[k].i = HALF_OF(f1k.i + tw.i);
        Fout[ncfft - k].r = HALF_OF(f1k.r - tw.r);
        Fout[ncfft - k].i = HALF_OF(tw.i - f1k.i);
    }

    return ret;
}

static kfft_return_t
eval_func(kfft_plan_sclr* stu, kfft_cpx* tmpbuf, const kfft_scalar* timedata, kfft_cpx* freqdata) {
    kfft_return_t ret = KFFT_RET_SUCCESS;
    uint32_t ncfft = stu->nfft;

    KFFT_TMP_ZEROMEM(tmpbuf, ncfft * sizeof(kfft_cpx));
    for (uint32_t i = 0; i < ncfft; i++) {
        freqdata[i].r = timedata[i];
        freqdata[i].i = 0;
    }
    ret = kfft_eval_cpx(stu->substate, freqdata, tmpbuf);
    if (ret == KFFT_RET_SUCCESS) {
        ret = eval_forward_internal(stu, tmpbuf, freqdata);
    }

    return ret;
}

kfft_return_t
kfft_eval_scalar_internal(kfft_plan_sclr* stu, const kfft_scalar* timedata, kfft_cpx* freqdata,
                          kfft_cpx* tmpbuf) {
    kfft_return_t ret = KFFT_RET_SUCCESS;
    /* input buffer timedata is stored row-wise */
    kfft_plan_sclr* st = (kfft_plan_sclr*)stu;

    if (st->substate->flags & KFFT_FLAG_INVERSE) {
        return KFFT_RET_IMPROPER_PLAN;
    }

    uint32_t ncfft = stu->nfft;
    if (tmpbuf == NULL) {
        kfft_cpx* tbuf = KFFT_TMP_ALLOC(sizeof(kfft_cpx) * ncfft, KFFT_PLAN_ALIGN(stu));
        if (tbuf) {
            ret = eval_func(stu, tbuf, timedata, freqdata);
            KFFT_TMP_FREE(tbuf, KFFT_PLAN_ALIGN(stu));
        } else {
            ret = KFFT_RET_BUFFER_FAIL;
        }
    } else {
        ret = eval_func(stu, tmpbuf, timedata, freqdata);
    }

    return ret;
}

KFFT_API kfft_return_t
kfft_eval_scalar(kfft_plan_sclr* stu, const kfft_scalar* timedata, kfft_cpx* freqdata) {
    return kfft_eval_scalar_internal(stu, timedata, freqdata, NULL);
}

static inline kfft_return_t
eval_inverse_internal(const kfft_plan_sclr* st, const kfft_cpx* Fin, kfft_cpx* Fout) {
    kfft_return_t ret = KFFT_RET_SUCCESS;
    kfft_cpx fk, fnkc, fek, fok, tmp;

    uint32_t k, ncfft = st->nfft;

    C_CPY(Fout[0], Fin[0]);
    C_MULBYSCALAR(Fout[0], 2);

    for (k = 1; k <= ncfft / 2; ++k) {
        fk = Fin[k];
        fnkc.r = Fin[ncfft - k].r;
        fnkc.i = -Fin[ncfft - k].i;

        C_ADD(fek, fk, fnkc);
        C_SUB(tmp, fk, fnkc);
        C_MUL(fok, tmp, SUPER_TWIDDLE(k - 1, st) /* st->super_twiddles[k - 1] */);
        C_ADD(Fout[k], fek, fok);
        C_SUB(Fout[ncfft - k], fek, fok);
        Fout[ncfft - k].i *= -1;
    }
    return ret;
}

static kfft_return_t
evali_func(kfft_plan_sclr* stu, const kfft_cpx* freqdata, kfft_scalar* timedata, kfft_cpx* tmpbuf) {
    kfft_return_t ret = KFFT_RET_SUCCESS;
    uint32_t ncfft = stu->nfft;

    KFFT_TMP_ZEROMEM(tmpbuf, sizeof(kfft_cpx) * ncfft);

    ret = eval_inverse_internal(stu, freqdata, tmpbuf);
    if (ret == KFFT_RET_SUCCESS) {
        ret = kfft_eval_cpx(stu->substate, tmpbuf, tmpbuf);
        if (ret == KFFT_RET_SUCCESS) {
            for (uint32_t i = 0; i < ncfft; i++) {
                timedata[i] = tmpbuf[i].r / 2;
            }
        }
    }
    return ret;
}

kfft_return_t
kfft_evali_scalar_internal(kfft_plan_sclr* stu, const kfft_cpx* freqdata, kfft_scalar* timedata,
                           kfft_cpx* tmpbuf) {
    kfft_return_t ret = KFFT_RET_SUCCESS;
    kfft_plan_sclr* st = (kfft_plan_sclr*)stu;

    if (!(st->substate->flags & KFFT_FLAG_INVERSE)) {
        return KFFT_RET_IMPROPER_PLAN;
    }

    uint32_t ncfft = st->nfft;

    if (tmpbuf == NULL) {
        kfft_cpx* tbuf = KFFT_TMP_ALLOC(sizeof(kfft_cpx) * ncfft, KFFT_PLAN_ALIGN(stu));
        if (tbuf) {
            ret = evali_func(stu, freqdata, timedata, tbuf);
            KFFT_TMP_FREE(tbuf, KFFT_PLAN_ALIGN(stu));
        } else {
            ret = KFFT_RET_BUFFER_FAIL;
        }
    } else {
        ret = evali_func(stu, freqdata, timedata, tmpbuf);
    }
    return ret;
}

KFFT_API kfft_return_t
kfft_evali_scalar(kfft_plan_sclr* stu, const kfft_cpx* freqdata, kfft_scalar* timedata) {
    /* input buffer timedata is stored row-wise */
    return kfft_evali_scalar_internal(stu, freqdata, timedata, NULL);
}

kfft_return_t
kfft_eval_scalar_norm_internal(kfft_plan_sclr* cfg, const kfft_scalar* fin, kfft_scalar* fout,
                               kfft_cpx* ftemp) {
    kfft_return_t ret = KFFT_RET_SUCCESS;
    size_t memneed = 2 * sizeof(kfft_cpx) * cfg->nfft;

    kfft_cpx* fbuf = (ftemp == NULL) ? KFFT_TMP_ALLOC(memneed, KFFT_MMGR_ALIGN(cfg)) : ftemp;
    if (fbuf) {
        KFFT_TMP_ZEROMEM(fbuf, memneed);
        kfft_cpx* ftmp = fbuf + cfg->nfft;

        ret = kfft_eval_scalar_internal((kfft_plan_sclr*)cfg, fin, fbuf, ftmp);
        if (ret == KFFT_RET_SUCCESS) {
            kfft_math_magnitude(fbuf, fout, cfg->nfft);
        }

        KFFT_TMP_FREE(fbuf, KFFT_MMGR_ALIGN(cfg));
    } else {
        ret = KFFT_RET_BUFFER_FAIL;
    }
    return ret;
}

KFFT_API kfft_return_t
kfft_eval_scalar_norm(kfft_plan_sclr* cfg, const kfft_scalar* timedata, kfft_scalar* data) {
    return kfft_eval_scalar_norm_internal(cfg, timedata, data, NULL);
}

/* ********************************************************************************
      TODO  Functionality
******************************************************************************** */
