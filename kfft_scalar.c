/* The guts header contains all the multiplication and addition macros that are defined for
 fixed or floating point complex numbers.  It also delares the kf_ internal functions.
 */
#include "kfft.h"

#include "kfft_math.h"
#include "kfft_trace.h"

static inline kfft_cpx
kfft_sclr_twiddle(uint32_t i, const kfft_sclr_t* P) {
    kfft_cpx ret;

    kfft_scalar phase = -KFFT_CONST_PI * ((kfft_scalar)(i + 1) / P->substate->nfft + .5);
    if (P->substate->flags & KFFT_FLAG_INVERSE)
        phase *= -1;
    kf_cexp(&ret, phase);
    return ret;
}

#if defined(KFFT_MEMLESS_MODE)
    #define SUPER_TWIDDLE(i, P) kfft_sclr_twiddle(i, P)
#else
    #define SUPER_TWIDDLE(i, P) P->super_twiddles[i]
#endif /* KFFT_MEMLESS_MODE */

static inline size_t
kfft_calculate(const uint32_t nfft, const uint32_t flags) {
    size_t ret = sizeof(kfft_sclr_t);
#if !defined(KFFT_MEMLESS_MODE)
    ret += sizeof(kfft_cpx) * (nfft / 2);
#endif
    size_t subsize = 0;

    kfft_config_cpx(nfft, flags, NULL, &subsize);
    return ret + subsize;
}

#ifdef KFFT_TRACE

static void
kfft_trace_plan(kfft_sclr_t* P) {
    kfft_trace_scalar("%s: %p", "Create KFFT scalar plan", (void*)P);
    kfft_trace("\n\t %s - %p", "Uses complex plan", (void*)(P->substate));
    kfft_trace("\n\t %s - %p\n", "scalar twiddles", (void*)(P->super_twiddles));
}

#endif /* KFFT_TRACE */
/* ********************************************************************************
      TODO  Functionality
******************************************************************************** */

KFFT_API kfft_sclr_t*
kfft_config_scalar(const uint32_t nfft, const uint32_t flags, const kfft_pool_t* A,
                   size_t* lenmem) {
    kfft_sclr_t* st = NULL;

    kfft_pool_t* mmgr = NULL;
    bool flag_create = false;

    if (lenmem == NULL) {
        if (A == 0) {
            size_t memneeded = kfft_calculate(nfft, flags);

            mmgr = kfft_allocator_create(memneeded);
            flag_create = true;

            kfft_trace_scalar("%s: %p\n", "Create new allocator and plan", (void*)mmgr);
        } else {
            mmgr = (kfft_pool_t*)A;
            kfft_trace_scalar("%s: %p\n", "Use allocator and create plan", (void*)mmgr);
        }

        if (mmgr)
            st = kfft_internal_alloc(mmgr, sizeof(kfft_sclr_t));
    } else {
        size_t memneeded = kfft_calculate(nfft, flags);
        if (A && *lenmem >= memneeded) {
            mmgr = (kfft_pool_t*)A;

            if (flags & KFFT_FLAG_RENEW)
                kfft_allocator_clear(mmgr);

            st = kfft_internal_alloc(mmgr, sizeof(kfft_sclr_t));

            kfft_trace_scalar("%s: %p\n", "Reuse allocator and create plan", (void*)mmgr);
        }
        *lenmem = memneeded;
    }

    if (!st) {
    bailout:
        if (mmgr && (flag_create == true)) {
            kfft_allocator_free(mmgr);
        }
        return 0;
    }

    st->object.mmgr = mmgr;

    st->substate = kfft_config_cpx(nfft, KFFT_CHECK_FLAGS(flags), st->object.mmgr, NULL);
    if (st->substate == NULL)
        goto bailout;

#if !defined(KFFT_MEMLESS_MODE)
    if (nfft > 1) {
        st->super_twiddles = kfft_internal_alloc(st->object.mmgr, sizeof(kfft_cpx) * (nfft / 2));
        if (st->super_twiddles == NULL)
            goto bailout;
    }
    for (uint32_t i = 0; i < nfft / 2; ++i) {
        st->super_twiddles[i] = kfft_sclr_twiddle(i, st);
    }
#endif /* not KFFT_MEMLESS_MODE */

#ifdef KFFT_TRACE
    kfft_trace_plan(st);
#endif

    return st;
}

static inline kfft_return_t
eval_forward_internal(const kfft_sclr_t* st, const kfft_cpx* Fin, kfft_cpx* Fout) {
    kfft_return_t ret = KFFT_RET_SUCCESS;
    kfft_cpx fpnk, fpk, f1k, f2k, tw;

    uint32_t k, ncfft = st->substate->nfft;

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
eval_func(kfft_sclr_t* stu, kfft_cpx* tmpbuf, const kfft_scalar* timedata, kfft_cpx* freqdata) {
    kfft_return_t ret = KFFT_RET_SUCCESS;
    uint32_t ncfft = stu->substate->nfft;

    KFFT_TMP_ZEROMEM(tmpbuf, stu->substate->nfft * sizeof(kfft_cpx));
    for (uint32_t i = 0; i < ncfft; i++) {
        freqdata[i].r = timedata[i];
    }
    ret = kfft_eval_cpx(stu->substate, freqdata, tmpbuf);
    if (ret == KFFT_RET_SUCCESS) {
        ret = eval_forward_internal(stu, tmpbuf, freqdata);
    }

    return ret;
}

kfft_return_t
kfft_eval_scalar_internal(kfft_sclr_t* stu, const kfft_scalar* timedata, kfft_cpx* freqdata,
                          kfft_cpx* tmpbuf) {
    kfft_return_t ret = KFFT_RET_SUCCESS;
    /* input buffer timedata is stored row-wise */
    kfft_sclr_t* st = (kfft_sclr_t*)stu;

    if (st->substate->flags & KFFT_FLAG_INVERSE) {
        return KFFT_RET_IMPROPER_PLAN;
    }

    uint32_t ncfft = stu->substate->nfft;
    if (tmpbuf == NULL) {
        kfft_cpx* tbuf = KFFT_TMP_ALLOC(sizeof(kfft_cpx) * ncfft);
        if (tbuf) {
            ret = eval_func(stu, tbuf, timedata, freqdata);
            KFFT_TMP_FREE(tbuf);
        } else {
            ret = KFFT_RET_BUFFER_FAIL;
        }
    } else {
        ret = eval_func(stu, tmpbuf, timedata, freqdata);
    }

    return ret;
}

KFFT_API kfft_return_t
kfft_eval_scalar(kfft_sclr_t* stu, const kfft_scalar* timedata, kfft_cpx* freqdata) {
    return kfft_eval_scalar_internal(stu, timedata, freqdata, NULL);
}

static inline kfft_return_t
eval_inverse_internal(const kfft_sclr_t* st, const kfft_cpx* Fin, kfft_cpx* Fout) {
    kfft_return_t ret = KFFT_RET_SUCCESS;
    kfft_cpx fk, fnkc, fek, fok, tmp;

    uint32_t k, ncfft = st->substate->nfft;

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
        // #ifdef KFFT_USE_SIMD
        //         Fout[ncfft - k].i *= _mm_set1_ps(-1.0);
        // #else
        Fout[ncfft - k].i *= -1;
        // #endif
    }
    return ret;
}

static kfft_return_t
evali_func(kfft_sclr_t* stu, const kfft_cpx* freqdata, kfft_scalar* timedata, kfft_cpx* tmpbuf) {
    kfft_return_t ret = KFFT_RET_SUCCESS;
    uint32_t ncfft = stu->substate->nfft;

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
kfft_evali_scalar_internal(kfft_sclr_t* stu, const kfft_cpx* freqdata, kfft_scalar* timedata,
                           kfft_cpx* tmpbuf) {
    kfft_return_t ret = KFFT_RET_SUCCESS;
    kfft_sclr_t* st = (kfft_sclr_t*)stu;

    if (!(st->substate->flags & KFFT_FLAG_INVERSE)) {
        return KFFT_RET_IMPROPER_PLAN;
    }

    uint32_t ncfft = st->substate->nfft;

    if (tmpbuf == NULL) {
        kfft_cpx* tbuf = KFFT_TMP_ALLOC(sizeof(kfft_cpx) * ncfft);
        if (tbuf) {
            ret = evali_func(stu, freqdata, timedata, tbuf);
            KFFT_TMP_FREE(tbuf);
        } else {
            ret = KFFT_RET_BUFFER_FAIL;
        }
    } else {
        ret = evali_func(stu, freqdata, timedata, tmpbuf);
    }
    return ret;
}

KFFT_API kfft_return_t
kfft_evali_scalar(kfft_sclr_t* stu, const kfft_cpx* freqdata, kfft_scalar* timedata) {
    /* input buffer timedata is stored row-wise */
    return kfft_evali_scalar_internal(stu, freqdata, timedata, NULL);
}
/* ********************************************************************************
      TODO  Functionality
******************************************************************************** */
