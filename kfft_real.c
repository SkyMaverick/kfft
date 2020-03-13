/* The guts header contains all the multiplication and addition macros that are defined for
 fixed or floating point complex numbers.  It also delares the kf_ internal functions.
 */
#include "kfft.h"

#include "kfft_math.h"
#include "kfft_trace.h"

#ifndef KFFT_MEMLESS_MODE
    #define SUPER_TWIDDLE(i, P) P->super_twiddles[i]
#else
static inline kfft_cpx
get_super_twiddle(uint32_t i, kfft_plan_t* P) {
    kfft_cpx ret;

    kfft_scalar phase = -KFFT_CONST_PI * ((kfft_scalar)(i + 1) / P->substate->nfft + .5);
    if (P->substate->flags & KFFT_FLAG_INVERSE)
        phase *= -1;
    kf_cexp(&ret, phase);
    return ret;
}
    #define SUPER_TWIDDLE(i, P) get_super_twiddle(i, P)
#endif

static inline size_t
kfft_calculate(const uint32_t nfft, const uint32_t flags) {
    size_t ret = sizeof(kfft_real_t) + sizeof(kfft_cpx) * (nfft * 3 / 2);
    size_t subsize = 0;

    kfft_config_cpx(nfft, flags, 0, NULL, &subsize);

    ret += subsize;

    return ret;
}

#ifdef KFFT_TRACE

static void
kfft_trace_plan(kfft_real_t* P) {
    kfft_trace("[REAL] %s: %p", "Create KFFT real plan", (void*)P);
    kfft_trace("\n\t %s - %p", "Uses complex plan", (void*)(P->substate));
    kfft_trace("\n\t %s - %p", "Temporary buffer", (void*)(P->tmpbuf));
    kfft_trace("\n\t %s - %p\n", "Real twiddles", (void*)(P->super_twiddles));
}

#endif
/* ********************************************************************************
      TODO  Functionality
******************************************************************************** */

KFFT_API kfft_real_t*
kfft_config_real(const uint32_t nfft, const uint32_t flags, const kfft_pool_t* A, size_t* lenmem) {
    kfft_real_t* st = NULL;

    kfft_pool_t* mmgr = NULL;
    bool flag_create = false;

    if (lenmem == NULL) {
        if (A == 0) {
            size_t memneeded = kfft_calculate(nfft, flags);

            mmgr = kfft_allocator_create(memneeded);
            flag_create = true;

            kfft_trace("[REAL] %s: %p\n", "Create new allocator and plan", (void*)mmgr);
        } else {
            mmgr = (kfft_pool_t*)A;
            kfft_trace("[REAL] %s: %p\n", "Use allocator and create plan", (void*)mmgr);
        }

        if (mmgr)
            st = kfft_internal_alloc(mmgr, sizeof(kfft_real_t));
    } else {
        size_t memneeded = kfft_calculate(nfft, flags);
        if (A && *lenmem >= memneeded) {
            mmgr = (kfft_pool_t*)A;

            if (flags & KFFT_FLAG_RENEW)
                kfft_allocator_clear(mmgr);

            st = kfft_internal_alloc(mmgr, sizeof(kfft_real_t));

            kfft_trace("[REAL] %s: %p\n", "Reuse allocator and create plan", (void*)mmgr);
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

    st->tmpbuf = kfft_internal_alloc(st->object.mmgr, sizeof(kfft_cpx) * nfft);
    if (st->tmpbuf == NULL)
        goto bailout;

    // TODO Maybe memless
    if (nfft > 1) {
        st->super_twiddles = kfft_internal_alloc(st->object.mmgr, sizeof(kfft_cpx) * (nfft / 2));
        if (st->super_twiddles == NULL)
            goto bailout;
    }

    st->substate = kfft_config_cpx(nfft, flags | (!(KFFT_FLAG_RENEW)), 0, st->object.mmgr, NULL);
    if (st->substate == NULL)
        goto bailout;

    for (uint32_t i = 0; i < nfft / 2; ++i) {
        double phase = -KFFT_CONST_PI * ((double)(i + 1) / nfft + .5);
        if (flags & KFFT_FLAG_INVERSE)
            phase *= -1;
        kf_cexp(st->super_twiddles + i, phase);
    }

#ifdef KFFT_TRACE
    kfft_trace_plan(st);
#endif

    return st;
}

KFFT_API void
kfft_eval_real(kfft_real_t* stu, const kfft_scalar* timedata, kfft_cpx* freqdata) {
    /* input buffer timedata is stored row-wise */
    uint32_t k, ncfft;
    kfft_cpx fpnk, fpk, f1k, f2k, tw, tdc;

    kfft_real_t* st = (kfft_real_t*)stu;

    if (st->substate->flags & KFFT_FLAG_INVERSE) {
        kfft_trace("[REAL] %s\n", "kiss fft usage error: improper alloc");
        exit(1);
    }

    ncfft = st->substate->nfft;

    for (uint32_t i = 0; i < ncfft; i++) {
        st->tmpbuf[i].r = timedata[i];
        st->tmpbuf[i].i = 0;
    }

    kfft_eval_cpx(st->substate, st->tmpbuf, st->tmpbuf);

    tdc.r = st->tmpbuf[0].r;
    tdc.i = st->tmpbuf[0].i;
    freqdata[0].r = tdc.r + tdc.i;
    freqdata[ncfft].r = tdc.r - tdc.i;
#ifdef KFFT_USE_SIMD
    freqdata[ncfft].i = freqdata[0].i = _mm_set1_ps(0);
#else
    freqdata[ncfft].i = freqdata[0].i = 0;
#endif

    for (k = 1; k <= ncfft / 2; ++k) {
        fpk = st->tmpbuf[k];
        fpnk.r = st->tmpbuf[ncfft - k].r;
        fpnk.i = -st->tmpbuf[ncfft - k].i;

        C_ADD(f1k, fpk, fpnk);
        C_SUB(f2k, fpk, fpnk);
        C_MUL(tw, f2k, SUPER_TWIDDLE(k - 1, st) /* st->super_twiddles[k - 1] */);

        freqdata[k].r = HALF_OF(f1k.r + tw.r);
        freqdata[k].i = HALF_OF(f1k.i + tw.i);
        freqdata[ncfft - k].r = HALF_OF(f1k.r - tw.r);
        freqdata[ncfft - k].i = HALF_OF(tw.i - f1k.i);
    }
}

KFFT_API void
kfft_evali_real(kfft_real_t* stu, const kfft_cpx* freqdata, kfft_scalar* timedata) {
    /* input buffer timedata is stored row-wise */
    uint32_t k, ncfft;

    kfft_real_t* st = (kfft_real_t*)stu;

    if (!(st->substate->flags & KFFT_FLAG_INVERSE)) {
        kfft_trace("[REAL] %s\n", "kiss fft usage error: improper alloc");
        exit(1);
    }

    ncfft = st->substate->nfft;

    st->tmpbuf[0].r = freqdata[0].r + freqdata[ncfft].r;
    st->tmpbuf[0].i = freqdata[0].r - freqdata[ncfft].r;

    for (k = 1; k <= ncfft / 2; ++k) {
        kfft_cpx fk, fnkc, fek, fok, tmp;
        fk = freqdata[k];
        fnkc.r = freqdata[ncfft - k].r;
        fnkc.i = -freqdata[ncfft - k].i;

        C_ADD(fek, fk, fnkc);
        C_SUB(tmp, fk, fnkc);
        C_MUL(fok, tmp, SUPER_TWIDDLE(k - 1, st) /* st->super_twiddles[k - 1] */);
        C_ADD(st->tmpbuf[k], fek, fok);
        C_SUB(st->tmpbuf[ncfft - k], fek, fok);
#ifdef KFFT_USE_SIMD
        st->tmpbuf[ncfft - k].i *= _mm_set1_ps(-1.0);
#else
        st->tmpbuf[ncfft - k].i *= -1;
#endif
    }
    kfft_eval_cpx(st->substate, st->tmpbuf, st->tmpbuf);

    for (uint32_t i = 0; i < ncfft; i++) {
        timedata[i] = st->tmpbuf[i].r / 2;
    }
}
/* ********************************************************************************
      TODO  Functionality
******************************************************************************** */
