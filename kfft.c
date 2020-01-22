/*
Copyright (c) 2003-2010, Mark Borgerding
              2018-2019, Alexander Smirnov

All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice, this list of conditions
and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice, this list of
conditions and the following disclaimer in the documentation and/or other materials provided with
the distribution.
    * Neither the author nor the names of any contributors may be used to endorse or promote
products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "kfft_guts.h"
#include "kfft_core.h"
#include "kfft_alloc.h"
#include "kfft_trace.h"
/* The guts header contains all the multiplication and addition macros that are defined for
 fixed or floating point complex numbers.  It also delares the kf_ internal functions.
 */

void
kfft_info(kfft_info_t* info) {

    info->vmajor = KFFT_MAJOR_VERSION;
    info->vminor = KFFT_MINOR_VERSION;
    info->vpatch = KFFT_PATCH_VERSION;

#if defined(KFFT_TRACE)
    info->flags |= KFFT_FLAG_INFO_TRACE;
#endif
#if defined(KFFT_USE_SIMD)
    info->flags |= KFFT_FLAG_INFO_USE_SIMD;
#endif
#if defined(KFFT_USE_ALLOCA)
    info->flags |= KFFT_FLAG_INFO_USE_ALLOCA;
#endif
#if defined(KFFT_USE_SYSMATH)
    info->flags |= KFFT_FLAG_INFO_USE_SYSMATH;
#endif
#if defined(KFFT_RADER_ALGO)
    info->flags |= KFFT_FLAG_INFO_RADER_ALGO;
#endif
#if defined(KFFT_MEMLESS_MODE)
    info->flags |= KFFT_FLAG_INFO_MEMLESS_MODE;
#endif
}

uint32_t
kfft_version(void) {
    return (((uint32_t)(KFFT_MAJOR_VERSION)) << 16) | (KFFT_MINOR_VERSION);
}

uint32_t
kfft_next_fast_size(uint32_t n) {
    while (1) {
        uint32_t m = n;
        while ((m % 2) == 0)
            m /= 2;
        while ((m % 3) == 0)
            m /= 3;
        while ((m % 5) == 0)
            m /= 5;
        if (m <= 1)
            break; /* n is completely factorable by twos, threes, and fives */
        n++;
    }
    return n;
}

static inline size_t
kfft_calculate(const uint32_t nfft, const bool inverse_fft) {
    size_t ret = sizeof(kfft_plan_t) + sizeof(kfft_cpx) * (nfft * 3 / 2);
    size_t subsize = 0;
    kfft_kconfig(nfft, inverse_fft, 0, NULL, &subsize);

    ret += subsize;

    return ret;
}

/* ********************************************************************************
      TODO  Functionality
******************************************************************************** */

KFFT_API uintptr_t
kfft_config(const uint32_t nfft, const bool inverse_fft, const uintptr_t A, size_t* lenmem) {
    kfft_plan_t* st = NULL;

    kfft_pool_t* mmgr = NULL;
    bool flag_create = false;

    if (lenmem == NULL) {
        if (A == 0) {
            size_t memneeded = kfft_calculate(nfft, inverse_fft);

            mmgr = kfft_allocator_create(memneeded);
            flag_create = true;

            kfft_trace("[REAL] %s: %p\n", "Create new allocator and plan", (void*)mmgr);
        } else {
            mmgr = (kfft_pool_t*)A;
            kfft_trace("[REAL] %s: %p\n", "Use allocator and create plan", (void*)mmgr);
        }

        if (mmgr)
            st = kfft_internal_alloc(mmgr, sizeof(kfft_plan_t));
    } else {
        size_t memneeded = kfft_calculate(nfft, inverse_fft);
        if (A && *lenmem >= memneeded) {
            mmgr = (kfft_pool_t*)A;

            kfft_allocator_clear(mmgr);
            st = kfft_internal_alloc(mmgr, sizeof(kfft_plan_t));

            kfft_trace("[REAL] %s: %p\n", "Reuse allocator and create plan", (void*)mmgr);
        }
        *lenmem = memneeded;
    }

    if (!st) {
    bailout:
        if (mmgr && (flag_create == true)) {
            kfft_allocator_free(&mmgr);
        }
        return 0;
    }

    st->mmgr = mmgr;
    st->substate = kfft_kconfig(nfft, inverse_fft, 0, mmgr, NULL);
    if (st->substate == NULL)
        goto bailout;

    st->tmpbuf = kfft_internal_alloc(st->mmgr, sizeof(kfft_cpx) * nfft);
    if (st->tmpbuf == NULL)
        goto bailout;

    // TODO Maybe memless
    if (nfft > 1) {
        st->super_twiddles = kfft_internal_alloc(st->mmgr, sizeof(kfft_cpx) * (nfft / 2));
        if (st->super_twiddles == NULL)
            goto bailout;
    }

    for (uint32_t i = 0; i < nfft / 2; ++i) {
        double phase = -KFFT_CONST_PI * ((double)(i + 1) / nfft + .5);
        if (inverse_fft)
            phase *= -1;
        kf_cexp(st->super_twiddles + i, phase);
    }

#if defined(KFFT_TRACE)
    kfft_trace("%s: ", "Factors");
    for (uint32_t i = 0; i < st->substate->fac_count; i++) {
        kfft_trace("%d ", st->substate->factors[i]);
    }
    kfft_trace("%s\n", "");
    #if defined(KFFT_RADER_ALGO) && !defined(KFFT_MEMLESS_MODE)
    kfft_trace("%s: ", "Prime roots");
    for (int i = 0; i < st->substate->prm_count; i++) {
        kfft_trace("%d ", st->substate->primes[i].prime);
    }
    #endif
    kfft_trace("%s\n", "");
#endif /* TRACE */
    return (uintptr_t)st;
}

KFFT_API void
kfft(uintptr_t stu, const kfft_scalar* timedata, kfft_cpx* freqdata) {
    /* input buffer timedata is stored row-wise */
    uint32_t k, ncfft;
    kfft_cpx fpnk, fpk, f1k, f2k, tw, tdc;

    kfft_plan_t* st = (kfft_plan_t*)stu;

    if (st->substate->inverse) {
        kfft_trace("%s\n", "kiss fft usage error: improper alloc");
        exit(1);
    }

    ncfft = st->substate->nfft;

    for (uint32_t i = 0; i < ncfft; i++) {
        st->tmpbuf[i].r = timedata[i];
        st->tmpbuf[i].i = 0;
    }

    __kfft(st->substate, st->tmpbuf, st->tmpbuf);

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
kffti(uintptr_t stu, const kfft_cpx* freqdata, kfft_scalar* timedata) {
    /* input buffer timedata is stored row-wise */
    uint32_t k, ncfft;

    kfft_plan_t* st = (kfft_plan_t*)stu;

    if (st->substate->inverse == 0) {
        kfft_trace("%s\n", "kiss fft usage error: improper alloc");
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
    __kfft(st->substate, st->tmpbuf, st->tmpbuf);

    for (uint32_t i = 0; i < ncfft; i++) {
        timedata[i] = S_DIV(st->tmpbuf[i].r, (2 * st->substate->nfft));
    }
}

KFFT_API void
kfft_free(uintptr_t* cfg) {
    if (cfg && *cfg) {
        kfft_trace("Cleanup plan: %p\n", (void*)(*cfg));
        kfft_plan_t* st = (kfft_plan_t*)(*cfg);
        if (st->mmgr != NULL) {
            kfft_allocator_free(&(st->mmgr));
        }
        *cfg = 0;
    }
}
/* ********************************************************************************
      TODO  Functionality
******************************************************************************** */
