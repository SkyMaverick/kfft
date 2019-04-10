/*
Copyright (c) 2003-2010, Mark Borgerding
              2018-2019, Alexander Smirnov

All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    * Neither the author nor the names of any contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#include "_kfft_guts.h"
#include "kfft_core.h"
/* The guts header contains all the multiplication and addition macros that are defined for
 fixed or floating point complex numbers.  It also delares the kf_ internal functions.
 */


int kiss_fft_next_fast_size(int n)
{
    while(1) {
        int m=n;
        while ( (m%2) == 0 ) m/=2;
        while ( (m%3) == 0 ) m/=3;
        while ( (m%5) == 0 ) m/=5;
        if (m<=1)
            break; /* n is completely factorable by twos, threes, and fives */
        n++;
    }
    return n;
}

/* ********************************************************************************
      TODO  Functionality
******************************************************************************** */

kiss_fft_cfg
kiss_fft_config  (int         nfft,
                  int         inverse_fft,
                  int         delta,
                  int         step,
                  void *      mem,
                  size_t *    lenmem)
{
    int i;
    kiss_fft_cfg st = NULL;
    size_t subsize = 0, memneeded = 0;

        __kiss_fft_config (nfft, inverse_fft, delta, step, 0, NULL, &subsize);
        memneeded = sizeof(struct kiss_fftr_state) + subsize + sizeof(kiss_fft_cpx) * ( nfft * 3 / 2);

        if (lenmem == NULL) {
            st = (kiss_fft_cfg) KISS_FFT_MALLOC (memneeded);
        } else {
            if (*lenmem >= memneeded) {
                st = (kiss_fft_cfg) mem;
            }
            *lenmem = memneeded;
        }
        if (!st)
            return NULL;

        st->substate = (__fft_cfg) (st + 1); /*just beyond kiss_fftr_state struct */
        st->tmpbuf = (kiss_fft_cpx *) (((char *) st->substate) + subsize);
        st->super_twiddles = st->tmpbuf + nfft;
        __kiss_fft_config(nfft, inverse_fft, delta, step, 0, st->substate, &subsize);

#if defined (TRACE)
        kfft_trace ("%s:\t%zu\n", "Memory allocate", memneeded);
        kfft_trace ("%s:\t", "Factors");
        for (i = 0; st->substate->factors[i] != 0; i++) {
            kfft_trace("%d ", st->substate->factors[i]);
        }
        kfft_trace ("\n", "");

#endif

        for (i = 0; i < nfft/2; ++i) {
            double phase =
                -3.14159265358979323846264338327 * ((double) (i+1) / nfft + .5);
            if (inverse_fft)
                phase *= -1;
            kf_cexp (st->super_twiddles+i,phase);
        }
    return st;
}

void
kiss_fft (kiss_fft_cfg               st,
          const kiss_fft_scalar     *timedata,
          kiss_fft_cpx              *freqdata)
{
    /* input buffer timedata is stored row-wise */
    int k,ncfft;
    kiss_fft_cpx fpnk,fpk,f1k,f2k,tw,tdc;

    if ( st->substate->inverse) {
        fprintf(stderr,"kiss fft usage error: improper alloc\n");
        exit(1);
    }

    ncfft = st->substate->nfft;

    if (st->substate->delta || st->substate->step) {
        for (int i=0; i < ncfft; i++) {
            st->tmpbuf[i].r = timedata [i * st->substate->step + st->substate->delta];
            st->tmpbuf[i].i = 0;
        }
    } else {
        for (int i=0; i < ncfft; i++) {
            st->tmpbuf[i].r = timedata [i];
            st->tmpbuf[i].i = 0;
        }
    }

    __kiss_fft( st->substate , st->tmpbuf, st->tmpbuf );

    tdc.r = st->tmpbuf[0].r;
    tdc.i = st->tmpbuf[0].i;
    freqdata[0].r = tdc.r + tdc.i;
    freqdata[ncfft].r = tdc.r - tdc.i;
#ifdef USE_SIMD
    freqdata[ncfft].i = freqdata[0].i = _mm_set1_ps(0);
#else
    freqdata[ncfft].i = freqdata[0].i = 0;
#endif

    for ( k=1;k <= ncfft/2 ; ++k ) {
        fpk    = st->tmpbuf[k];
        fpnk.r =   st->tmpbuf[ncfft-k].r;
        fpnk.i = - st->tmpbuf[ncfft-k].i;

        C_ADD( f1k, fpk , fpnk );
        C_SUB( f2k, fpk , fpnk );
        C_MUL( tw , f2k , st->super_twiddles[k-1]);

        freqdata[k].r = HALF_OF(f1k.r + tw.r);
        freqdata[k].i = HALF_OF(f1k.i + tw.i);
        freqdata[ncfft-k].r = HALF_OF(f1k.r - tw.r);
        freqdata[ncfft-k].i = HALF_OF(tw.i - f1k.i);
    }
}

void
kiss_ffti (kiss_fft_cfg         st,
           const kiss_fft_cpx    *freqdata,
           kiss_fft_scalar       *timedata)
{
    /* input buffer timedata is stored row-wise */
    int k, ncfft;

    if (st->substate->inverse == 0) {
        fprintf (stderr, "kiss fft usage error: improper alloc\n");
        exit (1);
    }

    ncfft = st->substate->nfft;

    st->tmpbuf[0].r = freqdata[0].r + freqdata[ncfft].r;
    st->tmpbuf[0].i = freqdata[0].r - freqdata[ncfft].r;

    for (k = 1; k <= ncfft / 2; ++k) {
        kiss_fft_cpx fk, fnkc, fek, fok, tmp;
        fk = freqdata[k];
        fnkc.r = freqdata[ncfft - k].r;
        fnkc.i = -freqdata[ncfft - k].i;

        C_ADD (fek, fk, fnkc);
        C_SUB (tmp, fk, fnkc);
        C_MUL (fok, tmp, st->super_twiddles[k-1]);
        C_ADD (st->tmpbuf[k],     fek, fok);
        C_SUB (st->tmpbuf[ncfft - k], fek, fok);
#ifdef USE_SIMD
        st->tmpbuf[ncfft - k].i *= _mm_set1_ps(-1.0);
#else
        st->tmpbuf[ncfft - k].i *= -1;
#endif
    }
    __kiss_fft (st->substate, st->tmpbuf, st->tmpbuf);

    if (st->substate->delta || st->substate->step) {
        for (int i=0; i < ncfft; i++) {
            timedata [i * st->substate->step + st->substate->delta] = S_DIV(st->tmpbuf[i].r ,(2 * st->substate->nfft));
        }
    } else {
        for (int i=0; i < ncfft; i++) {
            timedata[i] = S_DIV(st->tmpbuf[i].r, (2 * st->substate->nfft));
        }
    }

}

void
kiss_fft_free (kiss_fft_cfg* cfg) {
    if (cfg && *cfg) {
        free(*cfg);
        *cfg = NULL;
    }
}
/* ********************************************************************************
      TODO  Functionality
******************************************************************************** */
