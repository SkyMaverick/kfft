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
#include "_kfft_bf.h"
/* The guts header contains all the multiplication and addition macros that are defined for
 fixed or floating point complex numbers.  It also delares the kf_ internal functions.
 */

static
void kf_work(
        kiss_fft_cpx * Fout,
        const kiss_fft_cpx * f,
        const size_t fstride,
        int in_stride,
        int * factors,
        const __fft_cfg st
        )
{
    kiss_fft_cpx * Fout_beg=Fout;
    const int p=*factors++; /* the radix  */
    const int m=*factors++; /* stage's fft length/p */
    const kiss_fft_cpx * Fout_end = Fout + p*m;

#ifdef _OPENMP
    // use openmp extensions at the
    // top-level (not recursive)
    if (fstride==1 && p<=5)
    {
        int k;

        // execute the p different work units in different threads
#       pragma omp parallel for
        for (k=0;k<p;++k)
            kf_work( Fout +k*m, f+ fstride*in_stride*k,fstride*p,in_stride,factors,st);
        // all threads have joined by this point

        switch (p) {
            case 2: kf_bfly2(Fout,fstride,st,m); break;
            case 3: kf_bfly3(Fout,fstride,st,m); break;
            case 4: kf_bfly4(Fout,fstride,st,m); break;
            case 5: kf_bfly5(Fout,fstride,st,m); break;
            default: kf_bfly_generic(Fout,fstride,st,m,p); break;
        }
        return;
    }
#endif

    if (m==1) {
        do{
            *Fout = *f;
            f += fstride*in_stride;
        }while(++Fout != Fout_end );
    }else{
        do{
            // recursive call:
            // DFT of size m*p performed by doing
            // p instances of smaller DFTs of size m,
            // each one takes a decimated version of the input
            kf_work( Fout , f, fstride*p, in_stride, factors,st);
            f += fstride*in_stride;
        }while( (Fout += m) != Fout_end );
    }

    Fout=Fout_beg;

    // recombine the p smaller DFTs
    switch (p) {
        case 2: kf_bfly2(Fout,fstride,st,m); break;
        case 3: kf_bfly3(Fout,fstride,st,m); break;
        case 4: kf_bfly4(Fout,fstride,st,m); break;
        case 5: kf_bfly5(Fout,fstride,st,m); break;
        default: kf_bfly_generic(Fout,fstride,st,m,p); break;
    }
}

/*  facbuf is populated by p1,m1,p2,m2, ...
    where
    p[i] * m[i] = m[i-1]
    m0 = n                  */
static
void kf_factor(int n,int * facbuf)
{
    int p=4;
    double floor_sqrt;
    floor_sqrt = floor( sqrt((double)n) );

    /*factor out powers of 4, powers of 2, then any remaining primes */
    do {
        while (n % p) {
            switch (p) {
                case 4: p = 2; break;
                case 2: p = 3; break;
                default: p += 2; break;
            }
            if (p > floor_sqrt)
                p = n;          /* no more factors, skip to end */
        }
        n /= p;
        *facbuf++ = p;
        *facbuf++ = n;
    } while (n > 1);
}

static inline void
__kiss_fft_init (__fft_cfg   st,
                 int         nfft,
                 int         inverse_fft,
                 int         delta,
                 int         step)
{
        int i;
        st->nfft=nfft;
        st->inverse = inverse_fft;
        st->delta   = delta;
        st->step    = step;

        kf_factor(nfft,st->factors);

        for (i=0;i<nfft;++i) {
            const double pi=3.141592653589793238462643383279502884197169399375105820974944;
            double phase = -2*pi*i / nfft;
            if (st->inverse)
                phase *= -1;
            kf_cexp(st->twiddles+i, phase );
        }
}

/*
 *
 * User-callable function to allocate all necessary storage space for the fft.
 *
 * The return value is a contiguous block of memory, allocated with malloc.  As such,
 * It can be freed with free(), rather than a kiss_fft-specific function.
 * */
static __fft_cfg
__kiss_fft_config (int nfft,
                   int inverse_fft,
                   int delta,
                   int step,
                   void * mem,
                   size_t * lenmem )
{
    __fft_cfg st=NULL;
    size_t memneeded = sizeof(struct kiss_fft_state)
        + sizeof(kiss_fft_cpx)*(nfft-1); /* twiddle factors*/

    if ( lenmem==NULL) {
        st = ( __fft_cfg)KISS_FFT_MALLOC( memneeded );
    }else{
        if (mem != NULL && *lenmem >= memneeded)
            st = (__fft_cfg)mem;
        *lenmem = memneeded;
    }
    if (st) {
        __kiss_fft_init(st, nfft, inverse_fft, delta, step);
    }
    return st;
}


static inline void
__kiss_fft_stride(__fft_cfg st,const kiss_fft_cpx *fin,kiss_fft_cpx *fout,int in_stride)
{
    if (fin == fout) {
        //NOTE: this is not really an in-place FFT algorithm.
        //It just performs an out-of-place FFT into a temp buffer
        kiss_fft_cpx * tmpbuf = (kiss_fft_cpx*)KISS_FFT_TMP_ALLOC( sizeof(kiss_fft_cpx)*st->nfft);
        kf_work(tmpbuf,fin,1,in_stride, st->factors,st);
        memcpy(fout,tmpbuf,sizeof(kiss_fft_cpx)*st->nfft);
        KISS_FFT_TMP_FREE(tmpbuf);
    }else{
        kf_work( fout, fin, 1,in_stride, st->factors,st );
    }
}

static void
__kiss_fft(__fft_cfg cfg,const kiss_fft_cpx *fin,kiss_fft_cpx *fout)
{
    __kiss_fft_stride (cfg,fin,fout,1);
}

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

        __kiss_fft_config (nfft, inverse_fft, delta, step, NULL, &subsize);
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
        __kiss_fft_config(nfft, inverse_fft, delta, step, st->substate, &subsize);

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
