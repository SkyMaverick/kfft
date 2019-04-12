#pragma once

#include "_kfft_guts.h"
#include "kfft_core.h"

static inline __fft_cfg
__fft_recursive_alloc (const __fft_cfg st,
                       int p)
{
    return __kiss_fft_config ( p, st->inverse, st->delta, st->step, st->level + 1, NULL, NULL);
}

static void kf_rader (
        kiss_fft_cpx * Fout,
        const size_t fstride,
        const __fft_cfg st,
        int m,
        int p
        )
{
    kfft_trace ("%s:\t m: %d\tp: %d\tlvl: %d\n", "Generic FFT use Rader", m, p, st->level);
    // TODO Add Rader algo and strart recursive FFT with new buffer
    __fft_cfg st_req = __fft_recursive_alloc (st, p);
    
    // TODO Traditional FFT with new buffer
    
    __kiss_fft_free(st_req);
}

/* perform the butterfly for one stage of a mixed radix FFT */
static void kf_bfly_generic(
        kiss_fft_cpx * Fout,
        const size_t fstride,
        const __fft_cfg st,
        int m,
        int p
        )
{
    int u,k,q1,q;
    kiss_fft_cpx * twiddles = st->twiddles;
    kiss_fft_cpx t;
    int Norig = st->nfft;

    kfft_trace ("%s:\t m: %d\tp: %d\n", "Generic FFT", m, p);

    kiss_fft_cpx * scratch = (kiss_fft_cpx*)KISS_FFT_TMP_ALLOC(sizeof(kiss_fft_cpx)*p);

    for ( u=0; u<m; ++u ) {
        k=u;
        for ( q1=0 ; q1<p ; ++q1 ) {
            scratch[q1] = Fout[ k  ];
            k += m;
        }

        k=u;
        for ( q1=0 ; q1<p ; ++q1 ) {
            int twidx=0;
            Fout[ k ] = scratch[0];
            for (q=1;q<p;++q ) {
                twidx += fstride * k;
                if (twidx>=Norig) twidx-=Norig;
                C_MUL(t,scratch[q] , twiddles[twidx] );
                C_ADDTO( Fout[ k ] ,t);
            }
            k += m;
        }
    }
    KISS_FFT_TMP_FREE(scratch);
}

