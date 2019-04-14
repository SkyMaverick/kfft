#pragma once

#include "_kfft_guts.h"
#include "kfft_core.h"

#if defined (USE_RADER_ALGO) 

static inline int
__kfr_prime_root (int p)
{
    return 0;
}

#endif /* USE_RADER_ALGO */





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
    
    // TODO Maybe use Rader for FFT bufer
#if defined ( USE_RADER_ALGO )
    if (( p < KFFT_RADER_LIMIT) || (st->level > KFFT_RADER_LEVEL)) {
#endif
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
#if defined (USE_RADER_ALGO)
    } else {
        __fft_cfg tmp_cfg = __kiss_fft_config (p, st->inverse, 0, 0, st->level + 1, NULL, NULL);
        for ( u=0; u<m; ++u ) {
            k=u;
            for ( q1=0 ; q1<p ; ++q1 ) {
                scratch[q1] = Fout[ k  ];
                k += m;
            }

            k=u;
            // TODO Rader here
        }
        __kiss_fft_free(tmp_cfg);
    }
#endif
    KISS_FFT_TMP_FREE(scratch);
}

