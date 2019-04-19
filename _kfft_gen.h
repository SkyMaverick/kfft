#pragma once

#include "_kfft_guts.h"
#include "kfft_core.h"

#if defined (USE_RADER_ALGO)
    // TODO Rader specific functions here
#endif /* USE_RADER_ALGO */

/* perform the butterfly for one stage of a mixed radix FFT */
static void kf_bfly_generic(
        kfft_cpx * Fout,
        const size_t fstride,
        const __fft_cfg st,
        int m,
        int p
        )
{
    int u,k,q1,q;
    kfft_cpx * twiddles = st->twiddles;
    kfft_cpx t;
    int Norig = st->nfft;

    kfft_trace ("%s: m - %d | p - %d\n", "Generic FFT", m, p);

    kfft_cpx * scratch = (kfft_cpx*)KFFT_TMP_ALLOC(sizeof(kfft_cpx)*p);

    // TODO Maybe use Rader for FFT buffer

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
        __fft_cfg tmp_cfg = __kfft_config (p, st->inverse, st->level + 1, NULL, NULL);
        for ( u=0; u<m; ++u ) {
            k=u;
            for ( q1=0 ; q1<p ; ++q1 ) {
                scratch[q1] = Fout[ k  ];
                k += m;
            }

            k=u;
            // TODO Rader here
        }
        __kfft_free(tmp_cfg);
    }
#endif /* USE_RADER_ALGO */
    KFFT_TMP_FREE(scratch);
}

