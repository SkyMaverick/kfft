#pragma once

#include "_kfft_guts.h"
#include "kfft_core.h"

#if defined (USE_RADER_ALGO)

static inline int
_kfr_ispower ( unsigned a,
               unsigned b,
               unsigned p )
{
    uint64_t res = 1;     // Initialize result
    uint64_t x = (uint64_t)a;
    uint64_t y = (uint64_t)b;

    x = x % p; // Update x if it is more than or
    // equal to p

    while (y > 0)
    {
        // If y is odd, multiply x with result
        if (y & 1)
            res = (res * x) % p; // FIXME Maybe integer overflow

        // y must be even now
        y = y >> 1; // y = y/2
        x = (x * x) % p; // FIXME Maybe integer overflow

    }
    return (res != 1) ? 1 : 0;
}

// WARNING in kfft num - always prime. Don't check this
static unsigned
_kfr_find_root (unsigned num) {
    unsigned phi = num - 1;
    unsigned n = phi;
    for (unsigned i = 2; i*i <= n; i++) {
        if (n % i == 0){
            // TODO nums_push(i);
            while (n % i == 0)
                n /= i;
        }
    }
    if (n > 1)
        nums_push(n);

//    printf ("Primes root candidate:\t");
//    for (unsigned i= 0; i < primes.count; i++) {
//        printf ("%d\t", primes.prime[i]);
//    }
//    printf ("\n");
//
    for (unsigned res=2; res<=num; ++res) {
		int ok = 1;
		for (unsigned i = 0; primes.count > i && ok; ++i) {
            int k = nums_get(i);
            kfft_trace ("%d\t%d\t%d\n", res, i, k);
			ok &= _kfr_ispower (res, phi / k, num);
        }
        if (ok)  return res;
	}
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

