#include "kfft_core.h"

#ifdef KFFT_RADER_ALGO
    #include "kfft_rader.h"
// TODO Rader specific functions here
#endif /* KFFT_RADER_ALGO */

/* perform the butterfly for one stage of a mixed radix FFT */
static void
kf_bfly_generic(kfft_cpx* Fout, const size_t fstride, const kfft_kplan_t* st, int m, int p) {
    int u, k, q1, q;
    //    const kfft_cpx* twiddles = st->twiddles;
    kfft_cpx t;
    int Norig = st->nfft;

    kfft_trace("%s: m - %d | p - %d\n", "Generic FFT", m, p);

    kfft_cpx* scratch = (kfft_cpx*)KFFT_TMP_ALLOC(sizeof(kfft_cpx) * p);

    // TODO Maybe use Rader for FFT buffer

#if defined(KFFT_RADER_ALGO)
    if (!((p < KFFT_RADER_LIMIT) || (st->level > KFFT_RADER_LEVEL))) {
        for (u = 0; u < m; ++u) {
            k = u;
            for (q1 = 0; q1 < p; ++q1) {
                scratch[q1] = Fout[k];
                k += m;
            }

            k = u;
            kf_rader(Fout, fstride, st, m, q1);
            // TODO Rader here
        }
    } else {
#endif /* KFFT_RADER_ALGO */
        for (u = 0; u < m; ++u) {
            k = u;
            for (q1 = 0; q1 < p; ++q1) {
                scratch[q1] = Fout[k];
                k += m;
            }

            k = u;
            for (q1 = 0; q1 < p; ++q1) {
                int twidx = 0;
                Fout[k] = scratch[0];
                for (q = 1; q < p; ++q) {
                    twidx += fstride * k;
                    if (twidx >= Norig)
                        twidx -= Norig;
                    C_MUL(t, scratch[q], TWIDDLE(twidx, st) /*twiddles[twidx]*/);
                    C_ADDTO(Fout[k], t);
                }
                k += m;
            }
        }
#if defined(KFFT_RADER_ALGO)
    }
#endif /* KFFT_RADER_ALGO */
    KFFT_TMP_FREE(scratch);
}
