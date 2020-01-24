#include "kfft.h"

/* perform the butterfly for one stage of a mixed radix FFT */
static void
kf_bfly_generic(kfft_cpx* Fout, const size_t fstride, const kfft_comp_t* st, uint32_t m,
                uint32_t p) {
    uint32_t u, k, q1, q;
    //    const kfft_cpx* twiddles = st->twiddles;
    kfft_cpx t;
    uint32_t Norig = st->nfft;

    kfft_trace("%s: m - %d | p - %d\n", "Generic FFT", m, p);

    kfft_cpx* scratch = (kfft_cpx*)KFFT_TMP_ALLOC(sizeof(kfft_cpx) * p);

    // TODO Maybe use Rader for FFT buffer

#if defined(KFFT_RADER_ALGO)
    if (!(p < KFFT_RADER_LIMIT)) {
        for (u = 0; u < m; ++u) {
            k = u;

            for (q1 = 0; q1 < p; ++q1) {
                // TODO Create buffer
                k += m;
            }

            k = u;
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
                uint32_t twidx = 0;
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
