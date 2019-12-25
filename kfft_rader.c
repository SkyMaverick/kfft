#if defined(KFFT_RADER_ALGO)

#include "kfft_guts.h"

#include "kfft_core.h"
#include "kfft_rader.h"

void
kf_rader(kfft_cpx* Fout, const uint32_t fstride, const kfft_kplan_t* st, uint32_t m, uint32_t p) {
    kfft_cpx x0 = {.r = Fout->r, .i = Fout->i}; // save first element value
    kfft_cpx* X0 = Fout;                        // remember first element
    uint32_t pr, ipr;

    Fout++; // select N-1 size for sequence
#ifndef KFFT_MEMLESS_MODE
    // select prevalte root / inverse
    const uint32_t* roots = st->roots;
    while (*roots) {
        if (*roots == p) { // FIXME
            pr = *(roots + 1);
            ipr = *(roots + 2);
        }
        roots += 3;
    }
#else
    pr = kfft_prime_root(p);
    ipr = kfft_primei_root(pr, p);
#endif
    uint32_t nfft = p - 1;

    kfft_trace("[RADER] pr - %u | ipr - %u\n", pr, ipr);
}

#endif /* KFFT_READER_ALGO */
