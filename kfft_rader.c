#include "kfft_guts.h"

#if defined(KFFT_RADER_ALGO)
    #include "kfft_core.h"
    #include "kfft_rader.h"

void
kf_rader(kfft_cpx* Fout, const uint32_t fstride, const kfft_kplan_t* st, uint32_t p) {
    //    kfft_cpx x0 = {.r = Fout->r, .i = Fout->i}; // save first element value
    //    kfft_cpx* X0 = Fout;                        // remember first element
    //    *Fout = X0 * p;
    //
    //    uint32_t nfft = p - 1;
    //
    //    Fout++; // select N-1 size for sequence
    //
    //    // Recursive work procedure
    //
    // kfft_trace("[RADER] pr - %u | ipr - %u\n", pr, ipr);
}

#endif /* KFFT_READER_ALGO */
