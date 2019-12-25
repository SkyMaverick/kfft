#if defined(KFFT_RADER_ALGO)

#include "kfft_guts.h"

#include "kfft_core.h"
#include "kfft_rader.h"

// kfft_rplan_t*
// kfft_rconfig(int nfft, int inverse_fft, int level, void* mem, size_t* lenmem) {
//     kfft_trace("%s: %d\n", "[RADER] Create RADER plan level", level);
//     kfft_rplan_t* st = NULL;
//     size_t subsize = 0;
//
//     kfft_kconfig(nfft, inverse_fft, level, NULL, &subsize);
// #ifndef KFFT_MEMLESS_MODE
//     size_t memneeded = sizeof(kfft_rplan_t) * 2 + subsize + sizeof(kfft_cpx) * (nfft * 3 / 2);
// #else
//     size_t memneeded = sizeof(kfft_rplan_t) * 2 + subsize + sizeof(kfft_cpx) * nfft;
// #endif /* memless */
//
//     if (lenmem == NULL) {
//         st = (kfft_rplan_t*)KFFT_MALLOC(memneeded);
//     } else {
//         if (*lenmem >= memneeded) {
//             st = (kfft_rplan_t*)mem;
//         }
//         *lenmem = memneeded;
//     }
//     if (!st)
//         return NULL;
//
//     size_t delta = sizeof(unsigned) * 2 + 1;
//     st->substate = (kfft_kplan_t*)(st + delta);
//     kfft_kconfig(nfft, inverse_fft, level, st->substate, &subsize);
//
//     st->pr = _kfr_prime(nfft);
//     st->ipr = _kfr_inverse_prime(nfft, st->pr);
//
// #if defined(TRACE)
//     kfft_trace("[RADER] %s: %zu\n", "Memory allocate", memneeded);
//     kfft_trace("[RADER] %s: ", "Factors");
//     for (int i = 0; st->substate->factors[i] != 0; i++) {
//         kfft_trace("%d ", st->substate->factors[i]);
//     }
//     kfft_trace("%s\n", "");
//     kfft_trace("[RADER] prime - %d | inverse prime - %d\n", st->pr, st->ipr);
// #endif
//
//     return st;
// }
//
void
kf_rader(kfft_cpx* Fout, const uint32_t fstride, const kfft_kplan_t* st, uint32_t m, uint32_t p) {
    kfft_cpx x0 = {.r = Fout->r, .i = Fout->i}; // save first element value
    kfft_cpx* X0 = Fout;                        // remember first element
    unsigned pr, ipr;

    Fout++; // select N-1 size for sequence
#ifndef KFFT_MEMLESS_MODE
    // select prevalte root / inverse
    const unsigned* roots = st->roots;
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

// void
// kfft_rfree(kfft_rplan_t** cfg) {
//     if (cfg && *cfg) {
//         KFFT_FREE(*cfg);
//         *cfg = NULL;
//     }
// }

#endif /* KFFT_READER_ALGO */
