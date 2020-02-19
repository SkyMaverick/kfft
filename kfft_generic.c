/* perform the butterfly for one stage of a mixed radix FFT */

#if defined(KFFT_RADER_ALGO)

static inline void
rader_method_eval(kfft_cpx* Fout, kfft_cpx* Ftmp, const size_t fstride, const kfft_comp_t* st,
                  uint32_t u, uint32_t m, uint32_t p) {
    uint32_t k = u, q1, idx;
    kfft_cpx x0;

    // Find needed subplan
    const kfft_splan_t* sP = st->primes;
    while ((sP->prime > 0) && (sP->prime != p))
        sP++;

    // Create suffled buffer
    C_CPY(x0, Fout[k]);
    for (q1 = 1, idx = 0; q1 < p; ++q1) {
        idx = sP->ridx[q1 - 1];

        k += m;
        Ftmp[idx] = Fout[k];
        C_ADDTO(Ftmp[0], Fout[k]);
    }

    // Eval recursive subplan complex FFT
    kfft_eval_cpx(sP->splan, Ftmp, Ftmp);

    // Reshuffle buffer
    k = u;

    C_CPY(Fout[k], Ftmp[0]);
    for (q1 = 1; q1 < p; ++q1) {
        idx = sP->ridx[q1 - 1];

        k += m;

        C_ADDTO(Ftmp[q1], x0);
        Fout[idx] = Ftmp[q1];
    }
}

#endif /* KFFT_RADER_ALGO */

static inline void
std_method_eval(kfft_cpx* Fout, kfft_cpx* Ftmp, const size_t fstride, const kfft_comp_t* st,
                uint32_t u, uint32_t m, uint32_t p) {
    uint32_t k = u, q1, q;
    kfft_cpx t;

    for (q1 = 0; q1 < p; ++q1) {
        Ftmp[q1] = Fout[k];
        k += m;
    }

    k = u;
    for (q1 = 0; q1 < p; ++q1) {
        uint32_t twidx = 0;
        Fout[k] = Ftmp[0];
        for (q = 1; q < p; ++q) {
            twidx += fstride * k;
            if (twidx >= st->nfft)
                twidx -= st->nfft;
            C_MUL(t, Ftmp[q], TWIDDLE(twidx, st) /*twiddles[twidx]*/);
            C_ADDTO(Fout[k], t);
        }
        k += m;
    }
}

static void
kf_bfly_generic(kfft_cpx* Fout, const size_t fstride, const kfft_comp_t* st, uint32_t m,
                uint32_t p) {
    kfft_trace("%s: m - %d | p - %d\n", "Generic FFT", m, p);

    kfft_cpx* scratch = (kfft_cpx*)KFFT_TMP_ALLOC(sizeof(kfft_cpx) * p);

    if (scratch) {

        for (uint32_t u = 0; u < m; ++u) {
#if defined(KFFT_RADER_ALGO)
            if ((!(p < KFFT_RADER_LIMIT) && (!(st->flags & KFFT_FLAG_GENERIC)))) {
                kfft_trace("[CORE] %s: %u\n", "Use Rader algorithm for resolve", p);
                rader_method_eval(Fout, scratch, fstride, st, u, m, p);
            } else {
#endif /* KFFT_RADER_ALGO */
                kfft_trace("[CORE] %s: %u\n", "Use standart algorithm for resolve", p);
                std_method_eval(Fout, scratch, fstride, st, u, m, p);
#if defined(KFFT_RADER_ALGO)
            }
#endif /* KFFT_RADER_ALGO */
        }

        KFFT_TMP_FREE(scratch);
    } else {
        kfft_trace("[LEVEL %d] %s\n", st->level, "Temporary buffer create fail.");
    }
}
