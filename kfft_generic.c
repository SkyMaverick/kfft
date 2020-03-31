/* perform the butterfly for one stage of a mixed radix FFT */

#if defined(KFFT_RADER_ALGO)

static inline int
rader_method_eval(kfft_cpx* Fout, kfft_cpx* Ftmp, const size_t fstride, const kfft_comp_t* st,
                  uint32_t u, uint32_t m, uint32_t p) {

    kfft_trace("[CORE] (lvl.%d) %s\n", st->level, "Change RADER algorithm");
    uint32_t k = u, q1, idx;
    kfft_cpx x0 = {0, 0};

    // Find needed subplan
    const kfft_splan_t* sP = st->primes;
    while ((sP->prime > 0) && (sP->prime != p))
        sP++;

    // Create suffled buffers
    C_CPY(x0, Fout[k]);
    for (q1 = 1, idx = 0; q1 < p; ++q1) {
        idx = sP->qidx[q1 - 1];

        k += m;

        C_CPY(Ftmp[q1], Fout[idx]);
        C_ADDTO(Ftmp[0], Fout[k]);
    }

    kfft_part_convolution(&Ftmp[1], sP->shuffle_twiddles, sP->splan, sP->splani);
    // Reshuffle buffer
    k = u;

    C_ADDTO(Fout[k], Ftmp[0]);

    for (q1 = 1; q1 < p; ++q1) {
        idx = sP->pidx[q1 - 1];
        C_ADDTO(Ftmp[q1], x0);

        k = u + m * idx;

        C_CPY(Fout[k], Ftmp[q1]);
    }

    return 0;
}

#endif /* KFFT_RADER_ALGO */

static inline int
std_method_eval(kfft_cpx* Fout, kfft_cpx* Ftmp, const size_t fstride, const kfft_comp_t* st,
                uint32_t u, uint32_t m, uint32_t p) {
    kfft_trace("[CORE] (lvl.%d) %s\n", st->level, "Change GENERIC algorithm");

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
    return 0;
}

static void
kf_bfly_generic(kfft_cpx* Fout, const size_t fstride, const kfft_comp_t* st, uint32_t m,
                uint32_t p) {
    kfft_trace("%s: m - %d | p - %d | stride - %zu\n", "Generic FFT", m, p, fstride);

    kfft_cpx* scratch = (kfft_cpx*)KFFT_TMP_ALLOC(sizeof(kfft_cpx) * p);

    if (scratch) {

        for (uint32_t u = 0; u < m; ++u) {
#if defined(KFFT_RADER_ALGO)
            if ((p >= KFFT_RADER_LIMIT) &&
                (!((st->flags & KFFT_FLAG_GENERIC) || (st->flags & KFFT_FLAG_GENERIC_ONLY)))) {
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
