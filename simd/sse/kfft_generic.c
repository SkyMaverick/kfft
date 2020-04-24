/* perform the butterfly for one stage of a mixed radix FFT */

#if defined(KFFT_RADER_ALGO)

static inline kfft_return_t
FUNC_SSE(rader_method_eval)(kfft_cpx* Fout, kfft_cpx* Ftmp, const size_t fstride,
                            const kfft_comp_t* st, uint32_t u, uint32_t m, uint32_t p) {

    KFFT_UNUSED_VAR(Fout);
    KFFT_UNUSED_VAR(Ftmp);
    KFFT_UNUSED_VAR(fstride);
    KFFT_UNUSED_VAR(st);
    KFFT_UNUSED_VAR(u);
    KFFT_UNUSED_VAR(m);
    KFFT_UNUSED_VAR(p);

    kfft_return_t ret = KFFT_RET_SUCCESS;
    //
    //    uint32_t k = u, q1, idx;
    //    kfft_cpx x0 = {0, 0};
    //
    //    // Find needed subplan
    //    const kfft_splan_t* sP = st->primes;
    //    while ((sP->prime > 0) && (sP->prime != p))
    //        sP++;
    //
    //    // Create suffled buffers
    //    C_CPY(x0, Fout[k]);
    //    for (q1 = 1, idx = 0; q1 < p; ++q1) {
    //        idx = RAD_PRIME_IDX(q1 - 1, sP);
    //
    //        k += m;
    //
    //        C_CPY(Ftmp[q1], Fout[idx]);
    //        C_ADDTO(Ftmp[0], Fout[k]);
    //    }
    //
    //    ret = kfft_part_convolution(&Ftmp[1], sP->shuffle_twiddles, sP->splan, sP->splani);
    //    if (ret == KFFT_RET_SUCCESS) {
    //        // Reshuffle buffer
    //        k = u;
    //
    //        C_ADDTO(Fout[k], Ftmp[0]);
    //
    //        for (q1 = 1; q1 < p; ++q1) {
    //            idx = RAD_INVERSE_IDX(q1 - 1, sP);
    //            C_ADDTO(Ftmp[q1], x0);
    //
    //            k = u + m * idx;
    //
    //            C_CPY(Fout[k], Ftmp[q1]);
    //        }
    //    }
    //
    return ret;
}

#endif /* KFFT_RADER_ALGO */

static inline int
FUNC_SSE(std_method_eval)(kfft_cpx* Fout, kfft_cpx* Ftmp, const size_t fstride,
                          const kfft_comp_t* st, uint32_t u, uint32_t m, uint32_t p) {
    KFFT_UNUSED_VAR(Fout);
    KFFT_UNUSED_VAR(Ftmp);
    KFFT_UNUSED_VAR(fstride);
    KFFT_UNUSED_VAR(st);
    KFFT_UNUSED_VAR(u);
    KFFT_UNUSED_VAR(m);
    KFFT_UNUSED_VAR(p);

    kfft_return_t ret = KFFT_RET_SUCCESS;
    //
    //    uint32_t k = u, q1, q;
    //    kfft_cpx t;
    //
    //    for (q1 = 0; q1 < p; ++q1) {
    //        Ftmp[q1] = Fout[k];
    //        k += m;
    //    }
    //
    //    k = u;
    //    for (q1 = 0; q1 < p; ++q1) {
    //        uint32_t twidx = 0;
    //        Fout[k] = Ftmp[0];
    //        for (q = 1; q < p; ++q) {
    //            twidx += fstride * k;
    //            if (twidx >= st->nfft)
    //                twidx -= st->nfft;
    //            C_MUL(t, Ftmp[q], TWIDDLE(twidx, st) /*twiddles[twidx]*/);
    //            C_ADDTO(Fout[k], t);
    //        }
    //        k += m;
    //    }
    return ret;
}

static kfft_return_t
FUNC_SSE(kf_bfly_generic)(kfft_cpx* Fout, const size_t fstride, const kfft_comp_t* st, uint32_t m,
                          uint32_t p) {

    kfft_trace_core(st->level, "[Generic (SSE)] m - %d | p - %d | stride - %zu\n", m, p, fstride);

    KFFT_UNUSED_VAR(Fout);
    KFFT_UNUSED_VAR(fstride);
    KFFT_UNUSED_VAR(st);
    KFFT_UNUSED_VAR(m);
    KFFT_UNUSED_VAR(p);

    kfft_return_t ret = KFFT_RET_SUCCESS;
    //
    //    kfft_cpx* scratch = (kfft_cpx*)KFFT_TMP_ALLOC(sizeof(kfft_cpx) * p);
    //    if (scratch) {
    //        KFFT_ALLOCA_CLEAR(scratch, sizeof(kfft_cpx) * p);
    //
    //        for (uint32_t u = 0; u < m; ++u) {
    //#if defined(KFFT_RADER_ALGO)
    //            if ((p >= KFFT_RADER_LIMIT) &&
    //                (!((st->flags & KFFT_FLAG_GENERIC) || (st->flags & KFFT_FLAG_GENERIC_ONLY)))
    //                && st->prm_count) {
    //
    //                kfft_trace_core(st->level, "%s: %u\n", "Use Rader algorithm for resolve", p);
    //                ret = rader_method_eval(Fout, scratch, fstride, st, u, m, p);
    //
    //                if (ret != KFFT_RET_SUCCESS) {
    //                    break;
    //                }
    //            } else {
    //#endif /* KFFT_RADER_ALGO */
    //
    //                kfft_trace_core(st->level, "%s: %u\n", "Use standart algorithm for resolve",
    //                p); ret = std_method_eval(Fout, scratch, fstride, st, u, m, p);
    //
    //                if (ret != KFFT_RET_SUCCESS) {
    //                    break;
    //                }
    //#if defined(KFFT_RADER_ALGO)
    //            }
    //#endif /* KFFT_RADER_ALGO */
    //        }
    //        KFFT_TMP_FREE(scratch);
    //    } else {
    //        ret = KFFT_RET_BUFFER_FAIL;
    //    }
    return ret;
}