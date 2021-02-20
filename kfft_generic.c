/* perform the butterfly for one stage of a mixed radix FFT */

#if defined(KFFT_RADER_ALGO)

    #include "kfft_rader.c"

#endif /* KFFT_RADER_ALGO */

static inline kfft_return_t
std_method_eval(kfft_cpx* Fout, kfft_cpx* Ftmp, const size_t fstride, const kfft_plan_cpx* P,
                uint32_t u, uint32_t m, uint32_t p) {
    kfft_return_t ret = KFFT_RET_SUCCESS;

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
            if (twidx >= P->nfft)
                twidx -= P->nfft;
            C_MUL(t, Ftmp[q], TWIDDLE(twidx, P) /*twiddles[twidx]*/);
            C_ADDTO(Fout[k], t);
        }
        k += m;
    }
    return ret;
}

static kfft_return_t
kf_bfly_generic(kfft_cpx* Fout, const size_t fstride, const kfft_plan_cpx* plan, uint32_t m,
                uint32_t p) {
    kfft_trace_core(plan->level, "[Generic] m - %d | p - %d | stride - %zu\n", m, p, fstride);

    kfft_return_t ret = KFFT_RET_SUCCESS;

    kfft_cpx* scratch = (kfft_cpx*)KFFT_TMP_ALLOC(sizeof(kfft_cpx) * p, KFFT_PLAN_ALIGN(plan));
    if (__likely__(scratch)) {
        KFFT_ALLOCA_CLEAR(scratch, sizeof(kfft_cpx) * p);

        for (uint32_t u = 0; u < m; ++u) {
#if defined(KFFT_RADER_ALGO)
            if (__likely__((m == 1) &&                // only solid buffer
                           (p >= KFFT_RADER_LIMIT) && // over Rader limit barrier (build option)
                           (!((plan->flags & KFFT_FLAG_GENERIC) || // if NOT GENERIC flag enabled
                              (plan->flags & KFFT_FLAG_GENERIC_ONLY))) && // also ^
                           plan->prm_count)) { // if prime numbers plans created
                kfft_trace_core(plan->level, "%s: %u\n", "Use Rader algorithm for resolve", p);
                ret = rader_method_eval(Fout, scratch, fstride, plan, u, m, p);

                if (__unlikely__(ret != KFFT_RET_SUCCESS)) {
                    break;
                }
            } else {
#endif /* KFFT_RADER_ALGO */

                kfft_trace_core(plan->level, "%s: %u\n", "Use standart algorithm for resolve", p);
                ret = VEXFUNC(plan, std_method_eval, Fout, scratch, fstride, plan, u, m, p);

                if (__unlikely__(ret != KFFT_RET_SUCCESS)) {
                    break;
                }
#if defined(KFFT_RADER_ALGO)
            }
#endif /* KFFT_RADER_ALGO */
        }
        KFFT_TMP_FREE(scratch, KFFT_PLAN_ALIGN(plan));
    } else {
        ret = KFFT_RET_BUFFER_FAIL;
    }
    return ret;
}
