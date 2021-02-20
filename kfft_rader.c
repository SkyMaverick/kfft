static inline kfft_return_t
kfft_part_convolution(kfft_cpx* Fout, kfft_cpx* Fin, kfft_plan_cpx* P, kfft_plan_cpx* Pi) {
    kfft_cpx* Fbuf = KFFT_TMP_ALLOC(P->nfft * sizeof(kfft_cpx), KFFT_PLAN_ALIGN(P));
    if (__likely__(Fbuf)) {
        kfft_return_t ret = kfft_eval_cpx(P, Fout, Fbuf);
        if (__likely__(ret == KFFT_RET_SUCCESS)) {
            VEXFUNC(P, kfft_math_hadamard_cpx, Fbuf, Fin, P->nfft);
            ret = kfft_eval_cpx(Pi, Fbuf, Fout);
        }

        KFFT_TMP_FREE(Fbuf, KFFT_PLAN_ALIGN(P));
        return ret;
    }
    return KFFT_RET_BUFFER_FAIL;
}

static inline kfft_return_t
rader_method_eval(kfft_cpx* Fout, kfft_cpx* Ftmp, const size_t fstride, const kfft_plan_cpx* P,
                  uint32_t u, uint32_t m, uint32_t p) {
    (void)fstride; // disable unused parameter
    kfft_return_t ret = KFFT_RET_SUCCESS;

    uint32_t k = u, q1, idx;
    kfft_cpx x0 = {0, 0};

    // Find needed subplan
    const kfft_plan_rader* sP = P->primes;
    while ((sP->prime > 0) && (sP->prime != p))
        sP++;

    // Create suffled buffers
    C_CPY(x0, Fout[k]);
    // MSVC non-zero alloc workaround
    Ftmp[0].r = 0;
    Ftmp[0].i = 0;
    for (q1 = 1, idx = 0; q1 < p; ++q1) {
        idx = RAD_PRIME_IDX(q1 - 1, sP);

        k += m;

        C_CPY(Ftmp[q1], Fout[idx]);
        C_ADDTO(Ftmp[0], Fout[k]);
    }

    ret = kfft_part_convolution(&Ftmp[1], sP->shuffle_twiddles, sP->plan, sP->plan_inv);
    if (__likely__(ret == KFFT_RET_SUCCESS)) {
        // Reshuffle buffer
        k = u;

        C_ADDTO(Fout[k], Ftmp[0]);

        for (q1 = 1; q1 < p; ++q1) {
            idx = RAD_INVERSE_IDX(q1 - 1, sP);
            C_ADDTO(Ftmp[q1], x0);

            k = u + m * idx;

            C_CPY(Fout[k], Ftmp[q1]);
        }
    }

    return ret;
}
