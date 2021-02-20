kfft_return_t
FUNC_SSE(std_method_eval)(kfft_cpx* Fout, kfft_cpx* Ftmp, const size_t fstride,
                          const kfft_plan_cpx* st, uint32_t u, uint32_t m, uint32_t p) {
    kfft_trace_core(st->level, "[GenStd (SSE2)] fstride - %zu | m - %u | p - %u\n", fstride, m, p);
    kfft_return_t ret = KFFT_RET_SUCCESS;
    uint32_t k = u, q1, q;
    __m128d T;

    for (q1 = 0; q1 < p; ++q1, k += m)
        _mm_store_pd((double*)&Ftmp[q1], _mm_load_pd((double*)&Fout[k]));

    k = u;

    for (q1 = 0; q1 < p; ++q1, k += m) {
        uint32_t twidx = 0;
        __m128d FOK = _mm_load_pd((double*)&Ftmp[0]);
        for (q = 1; q < p; ++q) {
            twidx += fstride * k;
            if (twidx >= st->nfft)
                twidx -= st->nfft;
            kfft_cpx ctw = TWIDDLE(twidx, st);
            C_MUL_SSE(T, _mm_load_pd((double*)&Ftmp[q]), _mm_load_pd((double*)&ctw));
            C_ADD_SSE(FOK, FOK, T);
            _mm_store_pd((double*)&Fout[k], FOK);
        }
    }
    return ret;
}
