#include "kfft_simd.h"
#include "kfft_math_intern.h"

void
FUNC_SSE(kf_bfly2)(kfft_cpx* Fout, const uint32_t fstride, const kfft_comp_t* st, uint32_t m) {
    kfft_trace_core(st->level, "[BFLY2 (SSE)] fstride - %u | m - %u\n", fstride, m);

    kfft_cpx* Fout2 = Fout + m;
    uint32_t twidx = 0;
    do {
        __m128d t;
        __m128d mf = _mm_load_pd((double*)Fout);
        __m128d mf2 = _mm_load_pd((double*)Fout2);

        kfft_cpx tw = TWIDDLE(twidx, st);
#if defined(KFFT_HAVE_SSE3)
        t = _mm_loaddup_pd(&(tw.i));
        C_MULDUP_SSE(t, mf2, _mm_loaddup_pd(&(tw.r)));
#else
        C_MUL_SSE(t, mf2, _mm_load_pd((double*)&tw));
#endif /* HAVE_SSE3 */

        C_SUB_SSE(mf2, mf, t);
        C_ADD_SSE(mf, mf, t);

        _mm_store_pd((double*)Fout2, mf2);
        _mm_store_pd((double*)Fout, mf);

        twidx += fstride;
        ++Fout2;
        ++Fout;
    } while (--m);
}

void
FUNC_SSE(kf_bfly4)(kfft_cpx* Fout, const uint32_t fstride, const kfft_comp_t* st,
                   const uint32_t m) {
    kfft_trace_core(st->level, "[BFLY4 (SSE)] fstride - %u | m - %u\n", fstride, m);

    uint32_t tw1, tw2, tw3;
    tw1 = tw2 = tw3 = 0;

    uint32_t k = m;
    const uint32_t m2 = 2 * m;
    const uint32_t m3 = 3 * m;

    __m128d T0, T1, T2, T3, T4, T5;
    do {
        kfft_cpx ctw0 = TWIDDLE(tw1, st);
        kfft_cpx ctw1 = TWIDDLE(tw2, st);
        kfft_cpx ctw2 = TWIDDLE(tw3, st);

        __m128d F0, F1, F2;
        __m128d mfout = _mm_load_pd((double*)(Fout));
        F0 = _mm_load_pd((double*)(&Fout[m]));
        F1 = _mm_load_pd((double*)(&Fout[m2]));
        F2 = _mm_load_pd((double*)(&Fout[m3]));
#if defined(KFFT_HAVE_SSE3)
        T0 = _mm_loaddup_pd(&(ctw0.i));
        T1 = _mm_loaddup_pd(&(ctw1.i));
        T2 = _mm_loaddup_pd(&(ctw2.i));

        C_MULDUP_SSE(T0, F0, _mm_loaddup_pd(&(ctw0.r)));
        C_MULDUP_SSE(T1, F1, _mm_loaddup_pd(&(ctw1.r)));
        C_MULDUP_SSE(T2, F2, _mm_loaddup_pd(&(ctw2.r)));
#else  /* KFFT_HAVE_SSE3 */
        C_MUL_SSE(T0, F0, _mm_load_pd((double*)&ctw0));
        C_MUL_SSE(T1, F1, _mm_load_pd((double*)&ctw1));
        C_MUL_SSE(T2, F2, _mm_load_pd((double*)&ctw2));
#endif /* KFFT_HAVE_SSE3 */
        C_SUB_SSE(T5, mfout, T1);
        C_ADD_SSE(T3, T0, T2);
        C_ADD_SSE(mfout, mfout, T1);
        C_SUB_SSE(T4, T0, T2);
        C_SUB_SSE(F1, mfout, T3);
        C_ADD_SSE(mfout, mfout, T3);
#if defined(KFFT_HAVE_SSE3)
        if (st->flags & KFFT_FLAG_INVERSE) {
            F0 = _mm_addsub_pd(T5, _mm_shuffle_pd(T4, T4, 0x1));
            F2 = _mm_addsub_pd(_mm_shuffle_pd(T5, T5, 0x1), T4);
            F2 = _mm_shuffle_pd(F2, F2, 0x1);
        } else {
            F0 = _mm_addsub_pd(_mm_shuffle_pd(T5, T5, 0x1), T4);
            F0 = _mm_shuffle_pd(F0, F0, 0x1);
            F2 = _mm_addsub_pd(T5, _mm_shuffle_pd(T4, T4, 0x1));
        } /* KFFT_FLAG_INVERSE */
#else     /* KFFT_HAVE_SSE3 */
        __m128d T4i;
        if (st->flags & KFFT_FLAG_INVERSE) {
            T4i = _mm_mul_pd(T4, _mm_set_pd(-1, 1));
            F0 = _mm_add_pd(T5, _mm_shuffle_pd(T4i, T4i, 0x1));
            F2 = _mm_add_pd(_mm_shuffle_pd(T5, T5, 0x1), _mm_mul_pd(T4, _mm_set_pd(1, -1)));
            F2 = _mm_shuffle_pd(F2, F2, 0x1);
        } else {
            F0 = _mm_add_pd(_mm_shuffle_pd(T5, T5, 0x1), _mm_mul_pd(T4, _mm_set_pd(1, -1)));
            F0 = _mm_shuffle_pd(F0, F0, 0x1);
            T4i = _mm_mul_pd(T4, _mm_set_pd(-1, 1));
            F2 = _mm_add_pd(T5, _mm_shuffle_pd(T4i, T4i, 0x1));
        } /* KFFT_FLAG_INVERSE */
#endif    /* KFFT_HAVE_SSE3 */
        _mm_store_pd((double*)Fout, mfout);
        _mm_store_pd((double*)(&Fout[m]), F0);
        _mm_store_pd((double*)(&Fout[m2]), F1);
        _mm_store_pd((double*)(&Fout[m3]), F2);

        tw1 += fstride;
        tw2 += fstride * 2;
        tw3 += fstride * 3;

        ++Fout;
    } while (--k);
}

void
FUNC_SSE(kf_bfly3)(kfft_cpx* Fout, const uint32_t fstride, const kfft_comp_t* st, uint32_t m) {
    kfft_trace_core(st->level, "[BFLY3 (SSE)] fstride - %u | m - %u\n", fstride, m);

    uint32_t k = m;
    const uint32_t m2 = 2 * m;
    uint32_t tw1, tw2;
    tw1 = tw2 = 0;

    kfft_cpx epi3 = TWIDDLE(fstride * m, st);
    __m128d mepi3 = _mm_set1_pd(epi3.i);

    __m128d T0, T1, T2, T3;
    do {
        kfft_cpx ctw0 = TWIDDLE(tw1, st);
        kfft_cpx ctw1 = TWIDDLE(tw2, st);

        __m128d F0, F1;
        __m128d mfout = _mm_load_pd((double*)(Fout));
        F0 = _mm_load_pd((double*)(&Fout[m]));
        F1 = _mm_load_pd((double*)(&Fout[m2]));
#if defined(KFFT_HAVE_SSE3)
        T1 = _mm_loaddup_pd(&(ctw0.i));
        T2 = _mm_loaddup_pd(&(ctw1.i));

        C_MULDUP_SSE(T1, F0, _mm_loaddup_pd(&(ctw0.r)));
        C_MULDUP_SSE(T2, F1, _mm_loaddup_pd(&(ctw1.r)));
#else  /* KFFT_HAVE_SSE3 */
        C_MUL_SSE(T1, F0, _mm_load_pd((double*)&ctw0));
        C_MUL_SSE(T2, F1, _mm_load_pd((double*)&ctw1));
#endif /* KFFT_HAVE_SSE3 */
        C_ADD_SSE(T3, T1, T2);
        C_SUB_SSE(T0, T1, T2);

        C_SUB_SSE(F0, mfout, _mm_mul_pd(T3, _mm_set1_pd(0.5)));
        T0 = _mm_mul_pd(T0, mepi3);
        C_ADD_SSE(mfout, mfout, T3);

#if defined(KFFT_HAVE_SSE3)
        F1 = _mm_addsub_pd(_mm_shuffle_pd(F0, F0, 0x1), T0);
        F1 = _mm_shuffle_pd(F1, F1, 0x1);

        F0 = _mm_addsub_pd(F0, _mm_shuffle_pd(T0, T0, 0x1));
#else  /* KFFT_HAVE_SSE3 */
        __m128d tmp = _mm_mul_pd(T0, _mm_set_pd(1, -1));
        F1 = _mm_add_pd(F0, _mm_shuffle_pd(tmp, tmp, 0x1));

        tmp = _mm_mul_pd(T0, _mm_set_pd(-1, 1));
        F0 = _mm_add_pd(F0, _mm_shuffle_pd(tmp, tmp, 0x1));
#endif /* KFFT_HAVE_SSE3 */
        _mm_store_pd((double*)Fout, mfout);
        _mm_store_pd((double*)(&Fout[m]), F0);
        _mm_store_pd((double*)(&Fout[m2]), F1);

        tw1 += fstride;
        tw2 += fstride * 2;
        ++Fout;
    } while (--k);
}

void
FUNC_SSE(kf_bfly5)(kfft_cpx* Fout, const uint32_t fstride, const kfft_comp_t* st, uint32_t m) {
    kfft_trace_core(st->level, "[BFLY5 (SSE)] fstride - %u | m - %u\n", fstride, m);

    uint32_t u;

    __m128d T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12;

    kfft_cpx ya = TWIDDLE(fstride * m, st);
    __m128d YAR = _mm_set1_pd(ya.r);
    __m128d YAI = _mm_set1_pd(ya.i);
    kfft_cpx yb = TWIDDLE(fstride * 2 * m, st);
    __m128d YBR = _mm_set1_pd(yb.r);
    __m128d YBI = _mm_set1_pd(yb.i);

    kfft_cpx* Fout0 = Fout;
    kfft_cpx* Fout1 = Fout0 + m;
    kfft_cpx* Fout2 = Fout0 + 2 * m;
    kfft_cpx* Fout3 = Fout0 + 3 * m;
    kfft_cpx* Fout4 = Fout0 + 4 * m;

    for (u = 0; u < m; ++u) {

        kfft_cpx ctw1 = TWIDDLE(u * fstride, st);
        kfft_cpx ctw2 = TWIDDLE(2 * u * fstride, st);
        kfft_cpx ctw3 = TWIDDLE(3 * u * fstride, st);
        kfft_cpx ctw4 = TWIDDLE(4 * u * fstride, st);

        __m128d F0 = _mm_load_pd((double*)Fout0);
        __m128d F1 = _mm_load_pd((double*)Fout1);
        __m128d F2 = _mm_load_pd((double*)Fout2);
        __m128d F3 = _mm_load_pd((double*)Fout3);
        __m128d F4 = _mm_load_pd((double*)Fout4);

        T0 = F0;
#if defined(KFFT_HAVE_SSE3)
        T1 = _mm_loaddup_pd(&(ctw1.i));
        T2 = _mm_loaddup_pd(&(ctw2.i));
        T3 = _mm_loaddup_pd(&(ctw3.i));
        T4 = _mm_loaddup_pd(&(ctw4.i));

        C_MULDUP_SSE(T1, F1, _mm_loaddup_pd(&(ctw1.r)));
        C_MULDUP_SSE(T2, F2, _mm_loaddup_pd(&(ctw2.r)));
        C_MULDUP_SSE(T3, F3, _mm_loaddup_pd(&(ctw3.r)));
        C_MULDUP_SSE(T4, F4, _mm_loaddup_pd(&(ctw4.r)));
#else  /* KFFT_HAVE_SSE3 */
        C_MUL_SSE(T1, F1, _mm_load_pd((double*)&ctw1));
        C_MUL_SSE(T2, F2, _mm_load_pd((double*)&ctw2));
        C_MUL_SSE(T3, F3, _mm_load_pd((double*)&ctw3));
        C_MUL_SSE(T4, F4, _mm_load_pd((double*)&ctw4));
#endif /* KFFT_HAVE_SSE3 */

        C_ADD_SSE(T7, T1, T4);
        C_SUB_SSE(T10, T1, T4);
        C_ADD_SSE(T8, T2, T3);
        C_SUB_SSE(T9, T2, T3);

        C_ADD_SSE(F0, F0, _mm_add_pd(T7, T8));

        C_ADD_SSE(T5, T0, _mm_add_pd(_mm_mul_pd(T7, YAR), _mm_mul_pd(T8, YBR)));
        C_ADD_SSE(T6, _mm_mul_pd(T10, YAI), _mm_mul_pd(T9, YBI));
        T6 = _mm_mul_pd(_mm_shuffle_pd(T6, T6, 0x1), _mm_set_pd(-1, 1));

        C_SUB_SSE(F1, T5, T6);
        C_ADD_SSE(F4, T5, T6);

        C_ADD_SSE(T11, T0, _mm_add_pd(_mm_mul_pd(T7, YBR), _mm_mul_pd(T8, YAR)));

        C_SUB_SSE(T12, _mm_mul_pd(_mm_shuffle_pd(T9, T10, 0x1), _mm_shuffle_pd(YAI, YBI, 0x1)),
                  _mm_mul_pd(_mm_shuffle_pd(T10, T9, 0x1), _mm_shuffle_pd(YBI, YAI, 0x1)));

        C_ADD_SSE(F2, T11, T12);
        C_SUB_SSE(F3, T11, T12);

        _mm_store_pd((double*)Fout0, F0);
        _mm_store_pd((double*)Fout1, F1);
        _mm_store_pd((double*)Fout2, F2);
        _mm_store_pd((double*)Fout3, F3);
        _mm_store_pd((double*)Fout4, F4);

        ++Fout0;
        ++Fout1;
        ++Fout2;
        ++Fout3;
        ++Fout4;
    }
}
