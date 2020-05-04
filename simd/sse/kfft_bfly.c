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
        __m128d mtw = _mm_loaddup_pd(&(tw.r));
        t = _mm_loaddup_pd(&(tw.i));
        C_MULDUP_SSE(t, mf2, mtw);
#else
        __m128d mtw = _mm_load_pd((double*)&tw);
        C_MUL_SSE(t, mf2, mtw);
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
        __m128d T4i, invi = {1, -1}, invr = {-1, 1};
        if (st->flags & KFFT_FLAG_INVERSE) {
            T4i = _mm_mul_pd(T4, invi);
            F0 = _mm_add_pd(T5, _mm_shuffle_pd(T4i, T4i, 0x1));
            F2 = _mm_add_pd(_mm_shuffle_pd(T5, T5, 0x1), _mm_mul_pd(T4, invr));
            F2 = _mm_shuffle_pd(F2, F2, 0x1);
        } else {
            F0 = _mm_add_pd(_mm_shuffle_pd(T5, T5, 0x1), _mm_mul_pd(T4, invr));
            F0 = _mm_shuffle_pd(F0, F0, 0x1);
            T4i = _mm_mul_pd(T4, invi);
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
    KFFT_UNUSED_VAR(Fout);
    KFFT_UNUSED_VAR(fstride);
    KFFT_UNUSED_VAR(st);
    KFFT_UNUSED_VAR(m);
    //    uint32_t k = m;
    //    const uint32_t m2 = 2 * m;
    //    uint32_t tw1, tw2;
    //    kfft_cpx scratch[5];
    //    kfft_cpx epi3;
    //    epi3 = TWIDDLE(fstride * m, st);
    //
    //    tw1 = tw2 = 0;
    //
    //    do {
    //        C_MUL(scratch[1], Fout[m], TWIDDLE(tw1, st));
    //        C_MUL(scratch[2], Fout[m2], TWIDDLE(tw2, st));
    //
    //        C_ADD(scratch[3], scratch[1], scratch[2]);
    //        C_SUB(scratch[0], scratch[1], scratch[2]);
    //        tw1 += fstride;
    //        tw2 += fstride * 2;
    //
    //        Fout[m].r = Fout->r - HALF_OF(scratch[3].r);
    //        Fout[m].i = Fout->i - HALF_OF(scratch[3].i);
    //
    //        C_MULBYSCALAR(scratch[0], epi3.i);
    //
    //        C_ADDTO(*Fout, scratch[3]);
    //
    //        Fout[m2].r = Fout[m].r + scratch[0].i;
    //        Fout[m2].i = Fout[m].i - scratch[0].r;
    //
    //        Fout[m].r -= scratch[0].i;
    //        Fout[m].i += scratch[0].r;
    //
    //        ++Fout;
    //    } while (--k);
}

void
FUNC_SSE(kf_bfly5)(kfft_cpx* Fout, const uint32_t fstride, const kfft_comp_t* st, uint32_t m) {
    kfft_trace_core(st->level, "[BFLY5 (SSE)] fstride - %u | m - %u\n", fstride, m);
    KFFT_UNUSED_VAR(Fout);
    KFFT_UNUSED_VAR(fstride);
    KFFT_UNUSED_VAR(st);
    KFFT_UNUSED_VAR(m);
    //    kfft_cpx *Fout0, *Fout1, *Fout2, *Fout3, *Fout4;
    //    uint32_t u;
    //    kfft_cpx scratch[13];
    //    kfft_cpx ya, yb;
    //    ya = TWIDDLE(fstride * m, st);
    //    yb = TWIDDLE(fstride * 2 * m, st);
    //
    //    Fout0 = Fout;
    //    Fout1 = Fout0 + m;
    //    Fout2 = Fout0 + 2 * m;
    //    Fout3 = Fout0 + 3 * m;
    //    Fout4 = Fout0 + 4 * m;
    //
    //    for (u = 0; u < m; ++u) {
    //        scratch[0] = *Fout0;
    //
    //        C_MUL(scratch[1], *Fout1, TWIDDLE(u * fstride, st));
    //        C_MUL(scratch[2], *Fout2, TWIDDLE(2 * u * fstride, st));
    //        C_MUL(scratch[3], *Fout3, TWIDDLE(3 * u * fstride, st));
    //        C_MUL(scratch[4], *Fout4, TWIDDLE(4 * u * fstride, st));
    //
    //        C_ADD(scratch[7], scratch[1], scratch[4]);
    //        C_SUB(scratch[10], scratch[1], scratch[4]);
    //        C_ADD(scratch[8], scratch[2], scratch[3]);
    //        C_SUB(scratch[9], scratch[2], scratch[3]);
    //
    //        Fout0->r += scratch[7].r + scratch[8].r;
    //        Fout0->i += scratch[7].i + scratch[8].i;
    //
    //        scratch[5].r = scratch[0].r + S_MUL(scratch[7].r, ya.r) + S_MUL(scratch[8].r, yb.r);
    //        scratch[5].i = scratch[0].i + S_MUL(scratch[7].i, ya.r) + S_MUL(scratch[8].i, yb.r);
    //
    //        scratch[6].r = S_MUL(scratch[10].i, ya.i) + S_MUL(scratch[9].i, yb.i);
    //        scratch[6].i = -S_MUL(scratch[10].r, ya.i) - S_MUL(scratch[9].r, yb.i);
    //
    //        C_SUB(*Fout1, scratch[5], scratch[6]);
    //        C_ADD(*Fout4, scratch[5], scratch[6]);
    //
    //        scratch[11].r = scratch[0].r + S_MUL(scratch[7].r, yb.r) + S_MUL(scratch[8].r, ya.r);
    //        scratch[11].i = scratch[0].i + S_MUL(scratch[7].i, yb.r) + S_MUL(scratch[8].i, ya.r);
    //        scratch[12].r = -S_MUL(scratch[10].i, yb.i) + S_MUL(scratch[9].i, ya.i);
    //        scratch[12].i = S_MUL(scratch[10].r, yb.i) - S_MUL(scratch[9].r, ya.i);
    //
    //        C_ADD(*Fout2, scratch[11], scratch[12]);
    //        C_SUB(*Fout3, scratch[11], scratch[12]);
    //
    //        ++Fout0;
    //        ++Fout1;
    //        ++Fout2;
    //        ++Fout3;
    //        ++Fout4;
    //  }
}
