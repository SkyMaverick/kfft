#include "kfft_simd.h"
#include "kfft_math_intern.h"

#define M2CH(X) (*(kfft_cpx*)(&((X)[0])))
#define M2CL(X) (*(kfft_cpx*)(&((X)[2])))

#define C2M128(X) (*((__m128d*)(&(X))))

#define INVHI_AVX(X) _mm256_mul_pd((X), _mm256_setr_pd(-1.0, -1.0, 1.0, 1.0))
#define INVLO_AVX(X) _mm256_mul_pd((X), _mm256_setr_pd(1.0, 1.0, -1.0, -1.0))
#define INVHL_AVX(X) _mm256_mul_pd((X), _mm256_setr_pd(-1.0, -1.0, -1.0, -1.0))

void
FUNC_AVX(kf_bfly2)(kfft_cpx* Fout, const uint32_t fstride, const kfft_comp_t* st, uint32_t m) {
    kfft_trace_core(st->level, "[BFLY2 (AVX)] fstride - %u | m - %u\n", fstride, m);
    //
    //    kfft_cpx* Fout2 = Fout + m;
    //    uint32_t twidx = 0;
    //
    //    do {
    //        __m256d FT, Tw;
    //        __m256d F0 = _mm256_setr_pd(Fout->r, Fout->i, Fout->r, Fout->i);
    //        __m256d F2 = _mm256_setr_pd(Fout2->r, Fout2->i, Fout2->r, Fout2->i);
    //
    //        kfft_cpx tw = TWIDDLE(twidx, st);
    //        Tw = _mm256_setr_pd(tw.r, tw.i, tw.r, tw.i);
    //
    //        C_MUL_AVX1(FT, F2, Tw);
    //        FT = INVLO_AVX(FT);
    //        F0 = _mm256_add_pd(F0, FT);
    //
    //        C_CPY(*Fout, M2CH(F0));
    //        C_CPY(*Fout2, M2CL(F0));
    //
    //        twidx += fstride;
    //        ++Fout2;
    //        ++Fout;
    //    } while (--m);
}

void
FUNC_AVX(kf_bfly4)(kfft_cpx* Fout, const uint32_t fstride, const kfft_comp_t* st,
                   const uint32_t m) {
    kfft_trace_core(st->level, "[BFLY4 (AVX)] fstride - %u | m - %u\n", fstride, m);

    //    uint32_t tw1, tw2, tw3;
    //    tw1 = tw2 = tw3 = 0;
    //
    //    uint32_t k = m;
    //    const uint32_t m2 = 2 * m;
    //    const uint32_t m3 = 3 * m;
    //
    //    __m256d T0, T1, T2, T3, T4, T5;
    //    do {
    //        kfft_cpx ctw0 = TWIDDLE(tw1, st);
    //        kfft_cpx ctw1 = TWIDDLE(tw2, st);
    //        kfft_cpx ctw2 = TWIDDLE(tw3, st);
    //
    //        _mm256_zeroupper();
    //        __m128d fo = _mm_load_pd((double*)(Fout));
    //        __m128d f0 = _mm_load_pd((double*)(&Fout[m]));
    //        __m128d f1 = _mm_load_pd((double*)(&Fout[m2]));
    //        __m128d f2 = _mm_load_pd((double*)(&Fout[m3]));
    //
    //
    //        __m256d Fo = _mm256_broadcast_pd(&fo);
    //        __m256d F0 = _mm256_setr_m128d(f0,C2M128(ctw0));
    //        __m256d F1 = _mm256_setr_m128d(f1,C2M128(ctw1));
    //        __m256d F2 = _mm256_setr_m128d(f2,C2M128(ctw2));
    //
    //        C_MUL_AVX2 (T0, F0);
    //        C_MUL_AVX2 (T1, F1);
    //        C_MUL_AVX2 (T2, F2);
    //
    //        T5 = _mm256_sub_pd(Fo, T1);
    //
    //        __m256d T34 = _mm256_add_pd(T0, INVLO_AVX(T2));
    //        T4 = _mm256_permute_pd(_mm256_permute2f128_pd(T34,T34,0x11),0x15);
    //
    //        Fo = _mm256_add_pd(Fo, T1);
    //        Fo = _mm256_add_pd(Fo, INVHI_AVX(_mm256_permute2f128_pd(T34, T34, 0x2)));
    //
    //        if (st->flags & KFFT_FLAG_INVERSE) {
    //            F0 = _mm256_add_pd(T5, _mm256_mul_pd(T4, _mm256_setr_pd(-1.0,1.0,1.0,-1.0)));
    //        } else {
    //            F0 = _mm256_add_pd(T5, _mm256_mul_pd(T4, _mm256_setr_pd(1.0,-1.0,-1.0,1.0)));
    //        } /* KFFT_FLAG_INVERSE */
    //
    //        C_CPY(*Fout, M2CL(Fo));
    //        C_CPY(Fout[m2], M2CH(Fo));
    //        C_CPY(Fout[m], M2CH(F0));
    //        C_CPY(Fout[m3], M2CL(F0));
    //
    ////        _mm256_zeroupper();
    ////        _mm_store_pd((double*)Fout,        _mm256_castpd256_pd128(Fo));
    ////        _mm_store_pd((double*)(&Fout[m]),  _mm256_castpd256_pd128(F0));
    ////        _mm_store_pd((double*)(&Fout[m2]), _mm256_castpd256_pd128(F1));
    ////        _mm_store_pd((double*)(&Fout[m3]), _mm256_castpd256_pd128(F2));
    //
    //        tw1 += fstride;
    //        tw2 += fstride * 2;
    //        tw3 += fstride * 3;
    //
    //        ++Fout;
}
while (--k)
    ;
}

void
FUNC_AVX(kf_bfly3)(kfft_cpx* Fout, const uint32_t fstride, const kfft_comp_t* st, uint32_t m) {
    kfft_trace_core(st->level, "[BFLY3 (AVX)] fstride - %u | m - %u\n", fstride, m);
}

void
FUNC_AVX(kf_bfly5)(kfft_cpx* Fout, const uint32_t fstride, const kfft_comp_t* st, uint32_t m) {
    kfft_trace_core(st->level, "[BFLY5 (AVX)] fstride - %u | m - %u\n", fstride, m);
}
