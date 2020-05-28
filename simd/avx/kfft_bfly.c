#include "kfft_simd.h"
#include "kfft_math_intern.h"

#define M2CH(X) (*(kfft_cpx*)(&((X)[0])))
#define M2CL(X) (*(kfft_cpx*)(&((X)[2])))

void
FUNC_AVX(kf_bfly2)(kfft_cpx* Fout, const uint32_t fstride, const kfft_comp_t* st, uint32_t m) {
    kfft_trace_core(st->level, "[BFLY2 (AVX)] fstride - %u | m - %u\n", fstride, m);

    kfft_cpx* Fout2 = Fout + m;
    uint32_t twidx = 0;

    do {
        __m256d FT, Tw;
        __m256d F0 = _mm256_setr_pd(Fout->r, Fout->i, Fout->r, Fout->i);
        __m256d F2 = _mm256_setr_pd(Fout2->r, Fout2->i, Fout2->r, Fout2->i);

        kfft_cpx tw = TWIDDLE(twidx, st);
        Tw = _mm256_setr_pd(tw.r, tw.i, tw.r, tw.i);

        C_MUL_AVX1(FT, F2, Tw);
        FT = _mm256_mul_pd(FT, _mm256_setr_pd(1.0, 1.0, -1.0, -1.0));
        F0 = _mm256_add_pd(F0, FT);

        C_CPY(*Fout, M2CH(F0));
        C_CPY(*Fout2, M2CL(F0));

        twidx += fstride;
        ++Fout2;
        ++Fout;
    } while (--m);
}

void
FUNC_AVX(kf_bfly4)(kfft_cpx* Fout, const uint32_t fstride, const kfft_comp_t* st,
                   const uint32_t m) {
    kfft_trace_core(st->level, "[BFLY4 (AVX)] fstride - %u | m - %u\n", fstride, m);
}

void
FUNC_AVX(kf_bfly3)(kfft_cpx* Fout, const uint32_t fstride, const kfft_comp_t* st, uint32_t m) {
    kfft_trace_core(st->level, "[BFLY3 (AVX)] fstride - %u | m - %u\n", fstride, m);
}

void
FUNC_AVX(kf_bfly5)(kfft_cpx* Fout, const uint32_t fstride, const kfft_comp_t* st, uint32_t m) {
    kfft_trace_core(st->level, "[BFLY5 (AVX)] fstride - %u | m - %u\n", fstride, m);
}
