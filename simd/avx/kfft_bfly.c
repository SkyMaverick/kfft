#include "kfft_simd.h"
#include "kfft_math_intern.h"

void
FUNC_AVX(kf_bfly2)(kfft_cpx* Fout, const uint32_t fstride, const kfft_comp_t* st, uint32_t m) {
    kfft_trace_core(st->level, "[BFLY2 (AVX)] fstride - %u | m - %u\n", fstride, m);
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
