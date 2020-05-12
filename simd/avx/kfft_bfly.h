#pragma once
#include "kfft_simd_compat.h"
#include "kfft_simd.h"

void FUNC_AVX(kf_bfly2)(kfft_cpx* Fout, const uint32_t fstride, const kfft_comp_t* st, uint32_t m);
void FUNC_AVX(kf_bfly4)(kfft_cpx* Fout, const uint32_t fstride, const kfft_comp_t* st,
                        const uint32_t m);
void FUNC_AVX(kf_bfly3)(kfft_cpx* Fout, const uint32_t fstride, const kfft_comp_t* st, uint32_t m);
void FUNC_AVX(kf_bfly5)(kfft_cpx* Fout, const uint32_t fstride, const kfft_comp_t* st, uint32_t m);
