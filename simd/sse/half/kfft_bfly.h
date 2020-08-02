#pragma once
#include "kfft_simd_compat.h"
#include "kfft_simd.h"

void FUNC_SSE(kf_bfly2)(kfft_cpx* Fout, const uint32_t fstride, const kfft_plan_cpx* st,
                        uint32_t m);
void FUNC_SSE(kf_bfly4)(kfft_cpx* Fout, const uint32_t fstride, const kfft_plan_cpx* st,
                        const uint32_t m);
void FUNC_SSE(kf_bfly3)(kfft_cpx* Fout, const uint32_t fstride, const kfft_plan_cpx* st,
                        uint32_t m);
void FUNC_SSE(kf_bfly5)(kfft_cpx* Fout, const uint32_t fstride, const kfft_plan_cpx* st,
                        uint32_t m);
