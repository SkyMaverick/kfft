#pragma once
#include "kfft_simd_compat.h"
#include "kfft_simd.h"

#define FUNC_SSE(X) X##_sse

void FUNC_SSE(kf_bfly2)(kfft_cpx* Fout, const uint32_t fstride, const kfft_plan_cpx* st,
                        uint32_t m);
void FUNC_SSE(kf_bfly4)(kfft_cpx* Fout, const uint32_t fstride, const kfft_plan_cpx* st,
                        const uint32_t m);
void FUNC_SSE(kf_bfly3)(kfft_cpx* Fout, const uint32_t fstride, const kfft_plan_cpx* st,
                        uint32_t m);
void FUNC_SSE(kf_bfly5)(kfft_cpx* Fout, const uint32_t fstride, const kfft_plan_cpx* st,
                        uint32_t m);

kfft_return_t FUNC_SSE(std_method_eval)(kfft_cpx* Fout, kfft_cpx* Ftmp, const size_t fstride,
                                        const kfft_plan_cpx* st, uint32_t u, uint32_t m,
                                        uint32_t p);

void FUNC_SSE(kfft_math_hadamard_cpx)(kfft_cpx* Fout, const kfft_cpx* Fin, uint32_t size);
