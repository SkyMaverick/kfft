#pragma once
#include "kfft_simd_compat.h"
#include "kfft_simd.h"

kfft_return_t FUNC_SSE(std_method_eval)(kfft_cpx* Fout, kfft_cpx* Ftmp, const size_t fstride,
                                        const kfft_plan_cpx* st, uint32_t u, uint32_t m,
                                        uint32_t p);
