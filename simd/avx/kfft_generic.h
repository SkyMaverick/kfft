#pragma once
#include "kfft_simd_compat.h"
#include "kfft_simd.h"

kfft_return_t FUNC_AVX(std_method_eval)(kfft_cpx* Fout, kfft_cpx* Ftmp, const size_t fstride,
                                        const kfft_comp_t* st, uint32_t u, uint32_t m, uint32_t p);
