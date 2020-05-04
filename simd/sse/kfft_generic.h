#pragma once
#include "kfft_simd_compat.h"
#include "kfft_simd.h"

kfft_return_t FUNC_SSE(kf_bfly_generic)(kfft_cpx* Fout, const size_t fstride, const kfft_comp_t* st,
                                        uint32_t m, uint32_t p);
