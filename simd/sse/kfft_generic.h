#pragma once
#include "kfft_simd_compat.h"

kfft_return_t FUNC_SSE(rader_method_eval)(kfft_cpx* Fout, kfft_cpx* Ftmp, const size_t fstride,
                                          const kfft_comp_t* st, uint32_t u, uint32_t m,
                                          uint32_t p);

int FUNC_SSE(std_method_eval)(kfft_cpx* Fout, kfft_cpx* Ftmp, const size_t fstride,
                              const kfft_comp_t* st, uint32_t u, uint32_t m, uint32_t p);
