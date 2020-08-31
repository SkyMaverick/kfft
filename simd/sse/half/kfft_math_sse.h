#pragma once

#include "kfft_types.h"

void FUNC_SSE(kfft_math_hadamard_cpx)(kfft_cpx* Fout, const kfft_cpx* Fin, uint32_t size);
