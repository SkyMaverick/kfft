#pragma once

#include "kfft_cpx.h"

KFFT_API int
kfft_convolution (kfft_cpx* Fout, kfft_cpx* Fin, kfft_comp_t* P);
