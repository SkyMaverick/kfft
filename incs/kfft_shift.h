#pragma once
#include "kfft.h"

KFFT_API void
kfft_shift_cpx(kfft_cpx* buf, const uint32_t size, const bool is_inverse);
KFFT_API void
kfft_shift_scalar(kfft_scalar* buf, const uint32_t size, const bool is_inverse);
