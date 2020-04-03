#pragma once
#include "kfft.h"

KFFT_API void
kfft_shift(kfft_cpx* buf, const uint32_t size, const bool is_inverse);
