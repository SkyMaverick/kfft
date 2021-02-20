#pragma once

#include "kfft_config.h"

#if defined(KFFT_SIMD_SSE_SUPPORT)
    #include "kfft_sse.h"
#endif

#include "kfft_simd_check.h"
