#pragma once

#include "kfft_config.h"

#if defined(KFFT_SIMD_SSE_SUPPORT)
    #define FUNC_SSE(X) X##_sse

    #include "sse/kfft_math_sse.h"
    #include "sse/kfft_bfly.h"
    #include "sse/kfft_generic.h"
#endif
#if defined(KFFT_SIMD_AVX_SUPPORT)
    #define FUNC_AVX(X) X##_avx
#endif
#if defined(KFFT_SIMD_AVX2_SUPPORT)
    #define FUNC_AVX2(X) X##_avx2
#endif

#include "kfft_simd_check.h"
