#pragma once

#include "kfft_config.h"

#if defined(KFFT_SIMD_SSE_SUPPORT)
    #define FUNC_SSE(X) X##_sse
    #if defined(KFFT_HALF_SCALAR)
        #include "sse/half/kfft_math_sse.h"
        #include "sse/half/kfft_bfly.h"
        #include "sse/half/kfft_generic.h"
    #else
        #include "sse/norm/kfft_math_sse.h"
        #include "sse/norm/kfft_bfly.h"
        #include "sse/norm/kfft_generic.h"
    #endif
#endif
#if defined(KFFT_SIMD_AVX_SUPPORT)
    #define FUNC_AVX(X) X##_avx

    #include "avx/kfft_math_avx.h"
    #include "avx/kfft_bfly.h"
    #include "avx/kfft_generic.h"
#endif

#include "kfft_simd_check.h"
