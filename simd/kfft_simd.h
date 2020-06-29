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

#include "kfft_simd_check.h"
