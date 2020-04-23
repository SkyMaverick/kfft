#pragma once

#include "kfft_config.h"

#include <inttypes.h>
#include <stdbool.h>

typedef struct {
    uint8_t arch; // Architecture ID
    uint32_t ext; // HW extensions extensionse (with operation system correct)
} kfft_simd_t;

#if defined(KFFT_SIMD_SSE_SUPPORT)
    #define FUNC_SSE(X) X##_sse
    #include "sse/kfft_math_sse.h"
#endif
#if defined(KFFT_SIMD_AVX_SUPPORT)
    #define FUNC_AVX(X) X##_avx
#endif
#if defined(KFFT_SIMD_AVX2_SUPPORT)
    #define FUNC_AVX2(X) X##_avx2
#endif

#include "kfft_simd_check.h"
