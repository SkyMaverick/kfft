#pragma once

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#include "kfft_config.h"

#ifdef __cplusplus
extern "C" {
#endif

#include "incs/kfft_macro.h"
#include "incs/kfft_trace.h"
#include "incs/kfft_system.h"

#include "incs/kfft_types.h"
#if defined(KFFT_USE_SIMD)
    #include "kfft_simd.h"

    #define __VEXST(S) ((kfft_object_t*)(S))->vex
// clang-format off
    #if defined(KFFT_SIMD_AVX2_SUPPORT)
        #define VEXFUNC(S, F, ...)                                                              \
            ( kfft_simd_check(__VEXST((S)),HW_AVX2) )            ? FUNC_AVX2(F)(__VA_ARGS__) :  \
            ( kfft_simd_check(__VEXST((S)),HW_AVX) )             ? FUNC_AVX (F)(__VA_ARGS__) :  \
            ( kfft_simd_check(__VEXST((S)),(HW_SSE | HW_SSE2)) ) ? FUNC_SSE (F)(__VA_ARGS__) :  \
            F(__VA_ARGS__)
    #elif defined(KFFT_SIMD_AVX_SUPPORT)
        #define VEXFUNC(S, F, ...)                                                              \
            ( kfft_simd_check(__VEXST((S)),HW_AVX) )             ? FUNC_AVX (F)(__VA_ARGS__) :  \
            ( kfft_simd_check(__VEXST((S)),(HW_SSE | HW_SSE2)) ) ? FUNC_SSE (F)(__VA_ARGS__) :  \
            F(__VA_ARGS__)
    #elif defined(KFFT_SIMD_SSE_SUPPORT)
        #define VEXFUNC(S, F, ...)                                                              \
            ( kfft_simd_check(__VEXST((S)),(HW_SSE | HW_SSE2)) ) ? FUNC_SSE (F)(__VA_ARGS__) :  \
            F(__VA_ARGS__)
    #else
        #define VEXFUNC(S, F, ...)                                                              \
            F(__VA_ARGS__)
    #endif
// clang-format on

#else
    #define VEXFUNC(S, F, ...) F(__VA_ARGS__)
#endif

#include "incs/kfft_math.h"
#include "incs/kfft_alloc.h"
#include "incs/kfft_ext.h"
#include "incs/kfft_shift.h"

#include "incs/kfft_cpx.h"
#include "incs/kfft_scalar.h"

#if defined(KFFT_2D_ENABLE)
    #include "2d/kfft_cpx2.h"
    #include "2d/kfft_scalar2.h"
#endif

#define KFFT_PLAN_ALLOCATOR(X) (*((kfft_pool_t**)(X)))

#define kfft_free(X) kfft_cleanup((uintptr_t)(X))

/* Protecting nested plans from destructive operations */
#define KFFT_CHECK_FLAGS(X) ((X) & (~KFFT_FLAG_RENEW))

#ifdef __cplusplus
}
#endif
