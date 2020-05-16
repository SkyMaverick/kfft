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
#endif

#include "incs/kfft_math.h"
#include "incs/kfft_alloc.h"
#include "incs/kfft_algo.h"
#include "incs/kfft_ext.h"

#include "incs/kfft_cpx.h"
#include "incs/kfft_scalar.h"

#include "incs/kfft_shift.h"

#if defined(KFFT_2D_ENABLE)
    #include "2d/kfft_cpx2.h"
    #include "2d/kfft_scalar2.h"
#endif
#if defined(KFFT_SPARSE_ENABLE)
    #include "sparse/kfft_cpx_sparse.h"
    #include "sparse/kfft_scalar_sparse.h"
#endif

#define KFFT_PLAN_MMGR(X) (*((kfft_pool_t**)(X)))

#define KFFT_PLAN_ALIGN(X) KFFT_PLAN_MMGR((X))->align
#define KFFT_PLAN_VEX(X) KFFT_PLAN_MMGR((X))->vex

#define kfft_free(X) kfft_cleanup((uintptr_t)(X))

/* Protecting nested plans from destructive operations */
#define KFFT_CHECK_FLAGS(X) ((X) & (~KFFT_FLAG_RENEW))

#if defined(KFFT_USE_SIMD)
    #define __VEXST(S) KFFT_PLAN_VEX((S))
// clang-format off
    #define VEX_CHECK_AVX2(S)                                                                   \
        kfft_simd_check(__VEXST((S)),HW_AVX2)
    #define VEX_CHECK_AVX(S)                                                                    \
       kfft_simd_check(__VEXST((S)),HW_AVX)
    #if defined(KFFT_HAVE_SSE3)
        #define VEX_CHECK_SSE(S)                                                                \
            kfft_simd_check(__VEXST((S)),(HW_SSE2 | HW_SSE3))
    #else
        #define VEX_CHECK_SSE(S)                                                                \
           kfft_simd_check(__VEXST((S)),(HW_SSE2))
    #endif

    #if defined(KFFT_SIMD_AVX2_SUPPORT)
        #define VEXFUNC(S, F, ...)                                                              \
            ( VEX_CHECK_AVX2(S) ) ? FUNC_AVX2(F)(__VA_ARGS__) :                                 \
            ( VEX_CHECK_AVX(S)  ) ? FUNC_AVX (F)(__VA_ARGS__) :                                 \
            ( VEX_CHECK_SSE(S)  ) ? FUNC_SSE (F)(__VA_ARGS__) :                                 \
            F(__VA_ARGS__)
    #elif defined(KFFT_SIMD_AVX_SUPPORT)
        #define VEXFUNC(S, F, ...)                                                              \
            ( VEX_CHECK_AVX(S) ) ? FUNC_AVX (F)(__VA_ARGS__) :                                  \
            ( VEX_CHECK_SSE(S) ) ? FUNC_SSE (F)(__VA_ARGS__) :                                  \
            F(__VA_ARGS__)
    #elif defined(KFFT_SIMD_SSE_SUPPORT)
        #define VEXFUNC(S, F, ...)                                                              \
            ( VEX_CHECK_SSE(S) ) ? FUNC_SSE (F)(__VA_ARGS__) :                                  \
            F(__VA_ARGS__)
    #else
        #define VEXFUNC(S, F, ...)                                                              \
            F(__VA_ARGS__)
    #endif
// clang-format on

#else /* KFFT_USE_SIMD */
    #define VEXFUNC(S, F, ...) F(__VA_ARGS__)
#endif
#ifdef __cplusplus
}
#endif
