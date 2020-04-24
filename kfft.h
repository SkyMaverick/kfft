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

#if defined(KFFT_HALF_SCALAR)
typedef float kfft_scalar;
#else
typedef double kfft_scalar;
#endif

typedef struct {
    kfft_scalar r;
    kfft_scalar i;
} kfft_cpx;

enum {
    KFFT_INFO_TRACE = 1 << 0,
    KFFT_INFO_USE_SIMD = 1 << 1,
    KFFT_INFO_USE_ALLOCA = 1 << 2,
    KFFT_INFO_USE_SYSMATH = 1 << 3,
    KFFT_INFO_USE_OPENMP = 1 << 4,
    KFFT_INFO_RADER_ALGO = 1 << 5,
    KFFT_INFO_MEMLESS_MODE = 1 << 6,
    KFFT_INFO_HALF_SCALAR = 1 << 7,
};

enum {
    KFFT_FLAG_NORMAL = 0,
    KFFT_FLAG_INVERSE = 1 << 0,
    KFFT_FLAG_RENEW = 1 << 1,
    KFFT_FLAG_GENERIC = 1 << 2,
    KFFT_FLAG_GENERIC_ONLY = 1 << 3,
};

typedef struct {
    uint16_t vmajor;
    uint16_t vminor;
    uint16_t vpatch;

    uint16_t flags;
} kfft_info_t;

enum {
    KFFT_RET_SUCCESS = 0x0000,
    KFFT_RET_ALLOC_FAIL = 0x0001,
    KFFT_RET_BUFFER_FAIL = 0x0002,
    KFFT_RET_FREE_NULL = 0x0003,
    KFFT_RET_IMPROPER_PLAN = 0x0004,
    KFFT_RET_BADARGUMENTS = 0x0005
};
typedef unsigned kfft_return_t;

#include "incs/kfft_macro.h"
#include "incs/kfft_trace.h"
#include "incs/kfft_system.h"

#if defined(KFFT_USE_SIMD)
    #include "kfft_simd.h"

    #define __ACCST(S) ((kfft_object_t*)(S))->accel.ext
// clang-format off
    #if defined(KFFT_SIMD_AVX2_SUPPORT)
        #define ACCELIT(S, F, ...)                                             \
            (__ACCST((S)) & HW_AVX2) ? FUNC_AVX2(F)(__VA_ARGS__) :             \
            (__ACCST((S)) & HW_AVX)  ? FUNC_AVX(F)(__VA_ARGS__)  :             \
            (__ACCST((S)) & (HW_SSE | HW_SSE2)) ? FUNC_SSE(F)(__VA_ARGS__) :   \
            F(__VA_ARGS__)
    #elif defined(KFFT_SIMD_AVX_SUPPORT)
        #define ACCELIT(S, F, ...)                                             \
            (__ACCST((S)) & HW_AVX)  ? FUNC_AVX(F)(__VA_ARGS__)  :             \
            (__ACCST((S)) & (HW_SSE | HW_SSE2)) ? FUNC_SSE(F)(__VA_ARGS__) :   \
            F(__VA_ARGS__)
    #elif defined(KFFT_SIMD_SSE_SUPPORT)
        #define ACCELIT(S, F, ...)                                             \
            (__ACCST((S)) & (HW_SSE | HW_SSE2)) ? FUNC_SSE(F)(__VA_ARGS__) :   \
            F(__VA_ARGS__)
    #else
        #define ACCELIT(S, F, ...)                                             \
            F(__VA_ARGS__)
    #endif
// clang-format on

#else
    #define ACCELIT(S, F, ...) F(__VA_ARGS__)
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
