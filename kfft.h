/*
 *  Copyright (c) 2019 - 2020, Alexander Smirnov. All rights reserved.
 *  This file is part of KFFT - https://github.com/SkyMaverick/kfft
 *
 *  SPDX-License-Identifier: Zlib
 *  See COPYING file for more information.
 */

/*!
 * @mainpage
 *
 * \author Alexander Smirnov (https://github.com/SkyMaverick)
 * \version 0.7.1
 *
 * LibKFFT is a free/open-source (ZLib licensed) C library
 * providing fast fourier transform and other related operations
 * for comlex and scalar double or float sequences.
 * LibKFFT is deeply redesigned fork KissFFT project by Mark Borgerding
 * (https://github.com/mborgerding/kissfft).
 *
 * See the LibKFFT page on GitHub for downloads and other information,
 * or the source code (https://github.com/SkyMaverick/kfft).
 *
 * The features of LibKFFT include:
 * - complex and scalar fourier transformation operations (forward and inverse)
 *      - 64-bit float type (double) as default
 *      - 32-bit float type (float) if set option KFFT_HALF_SCALAR
 * - SIMD acceleration transform operations (now support SSE, SSE2, SSE3 (optional)) (build option)
 * - Recursive Rader algorithm (https://en.wikipedia.org/wiki/Rader's_FFT_algorithm) for
 * prime-number lenght sequenses (build option)
 * - Internal math functions (aka sqrt, sincos etc.) without system libm functions (build
 * option)
 * - Parallelization with OpenMP framework (build option)
 * - 2D, sparse, convolution(+2D) operations as build-time extensions
 * - Detail trace log in stdout for diagnostic (build option)
 * - Tiny API support for simplify code in trivial usage
 *
 * Test build on OS with CC:
 * - Linux (x86-64): GCC 9.3, Clang-10
 * - Windows (x86-64): MSVC 2015, 2017, Clang-10
 * - RPi3 (armhf): GCC 8.3
 */

/** @file */

#pragma once

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#include "kfft_config.h"

#ifdef __cplusplus
extern "C" {
#endif

#include "incs/kfft_macro.h"  // OS and CC predefined macroses
#include "incs/kfft_trace.h"  // internal traces function
#include "incs/kfft_system.h" // system abstractions layer

#include "incs/kfft_types.h" // core basis types and enums
#if defined(KFFT_USE_SIMD)
    #include "kfft_simd.h" // vector optimized operation proveder
#endif

#include "incs/kfft_math.h" // internal primitive math operations
#include "incs/kfft_pool.h" // plan simple allocator
#include "incs/kfft_algo.h" // standart macro algorithms
#include "incs/kfft_ext.h"  // service functions (API)

#include "incs/kfft_cpx.h"    // complex FFT analysis functions (API)
#include "incs/kfft_scalar.h" // scalar FFT analysis functions (API)

#include "incs/kfft_shift.h" // complex/scalar shift functions (API)

/* Functions provided by the extension [2d] */
#if defined(KFFT_2D_ENABLE)
    #include "2d/kfft_cpx2.h"    // complex 2D FFT analysis functions (API)
    #include "2d/kfft_scalar2.h" // scalar 2D FFT analysis functions (API)
#endif

/* Functions provided by the extension [sparse] */
#if defined(KFFT_SPARSE_ENABLE)
    #include "sparse/kfft_cpx_sparse.h"    // complex sparse FFT analysis functions (API)
    #include "sparse/kfft_scalar_sparse.h" // scalar sparse FFT analysis functions (API)
#endif

/* Functions provided by the extension [conv] */
#if defined(KFFT_CONV_ENABLE)
    #include "conv/kfft_cpx_conv.h"    // complex convolution functions (API)
    #include "conv/kfft_scalar_conv.h" // scalar convolution functions (API)
#endif

/* Functions provided by the extension [conv2d] */
#if defined(KFFT_CONV2D_ENABLE)
    #include "conv2d/kfft_cpx_conv.h"    // complex 2D convolution functions (API)
    #include "conv2d/kfft_scalar_conv.h" // scalar 2D convolution functions (API)
#endif

/*!
  Macro to get the memory manager of a plan object
  \param[in] X - plan object
  \return plan memory manager pointer (kfft_pool_t*)
 */
#define KFFT_PLAN_MMGR(X) (*((kfft_pool_t**)(X)))

/*!
  Macro to get the memory manager of a plan object if it's enabled
  or NULL if plan is NULL. Required for the assignment operation.
  \param[in] X - plan object
  \return plan memory manager pointer (kfft_pool_t*) or NULL
 */
#define KFFT_PLAN_MMGR_NULL(X) ((X != NULL) ? KFFT_PLAN_MMGR((X)) : NULL)

/*!
  Macro to get the align information of a plan object
  \param[in] X - plan object
  \return uint8_t number memory align for plan (0, 16, 32)
 */
#define KFFT_PLAN_ALIGN(X) ((X != NULL) ? KFFT_PLAN_MMGR((X))->align : 0)

/*!
  Macro to get the acceleration info of a plan object

  \param[in] X - plan object
  \warning argument (X) mustn't NULL

  \return kfft_simd_t vector extension information
 */
#define KFFT_PLAN_VEX(X) KFFT_PLAN_MMGR((X))->vex

/*!
  Macro to protecting nested plans from destructive operations
  \param[in] X - plan object
  \return None
 */
#define KFFT_CHECK_FLAGS(X) ((X) & (~KFFT_FLAG_RENEW))

/*!
  Macro for compile- and run-time selection
  vectorized optimization functions

    \param[in] S - plan object
    \param[in] F - function name
    \param[in] ... - __VA_ARGS__ function arguments

    \return function name (arg 2) return value
 */
#if defined(KFFT_USE_SIMD)
    #define __VEXST(S) KFFT_PLAN_VEX((S)) // compact version
    // clang-format off
    #if defined(KFFT_HAVE_SSE3) // if define SSE3 support function
        #if defined(KFFT_HALF_SCALAR)
            #define VEX_CHECK_SSE(S)                                                                \
                kfft_simd_check(__VEXST((S)),(HW_SSE2 | HW_SSE3))
        #else /* KFFT_HALF_SCALAR */
            #define VEX_CHECK_SSE(S)                                                                \
                kfft_simd_check(__VEXST((S)),(HW_SSE | HW_SSE3))
        #endif /* KFFT_HALF_SCALAR */
    #else /* KFFT_HAVE_SSE3 */
        #if defined(KFFT_HALF_SCALAR)
            #define VEX_CHECK_SSE(S)                                                                \
                kfft_simd_check(__VEXST((S)),(HW_SSE2))
        #else /* KFFT_HALF_SCALAR */
            #define VEX_CHECK_SSE(S)                                                                \
                kfft_simd_check(__VEXST((S)),(HW_SSE))
        #endif /*KFFT_HALF_SCALAR*/
    #endif

    #if defined(KFFT_SIMD_SSE_SUPPORT)
        #define VEXFUNC(S, F, ...)                                                              \
            ( VEX_CHECK_SSE(S) ) ? FUNC_SSE (F)(__VA_ARGS__) :                                  \
            F(__VA_ARGS__)
    #else
        #define VEXFUNC(S, F, ...)                                                              \
            F(__VA_ARGS__)
    #endif
// clang-format on

#else                                         /* KFFT_USE_SIMD */
    #define VEXFUNC(S, F, ...) F(__VA_ARGS__) // default function call
#endif
#ifdef __cplusplus
}
#endif
