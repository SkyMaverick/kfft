#pragma once

#include "kfft_config.h"
#include "kfft_simd_check.h"

#define __VEXST(S) KFFT_PLAN_VEX((S)) // compact version

#if defined(KFFT_SIMD_SSE_SUPPORT)
    #include "kfft_sse.h"

    #if defined(KFFT_HAVE_SSE3) // if define SSE3 support function
        #if defined(KFFT_HALF_SCALAR)
            #define VEX_CHECK_SSE(S) kfft_simd_check(__VEXST((S)), (HW_SSE2 | HW_SSE3))
        #else /* KFFT_HALF_SCALAR */
            #define VEX_CHECK_SSE(S) kfft_simd_check(__VEXST((S)), (HW_SSE | HW_SSE3))
        #endif /* KFFT_HALF_SCALAR */
    #else      /* KFFT_HAVE_SSE3 */
        #if defined(KFFT_HALF_SCALAR)
            #define VEX_CHECK_SSE(S) kfft_simd_check(__VEXST((S)), (HW_SSE2))
        #else /* KFFT_HALF_SCALAR */
            #define VEX_CHECK_SSE(S) kfft_simd_check(__VEXST((S)), (HW_SSE))
        #endif /*KFFT_HALF_SCALAR*/
    #endif

    #define __VFN_SEE(S, F, ...) (VEX_CHECK_SSE(S)) ? FUNC_SSE(F)(__VA_ARGS__)
#endif

/*!
  Macro for compile- and run-time selection
  vectorized optimization functions

    \param[in] S - plan object
    \param[in] F - function name
    \param[in] ... - __VA_ARGS__ function arguments

    \return function name (arg 2) return value
 */
#if defined(KFFT_USE_SIMD)
    #if defined(__VFN_SEE)
        #define VEXFUNC(S, F, ...) __VFN_SEE(S, F, __VA_ARGS__) : F(__VA_ARGS__)
    #else
        #define VEXFUNC(S, F, ...) F(__VA_ARGS__)
    #endif
#else                                         /* KFFT_USE_SIMD */
    #define VEXFUNC(S, F, ...) F(__VA_ARGS__) // default function call
#endif

//#undef __VEXST
//#undef __VFN_SEE
