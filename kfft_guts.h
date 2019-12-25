/*
Copyright (c) 2003-2010, Mark Borgerding

All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice, this list of conditions
and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice, this list of
conditions and the following disclaimer in the documentation and/or other materials provided with
the distribution.
    * Neither the author nor the names of any contributors may be used to endorse or promote
products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/* kfft.h
   defines kfft_scalar as either short or a float type
   and defines
   typedef struct { kfft_scalar r; kfft_scalar i; }kfft_cpx; */

#pragma once

#include <inttypes.h>
#include <limits.h>

#include "kfft.h"

#define MAX_FACTORS 32
#define MAX_ROOTS 32
#define MAX_BFLY_LEVEL 5

#ifndef USE_SYSMATH
    #define KFFT_CONST_PI 3.141592653589793238462643383279502884197169399375105820974944
#else
    #include <math.h>
    #define KFFT_CONST_PI M_PI
#endif

/* e.g. an fft of length 128 has 4 factors
 as far as kissfft is concerned
 4*4*4*2
 */

typedef struct kfft_kstate {
    uint32_t nfft;
    bool inverse;
    uint8_t level;
    uint32_t factors[2 * MAX_FACTORS];
#ifndef KFFT_MEMLESS_MODE
    #ifdef KFFT_RADER_ALGO
    uint32_t roots[3 * MAX_ROOTS];
    #endif
    kfft_cpx twiddles[1];
#endif /* memless */
} kfft_kplan_t;

typedef struct kfft_state {
    kfft_kplan_t* substate;
    kfft_cpx* tmpbuf;
#ifndef KFFT_MEMLESS_MODE
    kfft_cpx* super_twiddles;
#endif /* memless */
#ifdef USE_SIMD
    void* pad;
#endif
} kfft_plan_t;

/*
  Explanation of macros dealing with complex math:

   C_MUL(m,a,b)         : m = a*b
   C_SUB( res, a,b)     : res = a - b
   C_SUBFROM( res , a)  : res -= a
   C_ADDTO( res , a)    : res += a
 * */

#define S_MUL(a, b) ((a) * (b))
#define S_DIV(a, b) ((a) / (b))

#define C_MUL(m, a, b)                                                                             \
    do {                                                                                           \
        (m).r = (a).r * (b).r - (a).i * (b).i;                                                     \
        (m).i = (a).r * (b).i + (a).i * (b).r;                                                     \
    } while (0)
#define C_MULBYSCALAR(c, s)                                                                        \
    do {                                                                                           \
        (c).r *= (s);                                                                              \
        (c).i *= (s);                                                                              \
    } while (0)
#define C_DIVBYSCALAR(c, s)                                                                        \
    do {                                                                                           \
        (c).r /= (s);                                                                              \
        (c).i /= (s);                                                                              \
    } while (0)

#define C_ADD(res, a, b)                                                                           \
    do {                                                                                           \
        (res).r = (a).r + (b).r;                                                                   \
        (res).i = (a).i + (b).i;                                                                   \
    } while (0)
#define C_SUB(res, a, b)                                                                           \
    do {                                                                                           \
        (res).r = (a).r - (b).r;                                                                   \
        (res).i = (a).i - (b).i;                                                                   \
    } while (0)
#define C_ADDTO(res, a)                                                                            \
    do {                                                                                           \
        (res).r += (a).r;                                                                          \
        (res).i += (a).i;                                                                          \
    } while (0)

#define C_SUBFROM(res, a)                                                                          \
    do {                                                                                           \
        (res).r -= (a).r;                                                                          \
        (res).i -= (a).i;                                                                          \
    } while (0)

#if defined(USE_SIMD)
    #define KFFT_COS(phase) _mm_set1_ps(cos(phase))
    #define KFFT_SIN(phase) _mm_set1_ps(sin(phase))
    #define HALF_OF(x) ((x)*_mm_set1_ps(.5))
#else
    #define KFFT_COS(phase) (kfft_scalar) cos(phase)
    #define KFFT_SIN(phase) (kfft_scalar) sin(phase)
    #define HALF_OF(x) ((x)*.5)
#endif

#define kf_cexp(x, phase)                                                                          \
    do {                                                                                           \
        (x)->r = KFFT_COS(phase);                                                                  \
        (x)->i = KFFT_SIN(phase);                                                                  \
    } while (0)

/* a debugging function */
#define pcpx(c) fprintf(stderr, "%g + %gi\n", (double)((c)->r), (double)((c)->i))

#ifndef KFFT_MEMLESS_MODE
    #define TWIDDLE(i, P) P->twiddles[i]
#else
static inline kfft_cpx
get_kernel_twiddle(uint32_t i, const kfft_kplan_t* P) {
    kfft_cpx ret;

    kfft_scalar phase = -2 * KFFT_CONST_PI * i / P->nfft;
    if (P->inverse)
        phase *= -1;

    kf_cexp(&ret, phase);
    return ret;
}

    #define TWIDDLE(i, P) get_kernel_twiddle(i, P)
#endif /* memless */

#ifdef KFFT_USE_ALLOCA
    // define this to allow use of alloca instead of malloc for temporary buffers
    // Temporary buffers are used in two case:
    // 1. FFT sizes that have "bad" factors. i.e. not 2,3 and 5
    // 2. "in-place" FFTs.  Notice the quotes, since kissfft does not really do an in-place
    // transform.
    #include <alloca.h>
    #define KFFT_TMP_ALLOC(nbytes) alloca(nbytes)
    #define KFFT_TMP_FREE(ptr)
#else
    #define KFFT_TMP_ALLOC(nbytes) KFFT_MALLOC(nbytes)
    #define KFFT_TMP_FREE(ptr) KFFT_FREE(ptr)
#endif

#if (defined TRACE)
    #define kfft_trace(fmt, ...) fprintf(stdout, fmt, __VA_ARGS__)
#else
    #define kfft_trace(fmt, ...) // noop
#endif
