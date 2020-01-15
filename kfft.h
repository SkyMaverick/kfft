#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

// clang-format off
#ifdef KFFT_USE_SIMD
    #include <xmmintrin.h>
    #define kfft_scalar __m128
    
    #define KFFT_MALLOC(nbytes) _mm_malloc(nbytes, 16)
    #define KFFT_FREE _mm_free
#else
    #define KFFT_MALLOC(X) calloc(1,(X))
    #define KFFT_FREE(X) free(X)
#endif
#define KFFT_ZEROMEM(M,X)  memset((M),0,(X))
// clang-format on

#define KFFT_FREE_NULL(X)                                                                          \
    do {                                                                                           \
        KFFT_FREE(X);                                                                              \
        X = NULL;                                                                                  \
    } while (0)

#ifndef kfft_scalar
    #define kfft_scalar double
#endif

#ifndef KFFT_RADER_LEVEL
    #define KFFT_RADER_LEVEL 3
#endif
#ifndef KFFR_RADER_LIMIT
    #define KFFT_RADER_LIMIT 50
#endif

#define KFFT_PLAN_ALLOCATOR(X) (*((uintptr_t*)(X)))

typedef struct {
    kfft_scalar r;
    kfft_scalar i;
} kfft_cpx;

typedef uintptr_t kfft_t;

kfft_t
kfft_config(const uint32_t nfft, const bool inverse_fft, const uintptr_t A, size_t* lenmem);
void
kfft(kfft_t cfg, const kfft_scalar* timedata, kfft_cpx* freqdata);
void
kffti(kfft_t cfg, const kfft_cpx* freqdata, kfft_scalar* timedata);
uint32_t
kfft_next_fast_size(uint32_t n);

void
kfft_free(kfft_t* cfg);

static inline uint32_t
kfft_get_size(const uint32_t n) {
    size_t memneeded = 0;
    kfft_config(n, 0, 0, &memneeded);
    return memneeded;
}

static inline bool
kfft_isnull(kfft_t in) {
    return (in == 0) ? true : false;
}

#ifdef __cplusplus
}
#endif
