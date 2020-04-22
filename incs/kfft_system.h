#pragma once

// clang-format off
#if defined (KFFT_USE_SIMD)
    #include <immintrin.h>
//    #define kfft_scalar __m128
    
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

#ifdef KFFT_USE_ALLOCA
    // define this to allow use of alloca instead of malloc for temporary buffers
    // Temporary buffers are used in two case:
    // 1. FFT sizes that have "bad" factors. i.e. not 2,3 and 5
    // 2. "in-place" FFTs.  Notice the quotes, since kissfft does not really do an in-place
    // transform.
    #include <alloca.h>
    #if defined(KFFT_TRACE)
        #define KFFT_TMP_ALLOC(nbytes) alloca(nbytes)
        #define KFFT_TMP_FREE(ptr)
        #define KFFT_TMP_ALLOC(nbytes) alloca(nbytes)
        #define KFFT_TMP_FREE(ptr)
    #endif

    #define KFFT_TMP_ZEROMEM(M, X)                                                                 \
        do {                                                                                       \
            for (size_t i = 0; i < (X); i++)                                                       \
                ((char*)((M)))[i] = 0;                                                             \
        } while (0)
    #define KFFT_ALLOCA_CLEAR(M, X) KFFT_TMP_ZEROMEM((M), (X))
#else
    #if defined(KFFT_TRACE)
static inline void*
__trace_malloc(size_t nmem) {
    void* ret = KFFT_MALLOC(nmem);
    kfft_trace("[SYS] %s - %zu: %p\n", "Allocate temporary buffer", nmem, ret);
    return ret;
}
        #define KFFT_TMP_ALLOC(nbytes) __trace_malloc(nbytes)
        #define KFFT_TMP_FREE(ptr)                                                                 \
            do {                                                                                   \
                kfft_trace("[SYS] %s: %p\n", "Free temporary buffer", (void*)ptr);                 \
                KFFT_FREE(ptr);                                                                    \
            } while (0)
    #else
        #define KFFT_TMP_ALLOC(nbytes) KFFT_MALLOC(nbytes)
        #define KFFT_TMP_FREE(ptr) KFFT_FREE(ptr)
    #endif
    #define KFFT_TMP_ZEROMEM(M, X) KFFT_ZEROMEM(M, X)
    #define KFFT_ALLOCA_CLEAR(M, X)
#endif /* KFFT_USE_ALLOCA */

#define KFFT_TMP_FREE_NULL(X)                                                                      \
    do {                                                                                           \
        KFFT_TMP_FREE(X);                                                                          \
        X = NULL;                                                                                  \
    } while (0)
