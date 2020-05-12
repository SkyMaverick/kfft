#pragma once

#if defined(KFFT_USE_SIMD)
    #include <immintrin.h>

    #if defined(KFFT_OS_WINDOWS)
        #include <malloc.h>
        #define KFFT_MALLOC(X, A) _aligned_malloc((X), (A))
        #define KFFT_FREE(X) _aligned_free(X)
    #else
        #define KFFT_MALLOC(X, A) aligned_alloc((A), (X))
        #define KFFT_FREE(X) free(X)
    #endif
#else
    #define KFFT_MALLOC(X) calloc(1, (X))
    #define KFFT_FREE(X) free(X)
#endif
#define KFFT_ZEROMEM(M, X) memset((M), 0, (X))

#define KFFT_FREE_NULL(X)                                                                          \
    do {                                                                                           \
        KFFT_FREE(X);                                                                              \
        X = NULL;                                                                                  \
    } while (0)

#if (defined(KFFT_USE_ALLOCA)) && (!defined(KFFT_USE_SIMD))
    #include <alloca.h>
    #define KFFT_TMP_ALLOC(X) alloca((X))
    #define KFFT_TMP_FREE(X)

    #define KFFT_TMP_ZEROMEM(M, X)                                                                 \
        do {                                                                                       \
            for (size_t i = 0; i < (X); i++)                                                       \
                ((char*)((M)))[i] = 0;                                                             \
        } while (0)
    #define KFFT_ALLOCA_CLEAR(M, X) KFFT_TMP_ZEROMEM((M), (X))
#else /*not ALLOCA */
    #if defined(KFFT_USE_ALLOCA)
        #pragma message WARN("SIMD functions need aligned memory. Disable alloca functionality")
    #endif /* KFFT_USE_ALLOCA */
    #if defined(KFFT_USE_SIMD)
        #if defined(KFFT_TRACE)
static inline void*
__trace_malloc_aligned(size_t nmem, uint8_t align) {
    void* ret = KFFT_MALLOC(nmem, align);
    kfft_trace("[SYS] %s - %zu: %p\n", "Allocate temporary buffer", nmem, ret);
    return ret;
}
            #define KFFT_TMP_ALLOC(X, A) __trace_malloc_aligned((X), (A))
            #define KFFT_TMP_FREE(X)                                                               \
                do {                                                                               \
                    kfft_trace("[SYS] %s: %p\n", "Free temporary buffer", (void*)(X));             \
                    KFFT_FREE((X));                                                                \
                } while (0)
        #else
            #define KFFT_TMP_ALLOC(X, A) KFFT_MALLOC((X), (A))
            #define KFFT_TMP_FREE(X) KFFT_FREE((X))
        #endif /* KFFT_TRACE */
    #else      /* KFFT_USE_SIMD */
        #if defined(KFFT_TRACE)
static inline void*
__trace_malloc(size_t nmem) {
    void* ret = KFFT_MALLOC(nmem);
    kfft_trace("[SYS] %s - %zu: %p\n", "Allocate temporary buffer", nmem, ret);
    return ret;
}
            #define KFFT_TMP_ALLOC(X) __trace_malloc((X))
            #define KFFT_TMP_FREE(X)                                                               \
                do {                                                                               \
                    kfft_trace("[SYS] %s: %p\n", "Free temporary buffer", (void*)(X));             \
                    KFFT_FREE((X));                                                                \
                } while (0)
        #else
            #define KFFT_TMP_ALLOC(X) KFFT_MALLOC((X))
            #define KFFT_TMP_FREE(X) KFFT_FREE((X))
        #endif /* KFFT_TRACE */
    #endif     /* KFFT_USE_SIMD */

    #define KFFT_TMP_ZEROMEM(M, X) KFFT_ZEROMEM(M, X)
    #define KFFT_ALLOCA_CLEAR(M, X)

#endif /* KFFT_USE_ALLOCA */

#define KFFT_TMP_FREE_NULL(X)                                                                      \
    do {                                                                                           \
        KFFT_TMP_FREE(X);                                                                          \
        X = NULL;                                                                                  \
    } while (0)
