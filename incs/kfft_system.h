#pragma once

/*!
    \file
*/

#define KFFT_MALLOC_UNALIGN(X) calloc(1, (X))

#if defined(KFFT_USE_SIMD)
    #if defined(KFFT_ARCH_INTEL)
        #include <immintrin.h>
    #endif

    #if defined(KFFT_OS_WINDOWS)
        #include <malloc.h>
        #define KFFT_MALLOC(X, A) (((A) > 0) ? _aligned_malloc((X), (A)) : KFFT_MALLOC_UNALIGN((X)))
        #define KFFT_FREE(X, A) (((A) > 0) ? _aligned_free(X) : free(X))
    #else
        #define KFFT_MALLOC(X, A) (((A) > 0) ? aligned_alloc((A), (X)) : KFFT_MALLOC_UNALIGN((X)))
        #define KFFT_FREE(X, A) free(X)
    #endif
#else
    #define KFFT_MALLOC(X, A) KFFT_MALLOC_UNALIGN(X)
    #define KFFT_FREE(X, A) free(X)
#endif

#define KFFT_ZEROMEM(M, X) memset((M), 0, (X))

#define KFFT_FREE_NULL(X, A)                                                                       \
    do {                                                                                           \
        KFFT_FREE((X), (A));                                                                       \
        X = NULL;                                                                                  \
    } while (0)

#if (defined(KFFT_USE_ALLOCA)) && (!defined(KFFT_USE_SIMD))
    #include <alloca.h>
    #define KFFT_TMP_ALLOC(X, A) alloca((X))
    #define KFFT_TMP_FREE(X, A)

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
static inline void*
__trace_malloc(size_t nmem) {
    void* ret = KFFT_MALLOC_UNALIGN(nmem);
    kfft_trace("[SYS] %s - %zu: %p\n", "Allocate temporary buffer", nmem, ret);
    return ret;
}
    #if defined(KFFT_USE_SIMD)
static inline void*
__trace_malloc_aligned(size_t nmem, uint8_t align) {
    void* ret = KFFT_MALLOC(nmem, align);
    kfft_trace("[SYS] %s (%s %d) - %zu : %p\n", "Allocate temporary buffer", "aligned", align, nmem,
               ret);
    return ret;
}
        #if defined(KFFT_TRACE)
            #define KFFT_TMP_ALLOC(X, A)                                                           \
                (((A) > 0) ? __trace_malloc_aligned((X), (A)) : __trace_malloc((X)))
            #define KFFT_TMP_FREE(X, A)                                                            \
                do {                                                                               \
                    kfft_trace("[SYS] %s: %p\n", "Free temporary buffer", (void*)(X));             \
                    KFFT_FREE((X), (A));                                                           \
                } while (0)
        #else
            #define KFFT_TMP_ALLOC(X, A) KFFT_MALLOC((X), (A))
            #define KFFT_TMP_FREE(X, A) KFFT_FREE((X), (A))
        #endif /* KFFT_TRACE */

    #else /* KFFT_USE_SIMD */

        #if defined(KFFT_TRACE)
            #define KFFT_TMP_ALLOC(X, A) __trace_malloc((X))
            #define KFFT_TMP_FREE(X, A)                                                            \
                do {                                                                               \
                    kfft_trace("[SYS] %s: %p\n", "Free temporary buffer", (void*)(X));             \
                    KFFT_FREE((X), (A));                                                           \
                } while (0)
        #else
            #define KFFT_TMP_ALLOC(X, A) KFFT_MALLOC((X), (A))
            #define KFFT_TMP_FREE(X, A) KFFT_FREE((X), (A))
        #endif /* KFFT_TRACE */
    #endif     /* KFFT_USE_SIMD */

    #define KFFT_TMP_ZEROMEM(M, X) KFFT_ZEROMEM(M, X)
    #define KFFT_ALLOCA_CLEAR(M, X)

#endif /* KFFT_USE_ALLOCA */

#define KFFT_TMP_FREE_NULL(X, A)                                                                   \
    do {                                                                                           \
        KFFT_TMP_FREE((X), (A));                                                                   \
        X = NULL;                                                                                  \
    } while (0)

// OpenMP macroses

#if defined(KFFT_USE_OPENMP)
    #include <omp.h>
    #if (defined(_OPENMP) && (_OPENMP >= OMP_MINVER))
        #define KFFT_OMP(X) KFFT_PRAGMA(X)
    #else
        #define KFFT_OMP(X)
        #define KFFT_OMP_ISBLOCKED 1
    #endif
#else
    #define KFFT_OMP(X)
#endif
