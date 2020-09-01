#pragma once

/*!
    \file
    \brief Operation system abstraction layer (OSAL)

    Wrappers over the operating system functions required
    for the library to work.
*/

/***************************************************************************************************
 *                                   SYSTEM MEMORY OPERATIONS
 ***************************************************************************************************/

/*!
  Standart memory allocation macro.
  \note It is not recommended to use this macro directly.
  \param[in] X - memory size (bytes)
  \result unaligned memory area pointer
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
#else /* KFFT_USE_SIMD */
    /*!
       \brief Memory allocation macro.

       This macro allocate unaligned memory (as ::KFFT_MALLOC_UNALIGN) if don't use vector CPU
       extensions (VEX) or allocate aligned memory if use VEX.

       \param[in] X - memory size (bytes)
       \param[in] A - aligned bytes (maybe 0[unalign], 16, 32)
       \result aligned memory area pointer
     */
    #define KFFT_MALLOC(X, A) KFFT_MALLOC_UNALIGN(X)
    /*!
        Memory release macro

        Need param A (aligned) for OS Windows because _aligned_malloc don't use 0. Other
        OS ignore this paramer

       \param[in] X - memory area pointer
       \param[in] A - aligned bytes (maybe 0[unalign], 16, 32)
       \warning A must equivalent ::KFFT_MALLOC param A for this buffer
       \result None
     */
    #define KFFT_FREE(X, A) free(X)
#endif

/*!
    Clear memory area.
    \param[in] M - memory area pointer
    \param[in] S - size area pointer (bytes)
    \result None
 */
#define KFFT_ZEROMEM(M, S) memset((M), 0, (S))

/*!
    Memory release and NULL variable macro
    \param[in] X - memory area pointer
    \param[in] A - aligned bytes (maybe 0[unalign], 16, 32)
    \warning A must equivalent ::KFFT_MALLOC param A for this buffer
    \result None
 */
#define KFFT_FREE_NULL(X, A)                                                                       \
    do {                                                                                           \
        KFFT_FREE((X), (A));                                                                       \
        X = NULL;                                                                                  \
    } while (0)

/***************************************************************************************************
 *                                   TEMPORARY BUFFERS OPERATIONS
 ***************************************************************************************************/

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
            /*!
                Memory allocate for temporary buffers. Maybe use alloca() if set this option.
                @see KFFT_MALLOC
             */
            #define KFFT_TMP_ALLOC(X, A) KFFT_MALLOC((X), (A))
            /*!
                Memory release macro for temporary buffers.
                @see KFFT_FREE
             */
            #define KFFT_TMP_FREE(X, A) KFFT_FREE((X), (A))
        #endif /* KFFT_TRACE */
    #endif     /* KFFT_USE_SIMD */

    /*!
        Clear memory area for temporary buffers
        @see KFFT_ZEROMEM
     */
    #define KFFT_TMP_ZEROMEM(M, S) KFFT_ZEROMEM(M, S)
    /*!
        Clear memory area if it's area allocate with alloca() function.
        \param[in] M - memory area pointer
        \param[in] S - area size (bytes)
        \result None
     */
    #define KFFT_ALLOCA_CLEAR(M, S)

#endif /* KFFT_USE_ALLOCA */

/*!
    Memory release and NULL variable macro for temporary buffers
    @see KFFT_FREE_NULL
 */
#define KFFT_TMP_FREE_NULL(X, A)                                                                   \
    do {                                                                                           \
        KFFT_TMP_FREE((X), (A));                                                                   \
        X = NULL;                                                                                  \
    } while (0)

/***************************************************************************************************
 *                                   OPENMP OPERATIONS PRAGMA WRAPPERS
 ***************************************************************************************************/

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
    /*!
        Start OpenMP block if OpenMP functionality is available.
        \param[in] X - pragma text
     */
    #define KFFT_OMP(X)
#endif
