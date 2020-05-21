#include "kfft.h"

KFFT_API void
kfft_info(kfft_info_t* info) {

    info->vmajor = KFFT_MAJOR_VERSION;
    info->vminor = KFFT_MINOR_VERSION;
    info->vpatch = KFFT_PATCH_VERSION;

#if defined(KFFT_TRACE)
    info->flags |= KFFT_INFO_TRACE;
#endif
#if defined(KFFT_USE_SIMD)
    info->flags |= KFFT_INFO_USE_SIMD;
#endif
#if defined(KFFT_USE_ALLOCA)
    info->flags |= KFFT_INFO_USE_ALLOCA;
#endif
#if defined(KFFT_USE_SYSMATH)
    info->flags |= KFFT_INFO_USE_SYSMATH;
#endif
#if defined(KFFT_USE_OPENMP)
    info->flags |= KFFT_INFO_USE_OPENMP;
#endif
#if defined(KFFT_RADER_ALGO)
    info->flags |= KFFT_INFO_RADER_ALGO;
#endif
#if defined(KFFT_MEMLESS_MODE)
    info->flags |= KFFT_INFO_MEMLESS_MODE;
#endif
#if defined(KFFT_HALF_SCALAR)
    info->flags |= KFFT_INFO_HALF_SCALAR;
#endif
}

KFFT_API uint32_t
kfft_next_fast_size(uint32_t n) {
    while (1) {
        uint32_t m = n;
        while ((m % 2) == 0)
            m /= 2;
        while ((m % 3) == 0)
            m /= 3;
        while ((m % 5) == 0)
            m /= 5;
        if (m <= 1)
            break; /* n is completely factorable by twos, threes, and fives */
        n++;
    }
    return n;
}

KFFT_API void*
kfft_malloc(uint32_t sz) {
#if defined(KFFT_USE_SIMD)
    return KFFT_MALLOC(sz, kfft_simd_align(kfft_simd_analize()));
#else
    return KFFT_MALLOC(sz, 0);
#endif
}
KFFT_API void
kfft_free_null(void** mem) {
#if defined(KFFT_USE_SIMD)
    KFFT_FREE(*mem, kfft_simd_align(kfft_simd_analize()));
#else
    KFFT_FREE(*mem, 0);
#endif
    *mem = NULL;
}

KFFT_API kfft_return_t
kfft_cleanup(void* mem) {
    kfft_object_t* M = (kfft_object_t*)(mem);
    if (M->mmgr == NULL)
        return KFFT_RET_FREE_NULL;

    kfft_pool_free(M->mmgr);
    return KFFT_RET_SUCCESS;
}

static const char* errors[] = {"Success operation result",
                               "Internal allocation memory fail",
                               "Temporary buffer allocation fail",
                               "Invalide free() operation for NULL buffer",
                               "Use improper real function for this plan",
                               "Use missing argument for configure plan"};

KFFT_API const char*
kfft_strerr(const kfft_return_t code) {
    return errors[code];
}
