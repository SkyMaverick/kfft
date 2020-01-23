#include "kfft_guts.h"

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
#if defined(KFFT_RADER_ALGO)
    info->flags |= KFFT_INFO_RADER_ALGO;
#endif
#if defined(KFFT_MEMLESS_MODE)
    info->flags |= KFFT_INFO_MEMLESS_MODE;
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
