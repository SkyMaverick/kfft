#include "kfft_simd.h"

#if defined(KFFT_SIMD_SSE_SUPPORT)
    #if defined(KFFT_HALF_SCALAR)
        #include "sse/half/kfft_sse.c"
    #else
        #include "sse/norm/kfft_sse.c"
    #endif
#endif
