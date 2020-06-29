#include "kfft_simd.h"

#if defined(KFFT_SIMD_SSE_SUPPORT)
    #if defined(KFFT_HALF_SCALAR)
        #include "sse/half/kfft_bfly.c"
        #include "sse/half/kfft_generic.c"
    #else
        #include "sse/norm/kfft_bfly.c"
        #include "sse/norm/kfft_generic.c"
    #endif
#endif
