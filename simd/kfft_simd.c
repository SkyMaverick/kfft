#include "kfft_simd.h"

#if defined(KFFT_SIMD_SSE_SUPPORT)
    #include "sse/kfft_bfly.c"
    #include "sse/kfft_generic.c"
#endif
#if defined(KFFT_SIMD_AVX_SUPPORT)
    #include "avx/kfft_bfly.c"
    #include "avx/kfft_generic.c"
#endif
#if defined(KFFT_SIMD_AVX2_SUPPORT)
// TODO
//    #include "avx2/kfft_bfly.c"
//    #include "avx2/kfft_generic.c"
#endif
