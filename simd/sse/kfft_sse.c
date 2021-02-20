#include "kfft_sse.h"

#if defined(KFFT_HALF_SCALAR)
    #include "half/kfft_math.c"
    #include "half/kfft_bfly.c"
    #include "half/kfft_generic.c"
#else
    #include "norm/kfft_math.c"
    #include "norm/kfft_bfly.c"
    #include "norm/kfft_generic.c"
#endif
