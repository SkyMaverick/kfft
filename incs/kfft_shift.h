#pragma once
#include "kfft.h"

KFFT_API void
kfft_shift_cpx(kfft_cpx* buf, const uint32_t size, const bool is_inverse, kfft_pool_t* mmgr);
KFFT_API void
kfft_shift_scalar(kfft_scalar* buf, const uint32_t size, const bool is_inverse, kfft_pool_t* mmgr);

#if defined(KFFT_DYNAPI_ENABLE)
typedef void (*kfft_callback_shift_cpx)(kfft_cpx* buf, const uint32_t size, const bool is_inverse,
                                        kfft_pool_t* mmgr);
typedef void (*kfft_callback_shift_scalar)(kfft_scalar* buf, const uint32_t size,
                                           const bool is_inverse, kfft_pool_t* mmgr);
#endif /* KFFT_DYNAPI_ENABLE */
