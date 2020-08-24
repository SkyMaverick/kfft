#pragma once
#include "kfft.h"

/*!
    \file
    Shift zero-frequency component to center of spectrum function.
 */
/*!
    FFT shift function for 1D complex sequence
    \param[in] buf - complex buffer for process operation
    \param[in] size - buffer size (count)
    \param[in] is_inverse - forward/inverse operation flag
    \param[in] mmgr - plan memory manager (use ::KFFT_PLAN_MMGR)
    \return None
 */
KFFT_API void
kfft_shift_cpx(kfft_cpx* buf, const uint32_t size, const bool is_inverse, kfft_pool_t* mmgr);
/*!
    FFT shift function for 1D complex sequence
    \param[in] buf - scalar buffer for process operation
    \param[in] size - buffer size (count)
    \param[in] is_inverse - forward/inverse operation flag
    \param[in] mmgr - plan memory manager (use ::KFFT_PLAN_MMGR)
    \return None
 */
KFFT_API void
kfft_shift_scalar(kfft_scalar* buf, const uint32_t size, const bool is_inverse, kfft_pool_t* mmgr);

#if defined(KFFT_DYNAPI_ENABLE)
// clang-format off
typedef void
(*kfft_callback_shift_cpx)(kfft_cpx* buf, const uint32_t size, const bool is_inverse, kfft_pool_t* mmgr);
typedef void
(*kfft_callback_shift_scalar)(kfft_scalar* buf, const uint32_t size, const bool is_inverse, kfft_pool_t* mmgr);
// clang-format on
#endif /* KFFT_DYNAPI_ENABLE */
