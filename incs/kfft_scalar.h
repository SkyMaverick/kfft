#pragma once

#include "kfft_cpx.h"
/*!
    \file
    \brief Scalar (timedata) 1D-buffer evaluations functionality.
 */

// clang-format off
    #define kfft_trace_scalar(fmt, ...)                                                                  \
        kfft_trace("[SCLR]"" " fmt, __VA_ARGS__)
// clang-format on

/// 1D scalar operations plan type
typedef struct kfft_state {
    kfft_object_t object; ///< standart object structure

    uint32_t nfft;  ///< sequence lenght (count)
    uint32_t flags; ///< operation flags ::kfft_eval_flags

    kfft_plan_cpx* basis;     ///< internal core complex plan
    kfft_cpx* super_twiddles; ///< extended twiddles for scalar evaluation
} kfft_plan_sclr;

/*!
    \brief Scalar plan create or configuration functions
    \param[in] nfft - lenght input sequense
    \param[in] flags - operation flags
    \param[in] A - plan ainternal allocator structure (if need use KFFT_PLAN_MMGR macro)
    \param[in] lenmem - vaiable for memory get pointer
    \return scalar plan structure or NULL

    @see kfft_config_cpx for more information of usage this API
 */
KFFT_API kfft_plan_sclr*
kfft_config_scalar(const uint32_t nfft, const uint32_t flags, kfft_pool_t* A, size_t* lenmem);
/*!
  \brief Process forward evaluation function.
  \param[in] plan - scalar plan ::kfft_plan_sclr pointer
  \param[in] fin - input ::kfft_scalar buffer (don't changed)
  \param[in] fout - output ::kfft_cpx buffer
  \result standart ::kfft_ret_flags return status

  \warning Function NOT control input and output buffer (such as the NULL, overflow, bad
  buffers-size etc.). Developer must control this manualy.
 */
KFFT_API kfft_return_t
kfft_eval_scalar(kfft_plan_sclr* plan, const kfft_scalar* fin, kfft_cpx* fout);
/*!
  \brief Process forward evaluation function.
  \param[in] plan - scalar plan ::kfft_plan_sclr pointer
  \param[in] fin - input ::kfft_scalar buffer (don't changed)
  \param[in] fout - output ::kfft_scalar buffer
  \result standart ::kfft_ret_flags return status

  \note Temporary use normalization with formula
    ![Absolute complex value](norm_abs.svg)

  \warning Function NOT control input and output buffer (such as the NULL, overflow, bad
  buffers-size etc.). Developer must control this manualy.
 */
KFFT_API kfft_return_t
kfft_eval_scalar_norm(kfft_plan_sclr* plan, const kfft_scalar* fin, kfft_scalar* fout);
/*!
  \brief Process inverse evaluation function.
  \param[in] plan - scalar plan ::kfft_plan_sclr pointer
  \param[in] fin - input ::kfft_cpx buffer (don't changed)
  \param[in] fout - output ::kfft_scalar buffer
  \result standart ::kfft_ret_flags return status

  \warning Function NOT control input and output buffer (such as the NULL, overflow, bad
  buffers-size etc.). Developer must control this manualy.
 */
KFFT_API kfft_return_t
kfft_evali_scalar(kfft_plan_sclr* plan, const kfft_cpx* fin, kfft_scalar* fout);

/*Internal use functions (NOT provide KFFT_API) for optimize extensions */

/*!
  \brief Internal usage function @see kfft_eval_scalar
  Process forward evaluation function.
  \param[in] plan - scalar plan ::kfft_plan_sclr pointer
  \param[in] fin - input ::kfft_scalar buffer (don't changed)
  \param[in] fout - output ::kfft_cpx buffer
  \param[in] ftmp - temporary ::kfft_cpx buffer

  \result standart ::kfft_ret_flags return status
*/
kfft_return_t
kfft_eval_scalar_internal(kfft_plan_sclr* plan, const kfft_scalar* fin, kfft_cpx* fout,
                          kfft_cpx* ftmp);
/*!
  \brief Internal usage function @see kfft_evali_scalar
  Process inverse evaluation function.
  \param[in] plan - scalar plan ::kfft_plan_sclr pointer
  \param[in] fin - input ::kfft_cpx buffer (don't changed)
  \param[in] fout - output ::kfft_scalar buffer
  \param[in] ftmp - temporary ::kfft_cpx buffer

  \result standart ::kfft_ret_flags return status
*/
kfft_return_t
kfft_evali_scalar_internal(kfft_plan_sclr* plan, const kfft_cpx* fin, kfft_scalar* fout,
                           kfft_cpx* ftmp);
/*!
  \brief Internal usage function @see kfft_eval_scalar_norm
  Process forward evaluation function with normalization.

  \param[in] plan - scalar plan ::kfft_plan_sclr pointer
  \param[in] fin - input ::kfft_scalar buffer (don't changed)
  \param[in] fout - output ::kfft_scalar buffer
  \param[in] ftmp - temporary ::kfft_cpx buffer

  \result standart ::kfft_ret_flags return status
*/
kfft_return_t
kfft_eval_scalar_norm_internal(kfft_plan_sclr* plan, const kfft_scalar* fin, kfft_scalar* fout,
                               kfft_cpx* ftmp);

#if defined(KFFT_DYNAPI_ENABLE)
// clang-format off
typedef kfft_plan_sclr*
(*kfft_callback_config_scalar)(const uint32_t nfft, const uint32_t flags,
                               kfft_pool_t* A, size_t* lenmem);
typedef kfft_return_t
(*kfft_callback_eval_scalar)(kfft_plan_sclr* plan, const kfft_scalar* fin, kfft_cpx* fout);
typedef kfft_return_t
(*kfft_callback_eval_scalar_norm)(kfft_plan_sclr* plan, const kfft_scalar* fin, kfft_scalar* fout);
typedef kfft_return_t
(*kfft_callback_evali_scalar)(kfft_plan_sclr* plan, const kfft_cpx* fin, kfft_scalar* fout);
// clang-format on
#endif /* KFFT_DYNAPI_ENABLE */
