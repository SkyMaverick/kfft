#pragma once

/*!
    \file
    \brief 1D complex evaluation functions

    Complex 1D-buffer evaluations functionality.
    This is basis fuctions for all library functionality.
 */
// clang-format off
#define kfft_trace_core(level, fmt, ...)                                                           \
    kfft_trace("[CORE (L%d)]"" " fmt,level, __VA_ARGS__)
// clang-format on

/// Maximum number of divisors
#define MAX_FACTORS 32
/// Maximum primitive roots by modulo for Rader algo
#define MAX_ROOTS 32

/// plan object abstraction and unification structure
typedef struct {
    kfft_pool_t* mmgr; ///< allocator for internal-plan fast memory allocations
} kfft_object_t;

/*!
  Internal usage smallest plan type for prime-lenght sequence evaluation
  with Rader algorithm. Use as part ::kfft_plan_cpx only.
  */
typedef struct {
    uint32_t prime; ///< plan lenght
    uint32_t p, q;  ///< forward and multiplicative inverse group generators

    uint32_t* qidx; ///< forward group indexes map
    uint32_t* pidx; ///< inverse group indexes map

    kfft_cpx* shuffle_twiddles; ///< reshuffle twiddles buffer

    struct kfft_kstate* plan;     ///< plan for recursive evaluate prime - 1 sequense
    struct kfft_kstate* plan_inv; ///< plan for evaluate result sequense
} kfft_plan_rader;

/*!
    1D complex operations plan type.
 */
typedef struct kfft_kstate {
    kfft_object_t object; ///< standart object structure

    uint32_t nfft; ///< sequence lenght (count)
    uint8_t level;
    /*!<
        \brief Plan recursive level.
        \note Because libkfft use reqursive evaluation we must separate top-level plan
        and sub-level plan we must separate so as not to get confused with the allocation
        of internal memory (for example)
     */
    uint32_t flags; ///< operation flags ::kfft_eval_flags

    uint8_t fac_count;                 ///< factors count
    uint32_t factors[2 * MAX_FACTORS]; ///< factor map

    uint8_t prm_count;                 ///< prime roots by modulo (PRbM) count
    kfft_plan_rader primes[MAX_ROOTS]; ///< PRbM map

    kfft_cpx* twiddles; ///< twiddles buffer
} kfft_plan_cpx;

/*!
  \param[in] nfft - lenght input sequense
  \param[in] flags - operation flags
  \param[in] A - plan ainternal allocator structure (if need use KFFT_PLAN_MMGR macro)
  \param[in] lenmem - vaiable for memory get pointer
  \return complex plan structure
 */
KFFT_API kfft_plan_cpx*
kfft_config_cpx(const uint32_t nfft, const uint32_t flags, kfft_pool_t* A, size_t* lenmem);
/*!
 */
KFFT_API kfft_return_t
kfft_eval_cpx(kfft_plan_cpx* plan, const kfft_cpx* fin, kfft_cpx* fout);

#if defined(KFFT_DYNAPI_ENABLE)
// clang-format off
typedef kfft_plan_cpx*
(*kfft_callback_config_cpx)(const uint32_t nfft, const uint32_t flags, kfft_pool_t* A, size_t* lenmem);
typedef kfft_return_t
(*kfft_callback_eval_cpx)(kfft_plan_cpx* plan, const kfft_cpx* fin, kfft_cpx* fout);
// clang-format on
#endif /* KFFT_DYNAPI_ENABLE */
