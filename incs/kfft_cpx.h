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
  \brief Complex plan create and configuration function

  Function create new, configure, reconfigure exists complex plan.
  It's universal function for prepare fft process operation.

  \param[in] nfft - lenght input sequense
  \param[in] flags - operation flags
  \param[in] A - plan a internal allocator structure (if need use KFFT_PLAN_MMGR macro)
  \param[in] lenmem - vaiable for memory get pointer
  \return complex plan structure or NULL

  \note Typical function usage:
    - kfft_config_cpx(<lenght>, <flags>, NULL, NULL) : create new standart complex plan (return
  ::kfft_plan_cpx or NULL if error)
    - kfft_config_cpx(<lenght>, <flags>, NULL, <valptr>) : return memory needed size in valptr value
  (return NULL)
    - kfft_config_cpx(<lenght>, <flags>, <MMGR>, NULL) : create new complex plan object
  ::kfft_plan_cpx into MMGR (return ::kfft_plan_cpx or NULL if error)
    - kfft_config_cpx(<lenght>, <flags>, <MMGR>, <valptr>)
        - if *valptr value >= needed memory: create new complex plan object ::kfft_plan_cpx into
  MMGR (return ::kfft_plan_cpx or NULL if error and needed memory in *valptr value)
        - if *valptr value < needed memory: return NULL and needed memory size as *valptr value
        - if flags ::KFFT_FLAG_RENEW enabled: clear input MMGR and create new plan if the conditions
  for *valptr

  \warning Function may be allocate memory. Use ::kfft_cleanup for correctly clean and free this
  memory.
 */
KFFT_API kfft_plan_cpx*
kfft_config_cpx(const uint32_t nfft, const uint32_t flags, kfft_pool_t* A, size_t* lenmem);
/*!
  \brief Process evaluation function.
  \param[in] plan - complex plan ::kfft_plan_cpx pointer
  \param[in] fin - input ::kfft_cpx buffer (don't changed)
  \param[in] fout - output ::kfft_cpx buffer
  \result standart ::kfft_ret_flags return status

  \note Is a valide using fin == fout. Despite the fact that the function is NOT in-place
  evaluation, it can create a temporary buffer for the operation.
  \warning Function NOT control input and output buffer (such as the NULL, overflow, bad
  buffers-size etc.). Developer must control this manualy.
 */
KFFT_API kfft_return_t
kfft_eval_cpx(kfft_plan_cpx* plan, const kfft_cpx* fin, kfft_cpx* fout);

/*!
  Macro to get the memory manager of a plan object
  \param[in] X - plan object
  \return plan memory manager pointer (kfft_pool_t*)
 */
#define KFFT_PLAN_MMGR(X) (*((kfft_pool_t**)(X)))

/*!
  Macro to get the memory manager of a plan object if it's enabled
  or NULL if plan is NULL. Required for the assignment operation.
  \param[in] X - plan object
  \return plan memory manager pointer (kfft_pool_t*) or NULL
 */
#define KFFT_PLAN_MMGR_NULL(X) ((X != NULL) ? KFFT_PLAN_MMGR((X)) : NULL)

/*!
  Macro to get the align information of a plan object
  \param[in] X - plan object
  \return uint8_t number memory align for plan (0, 16, 32)
 */
#define KFFT_PLAN_ALIGN(X) ((X != NULL) ? KFFT_PLAN_MMGR((X))->align : 0)

/*!
  Macro to get the acceleration info of a plan object

  \param[in] X - plan object
  \warning argument (X) mustn't NULL

  \return kfft_simd_t vector extension information
 */
#define KFFT_PLAN_VEX(X) KFFT_PLAN_MMGR((X))->vex

#if defined(KFFT_DYNAPI_ENABLE)
// clang-format off
typedef kfft_plan_cpx*
(*kfft_callback_config_cpx)(const uint32_t nfft, const uint32_t flags, kfft_pool_t* A, size_t* lenmem);
typedef kfft_return_t
(*kfft_callback_eval_cpx)(kfft_plan_cpx* plan, const kfft_cpx* fin, kfft_cpx* fout);
// clang-format on
#endif /* KFFT_DYNAPI_ENABLE */
