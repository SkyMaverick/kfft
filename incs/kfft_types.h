#pragma once

/*!
  \file
  \brief Base primitives and enums defined here

  Standart module defining fundamental types, constants and enumerations
 */

#include <inttypes.h>

#if defined(KFFT_OS_WINDOWS)
    #include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#endif

/*!
  Defined base scalar value type
  (64-bit defualt, 32-bit - KFFT_HALF_SCALAR option)
  */
#if defined(KFFT_HALF_SCALAR)
typedef float kfft_scalar;
#else
typedef double kfft_scalar;
#endif

/// Defined base complex value type
typedef struct {
    kfft_scalar r; ///< real part
    kfft_scalar i; ///< imaginary part
} kfft_cpx;

/// Build info bit-mask flags. @see kfft_info - function @see kfft_info_t - type
enum kfft_info_flags {
    KFFT_INFO_TRACE = 1 << 0,        ///< Trace info enabled
    KFFT_INFO_USE_SIMD = 1 << 1,     ///< SIMD extensions enabled
    KFFT_INFO_USE_ALLOCA = 1 << 2,   ///< Use alloca for temporary buffers
    KFFT_INFO_USE_SYSMATH = 1 << 3,  ///< Use system math (0 - internal)
    KFFT_INFO_USE_OPENMP = 1 << 4,   ///< Use OpenMP parallelization
    KFFT_INFO_RADER_ALGO = 1 << 5,   ///< Use Rader or generic algorithm
    KFFT_INFO_MEMLESS_MODE = 1 << 6, ///< Use lesser memory mode (very slowly)
    KFFT_INFO_HALF_SCALAR = 1 << 7,  ///< Use 32-bit scalar
};

/// Configuration plan modificators. @see kfft_eval_cpx
enum kfft_eval_flags {
    KFFT_FLAG_NORMAL = 0,       ///< Standart evaluation flag
    KFFT_FLAG_INVERSE = 1 << 0, ///< Inverse evaluation flag
    KFFT_FLAG_RENEW = 1 << 1,
    /*!< Reuse plan object for new plan if the size allows it
      \warning Dangerous flag blocking for use by nested structures with ::KFFT_CHECK_FLAGS
     */
    KFFT_FLAG_GENERIC = 1 << 2,       ///< Manual use of the generic algorithm for prime-size blocks
    KFFT_FLAG_GENERIC_ONLY = 1 << 3,  ///< Use ONLY generic algorithm for all sequence
    KFFT_FLAG_EXPAND_SCALAR = 1 << 4, ///< expand scalar buffer to full lenght (1D scalar only)
    KFFT_FLAG_DISABLE_NORM = 1 << 5,  ///< disable 1/N normalization for inverse transform
};

/*!
  Macro to protecting nested plans from destructive operations
  \param[in] X - plan object
  \return None
 */

#define KFFT_CHECK_FLAGS(X) ((X) & (~KFFT_FLAG_RENEW))


/*!
  \brief Library build information structure

  Structure is used to pass information about the library assembly configuration
  @see kfft_info()
 */
typedef struct {
    uint16_t vmajor; ///< Version info major
    uint16_t vminor; ///< Version info minor
    uint16_t vpatch; ///< Version info patch
    uint16_t flags;  ///< Configuration word. @see kfft_info_flags
} kfft_info_t;

/// Functions return information. @see kfft_strerr
enum kfft_ret_flags {
    KFFT_RET_SUCCESS = 0x0000,       ///< Normal return
    KFFT_RET_ALLOC_FAIL = 0x0001,    /*!<
           Return this code if internal allocator return fail
           \note Return this if kfft_config_* function return NULL */
    KFFT_RET_BUFFER_FAIL = 0x0002,   ///< Return this if buffer allocation is fail
    KFFT_RET_FREE_NULL = 0x0003,     ///< Return this if ::kfft_free_null give NULL argument
    KFFT_RET_IMPROPER_PLAN = 0x0004, /*!<
        Return this if use ::kfft_eval_scalar (or similar) with KFFT_FLAG_INVERSE plan.
        \note Use kfft_evali_scalar or kfft_eval_norm (or similar) for inverse plan.
        Remake plan if necessary. Use ::kfft_config_scalar() with ::KFFT_FLAG_RENEW*/
    KFFT_RET_BADARGUMENTS = 0x0005   ///< Return this if configuration arguments is bad
};

/// Type to define the return value
typedef unsigned kfft_return_t;

/*!
    This structure containe hardware VEX acceleration info
  */
typedef struct {
    uint8_t arch; ///< Architecture ID
    uint32_t ext; ///< HW extensions extensionse (with operation system correct)
} kfft_simd_t;
