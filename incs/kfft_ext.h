#pragma once

/*!
    \file
    \brief Service functions API

    Service functions wich don't realize evaluation or math functions,
    but needed for library functionality.
 */

/*!
    Get all build and runtime information of library
    \param[in] info - pointer of ::kfft_info_t structure
    \result None
 */
KFFT_API void
kfft_info(kfft_info_t* info);
/*!
    Get optimal lenght sequense (>= N) for faster evaluation.
    \param[in] n - sequense lenght
    \return - next fast sequense lenght
 */
KFFT_API uint32_t
kfft_next_fast_size(uint32_t n);
/*!
    Allocation memory for evaluation data. Optimized analog malloc().
    \param[in] sz - size of allocated buffer (bytes)
    \return - allocated memory pointer or NULL if error.

    \note Use this feature for custom applications (replace malloc())
    because it respects memory alignment.
    \warning Don't use this function if you create libkfft extensions or internal structures.
    Use ::KFFT_MALLOC (or ::KFFT_TMP_ALLOC for temporary buffer) for it.

 */
KFFT_API void*
kfft_malloc(uint32_t sz);
/*!
    Freeing memory block what did allocate with ::kfft_malloc function. NULL input variable.
    \param[in] mem - pointer memory pointer
    \return None
 */
KFFT_API void
kfft_free_null(void** mem);
/// Service macro for ::kfft_free_null
#define kfft_free(X) (kfft_free_null((void**)(&(X))))
/*!
    Cleanup and free plan allocation after kfft_config_* function.
    \param[in] mem - memory area pointer
    \return kfft operation return status ::kfft_ret_flags
 */
KFFT_API kfft_return_t
kfft_cleanup(void* mem);
/// Service macro for ::kfft_cleanup
#define kfft_release(X)                                                                            \
    do {                                                                                           \
        if (kfft_cleanup((void*)(X)) == KFFT_RET_SUCCESS)                                          \
            (X) = NULL;                                                                            \
    } while (0)
/*!
    Return text interpratation ::kfft_ret_flags return code.
    \param[in] code - return code
    \return Constant message string
 */
KFFT_API const char*
kfft_strerr(const kfft_return_t code);

#if defined(KFFT_DYNAPI_ENABLE)
// clang-format off
typedef void
(*kfft_callback_info)(kfft_info_t* info);
typedef uint32_t
(*kfft_callback_next_fast_size)(uint32_t n);
typedef void*
(*kfft_callback_malloc)(uint32_t sz);
typedef void
(*kfft_callback_free_null)(void** mem);
typedef kfft_return_t
(*kfft_callback_cleanup)(void* mem);
typedef const char*
(*kfft_callback_strerr)(const kfft_return_t code);
// clang-format on
#endif /* KFFT_DYNAPI_ENABLE */
