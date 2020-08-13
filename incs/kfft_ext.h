#pragma once

KFFT_API void
kfft_info(kfft_info_t* info);
KFFT_API uint32_t
kfft_next_fast_size(uint32_t n);
KFFT_API void*
kfft_malloc(uint32_t sz);
KFFT_API void
kfft_free_null(void** mem);
#define kfft_free(X) (kfft_free_null((void**)(&(X))))
KFFT_API kfft_return_t
kfft_cleanup(void* mem);
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
