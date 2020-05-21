#pragma once

KFFT_API void
kfft_info(kfft_info_t* info);
KFFT_API uint32_t
kfft_next_fast_size(uint32_t n);
KFFT_API void*
kfft_malloc(uint32_t sz);
KFFT_API void
kfft_free_null(void** mem);
KFFT_API kfft_return_t
kfft_cleanup(void* mem);
KFFT_API const char*
kfft_strerr(const kfft_return_t code);
