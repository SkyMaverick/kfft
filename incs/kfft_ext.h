#pragma once

KFFT_API void
kfft_info(kfft_info_t* info);
KFFT_API uint32_t
kfft_next_fast_size(uint32_t n);
KFFT_API kfft_return_t
kfft_cleanup(uintptr_t mem);
KFFT_API const char*
kfft_strerr(const kfft_return_t code);
