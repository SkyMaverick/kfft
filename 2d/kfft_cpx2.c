#include "kfft.h"

// clang-format off
#define kfft_trace_2d(level, fmt, ...)                                                           \
    kfft_trace("[CORE (L%d)]"" " fmt,level, __VA_ARGS__)
// clang-format on

KFFT_API kfft_comp_t*
kfft_config_cpx2(const uint32_t x_size, const uint32_t y_size, const uint32_t flags, kfft_pool_t* A,
                 size_t* lenmem) {
    return NULL;
}

KFFT_API kfft_return_t
kfft_eval_cpx2(kfft_comp2_t* cfg, const kfft_cpx* fin, kfft_cpx* fout) {
    return KFFT_RET_SUCCESS;
}

#undef kfft_trace_2d
