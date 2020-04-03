#pragma once

typedef struct {
    kfft_object_t object;

    uint32_t nfft, x, y;
    uint32_t flags;

    kfft_comp_t* plan_x;
    kfft_comp_t* plan_y;
} kfft_comp2_t;

KFFT_API kfft_comp2_t*
kfft_config_cpx2(const uint32_t x_size, const uint32_t y_size, const uint32_t flags, kfft_pool_t* A,
                 size_t* lenmem);
KFFT_API kfft_return_t
kfft_eval_cpx2(kfft_comp2_t* cfg, const kfft_cpx* fin, kfft_cpx* fout);
