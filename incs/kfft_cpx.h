#pragma once

// clang-format off
#define kfft_trace_core(level, fmt, ...)                                                           \
    kfft_trace("[CORE (L%d)]"" " fmt,level, __VA_ARGS__)
// clang-format on

#define MAX_FACTORS 32
#define MAX_ROOTS 32

typedef struct {
    kfft_pool_t* mmgr;
} kfft_object_t;

typedef struct {
    uint32_t prime, p, q;

    uint32_t* qidx;
    uint32_t* pidx;

    kfft_cpx* shuffle_twiddles;

    struct kfft_kstate* splan;
    struct kfft_kstate* splani;
} kfft_splan_t;

typedef struct kfft_kstate {
    kfft_object_t object;

    uint32_t nfft;

    uint32_t level;
    uint32_t flags;

    uint8_t fac_count;                 // factors count
    uint32_t factors[2 * MAX_FACTORS]; // factor values

    uint8_t prm_count;
    kfft_splan_t primes[MAX_ROOTS]; //

    kfft_cpx* twiddles; // twiddles
} kfft_comp_t;

static inline uint32_t
_kfr_power(uint32_t x, uint32_t y, uint32_t m) {
    if (y == 0)
        return 1;
    uint64_t p = _kfr_power(x, y / 2, m) % m;
    p = (p * p) % m;

    return (y % 2 == 0) ? (uint32_t)p : (uint32_t)((x * p) % m);
}

KFFT_API kfft_comp_t*
kfft_config_cpx(const uint32_t nfft, const uint32_t flags, const uint8_t level, kfft_pool_t* A,
                size_t* lenmem);

KFFT_API kfft_return_t
kfft_eval_cpx(kfft_comp_t* cfg, const kfft_cpx* fin, kfft_cpx* fout);

KFFT_API kfft_return_t
kfft_convolution(kfft_cpx* Fout, kfft_cpx* Fin, kfft_comp_t* P, kfft_comp_t* Pi);

#define kfft_kfree(X)                                                                              \
    do {                                                                                           \
        free(X);                                                                                   \
        X = NULL;                                                                                  \
    } while (0)
