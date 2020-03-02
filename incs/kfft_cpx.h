#pragma once

#define MAX_FACTORS 32
#define MAX_ROOTS 32

typedef struct {
    uint32_t prime, p, q;

    uint32_t* ridx;
    uint32_t* rtidx;

    struct kfft_kstate* splan;
} kfft_splan_t;

typedef struct kfft_kstate {
    kfft_pool_t* mmgr;

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

KFFT_API void
kfft_eval_cpx(kfft_comp_t* cfg, const kfft_cpx* fin, kfft_cpx* fout);

#define kfft_kfree(X)                                                                              \
    do {                                                                                           \
        free(X);                                                                                   \
        X = NULL;                                                                                  \
    } while (0)
