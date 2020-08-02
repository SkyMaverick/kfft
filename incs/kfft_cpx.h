#pragma once

// clang-format off
#define kfft_trace_core(level, fmt, ...)                                                           \
    kfft_trace("[CORE (L%d)]"" " fmt,level, __VA_ARGS__)
// clang-format on

#define MAX_FACTORS 32
#define MAX_ROOTS 32

typedef struct {
    kfft_pool_t* mmgr; // allocator for internal-plan fast memory allocations
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

    uint8_t level;
    uint32_t flags;

    uint8_t fac_count;                 // factors count
    uint32_t factors[2 * MAX_FACTORS]; // factor values

    uint8_t prm_count;
    kfft_splan_t primes[MAX_ROOTS]; //

    kfft_cpx* twiddles; // twiddles
} kfft_plan_cpx;

KFFT_API kfft_plan_cpx*
kfft_config_cpx(const uint32_t nfft, const uint32_t flags, kfft_pool_t* A, size_t* lenmem);
KFFT_API kfft_return_t
kfft_eval_cpx(kfft_plan_cpx* cfg, const kfft_cpx* fin, kfft_cpx* fout);

#if defined(KFFT_DYNAPI_ENABLE)
// clang-format off
typedef kfft_plan_cpx*
(*kfft_callback_config_cpx)(const uint32_t nfft, const uint32_t flags, kfft_pool_t* A, size_t* lenmem);
typedef kfft_return_t
(*kfft_callback_eval_cpx)(kfft_plan_cpx* cfg, const kfft_cpx* fin, kfft_cpx* fout);
// clang-format on
#endif /* KFFT_DYNAPI_ENABLE */
