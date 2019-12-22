#include "kfft_guts.h"

#include "kfft_core.h"
#include "kfft_rader.h"

static inline int
_kfr_ispower(unsigned a, unsigned b, unsigned p) {
    uint64_t res = 1;         // Initialize result
    uint64_t x = (uint64_t)a; // Potential overflow fix
    uint64_t y = (uint64_t)b; // Potential overflow fix

    x = x % p; // Update x if it is more than or equal to p

    while (y > 0) {
        // If y is odd, multiply x with result
        if (y & 1)
            res = (res * x) % p; // FIXME Maybe integer overflow
        // y must be even now
        y = y >> 1;      // y = y/2
        x = (x * x) % p; // FIXME Maybe integer overflow
    }
    return (res != 1) ? 1 : 0;
}

// WARNING in kfft num - always prime. Don't check this
static unsigned
_kfr_find_root(unsigned num) {
    unsigned phi = num - 1;
    unsigned n = phi;

    unsigned primes[MAX_ROOTS];
    unsigned count = 0;

    for (unsigned i = 2; i * i <= n; i++) {
        if (n % i == 0) {
            primes[count++] = i;
            while (n % i == 0)
                n /= i;
        }
    }
    if (n > 1)
        primes[count++] = n;

    for (unsigned res = 2; res <= num; ++res) {
        int ok = 1;
        for (unsigned i = 0; count > i && ok; ++i) {
            ok &= _kfr_ispower(res, phi / primes[i], num);
        }
        if (ok)
            return res;
    }
    return 0;
}

kfft_rplan
kfft_rconfig(int nfft, int inverse_fft, int level, void* mem, size_t* lenmem) {
    kfft_trace("%s: %d\n", "[RADER] Create RADER plan level", level);
    kfft_rplan st = NULL;
    size_t subsize = 0, memneeded = 0;

    kfft_kconfig(nfft, inverse_fft, level, NULL, &subsize);
    memneeded = sizeof(kfft_rplan_t) * 2 + subsize + sizeof(kfft_cpx) * (nfft * 3 / 2);

    if (lenmem == NULL) {
        st = (kfft_rplan)KFFT_MALLOC(memneeded);
    } else {
        if (*lenmem >= memneeded) {
            st = (kfft_rplan)mem;
        }
        *lenmem = memneeded;
    }
    if (!st)
        return NULL;

    st->substate = (kfft_kplan_t*)(st + 1);
    kfft_kconfig(nfft, inverse_fft, level, st->substate, &subsize);

#if defined(TRACE)
    kfft_trace("[RADER] %s: %zu\n", "Memory allocate", memneeded);
    kfft_trace("[RADER] %s: ", "Factors");
    for (int i = 0; st->substate->factors[i] != 0; i++) {
        kfft_trace("%d ", st->substate->factors[i]);
    }
    kfft_trace("%s\n", "");
#endif

    return st;
}

void
kfft_rfree(kfft_rplan* cfg) {
    if (cfg && *cfg) {
        free(*cfg);
        *cfg = NULL;
    }
}
