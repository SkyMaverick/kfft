#include "kfft.h"

#if defined(KFFT_RADER_ALGO)
// WARNING in kfft num - always prime. Don't check this
uint32_t
kfft_math_prmn(uint32_t num) {
    uint32_t phi = num - 1;
    uint32_t n = phi;

    uint32_t primes[MAX_ROOTS];
    uint32_t count = 0;

    for (uint32_t i = 2; i * i <= n; i++) {
        if (n % i == 0) {
            primes[count++] = i;
            while (n % i == 0)
                n /= i;
        }
    }
    if (n > 1)
        primes[count++] = n;

    for (uint32_t res = 2; res <= num; ++res) {
        bool ok = true;
        for (uint32_t i = 0; count > i && ok; ++i) {
            if (kfft_math_modpow(res, phi / primes[i], num) == 1)
                ok = false;
        }
        if (ok)
            return res;
    }
    return 0;
}

uint32_t
kfft_math_prmni(uint32_t a, uint32_t m) {
    return (kfft_math_gcd(a, m) != 1) ? 0 : kfft_math_modpow(a, m - 2, m);
}

#endif /* KFFT_RADER_ALGO */

void
kfft_math_adamar_cpx(kfft_cpx* Fout, kfft_cpx* Fin, uint32_t size) {
    // FIXME Now primitive algorithm
    kfft_cpx tmp;
    for (uint32_t i = 0; i < size; i++) {
        C_MUL(tmp, Fout[i], Fin[i]);
        C_CPY(Fout[i], tmp);
    }
}
