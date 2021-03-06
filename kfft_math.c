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
    return (kfft_math_gcd(a, m) != 1) ? 0 : kfft_math_modpow(m, a - 2, a);
}

#endif /* KFFT_RADER_ALGO */

KFFT_API void
kfft_math_hadamard_cpx(kfft_cpx* Fout, const kfft_cpx* Fin, uint32_t size) {
    kfft_cpx tmp;
    KFFT_OMP( omp parallel for schedule(static) private(tmp))
    for (uint32_t i = 0; i < size; i++) {
        C_MUL(tmp, Fout[i], Fin[i]);
        C_CPY(Fout[i], tmp);
    }
}

KFFT_API void
kfft_math_hadamard_scalar(kfft_scalar* Fout, const kfft_scalar* Fin, uint32_t size) {
    KFFT_OMP( omp parallel for schedule(static))
    for (uint32_t i = 0; i < size; i++) {
        Fout[i] *= Fin[i];
    }
}

KFFT_API void
kfft_math_transpose_cpx(const kfft_cpx* Fin, kfft_cpx* Fout, const uint32_t x, const uint32_t y) {
    KFFT_OMP( omp parallel for schedule(static))
    for (uint64_t n = 0; n < x * y; n++) {
        uint64_t i = n / y;
        uint64_t j = n % y;
        C_CPY(Fout[n], Fin[x * j + i]);
    }
}
KFFT_API void
kfft_math_transpose_scalar(const kfft_scalar* Fin, kfft_scalar* Fout, const uint32_t x,
                           const uint32_t y) {
    KFFT_OMP( omp parallel for schedule(static))
    for (uint64_t n = 0; n < x * y; n++) {
        uint64_t i = n / y;
        uint64_t j = n % y;
        Fout[n] = Fin[x * j + i];
    }
}

KFFT_API void
kfft_math_transpose_ip_cpx(kfft_cpx* Fin, const uint32_t x, const uint32_t y) {
    uint32_t r = y;
    uint32_t c = x;

    if (r <= 1 || c <= 1)
        return;

    for (uint32_t ind_last = r * c - 2; c > 1; --ind_last, --c) {
        for (int64_t i = r - 2; i >= 0; --i) {
            uint64_t ind = i * c + (c - 1);
            kfft_cpx buf = {0, 0};
            C_CPY(buf, Fin[ind]);
            while (ind < ind_last) {
                C_CPY(Fin[ind], Fin[ind + 1]);
                ++ind;
            }
            C_CPY(Fin[ind_last], buf);
            --ind_last;
        }
    }
}
KFFT_API void
kfft_math_transpose_ip_scalar(kfft_scalar* Fin, const uint32_t x, const uint32_t y) {
    uint32_t c = x, r = y;

    if (r <= 1 || c <= 1)
        return;

    for (uint32_t ind_last = r * c - 2; c > 1; --ind_last, --c) {
        for (int64_t i = r - 2; i >= 0; --i) {
            uint64_t ind = i * c + (c - 1);
            kfft_scalar buf = Fin[ind];
            while (ind < ind_last) {
                Fin[ind] = Fin[ind + 1];
                ++ind;
            }
            Fin[ind_last] = buf;
            --ind_last;
        }
    }
}

KFFT_API void
kfft_math_magnitude(const kfft_cpx* Fin, kfft_scalar* Fout, uint32_t size) {
    KFFT_OMP( omp parallel for schedule(static))
    for (uint32_t i = 0; i < size; i++) {
        Fout[i] = kfft_math_mgnt(&Fin[i]);
    }
}

KFFT_API void
kfft_math_magnitude_ip(kfft_cpx* Fin, uint32_t size) {
    KFFT_OMP( omp parallel for schedule(static))
    for (uint32_t i = 0; i < size; i++) {
        ((kfft_scalar*)Fin)[i] = kfft_math_mgnt(&Fin[i]);
    }
}

KFFT_API kfft_scalar
kfft_math_average(const kfft_scalar* Fin, uint32_t size) {
    kfft_scalar ret = 0;

    KFFT_OMP( omp parallel for schedule(static) shared(ret))
    for (uint32_t i = 0; i < size; i++) {
        ret += Fin[i];
    }
    return ret / size;
}
