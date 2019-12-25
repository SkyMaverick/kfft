#include "kfft_guts.h"

#include "kfft_generic.c"
#include "kfft_bfly.c"

#ifdef KFFT_RADER_ALGO
static inline unsigned
_kfr_power(unsigned x, unsigned y, unsigned m) {
    if (y == 0)
        return 1;
    unsigned p = _kfr_power(x, y / 2, m) % m;
    p = (p * p) % m;

    return (y % 2 == 0) ? p : (x * p) % m;
}

static unsigned
_kfr_gcd(unsigned a, unsigned b) {
    if (a == 0)
        return b;
    return _kfr_gcd(b % a, a);
}

// WARNING in kfft num - always prime. Don't check this
unsigned
kfft_prime_root(unsigned num) {
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
            if (_kfr_power(res, phi / primes[i], num) == 1)
                ok = 0;
        }
        if (ok)
            return res;
    }
    return 0;
}

unsigned
kfft_primei_root(unsigned a, unsigned m) {
    return (_kfr_gcd(a, m) != 1) ? 0 : _kfr_power(a, m - 2, m);
}

#endif /* KFFT_RADER_ALGO */

static void
kf_work(kfft_cpx* Fout, const kfft_cpx* f, const size_t fstride, int in_stride, int* factors,
        const kfft_kplan_t* st) {
    kfft_cpx* Fout_beg = Fout;
    const int p = *factors++; /* the radix  */
    const int m = *factors++; /* stage's fft length/p */
    const kfft_cpx* Fout_end = Fout + p * m;

    kfft_trace("Work: in - %4.1fi%4.1f | end -  %4.1fi%4.1f | f - %4.1fi%4.1f | fstride - %zu | "
               "in_stride - %d \n",
               Fout->r, Fout->i, Fout_end->r, Fout_end->i, f->r, f->i, fstride, in_stride);
    kfft_trace("      p - %d | m - %d\n", p, m);

    if (m == 1) {
        do {
            *Fout = *f;
            f += fstride * in_stride;
        } while (++Fout != Fout_end);
    } else {
        do {
            // recursive call:
            // DFT of size m*p performed by doing
            // p instances of smaller DFTs of size m,
            // each one takes a decimated version of the input
            kf_work(Fout, f, fstride * p, in_stride, factors, st);
            f += fstride * in_stride;
        } while ((Fout += m) != Fout_end);
    }

    Fout = Fout_beg;

    // recombine the p smaller DFTs
    switch (p) {
    case 2:
        kf_bfly2(Fout, fstride, st, m);
        break;
    case 3:
        kf_bfly3(Fout, fstride, st, m);
        break;
    case 4:
        kf_bfly4(Fout, fstride, st, m);
        break;
    case 5:
        kf_bfly5(Fout, fstride, st, m);
        break;
    default:
        kf_bfly_generic(Fout, fstride, st, m, p);
        break;
    }
}

/*  facbuf is populated by p1,m1,p2,m2, ...
    where
    p[i] * m[i] = m[i-1]
    m0 = n                  */
static void
kf_factor(int n, int* facbuf, unsigned* root_buf) {
    unsigned p = 4;
    double floor_sqrt;
    floor_sqrt = floor(sqrt((double)n));

    unsigned* rbufs = root_buf;

    /*factor out powers of 4, powers of 2, then any remaining primes */
    do {
        while (n % p) {
            switch (p) {
            case 4:
                p = 2;
                break;
            case 2:
                p = 3;
                break;
            default:
                p += 2;
                break;
            }
            if (p > floor_sqrt)
                p = n; /* no more factors, skip to end */
        }
#if defined(KFFT_RADER_ALGO) && !defined(KFFT_MEMLESS_MODE)
        if (p > MAX_BFLY_LEVEL && p > KFFT_RADER_LIMIT) {
            unsigned* root_buf = rbufs;
            while (*root_buf != p) {
                if (*root_buf == 0) {
                    *root_buf = p;

                    unsigned root = kfft_prime_root(p);
                    *(root_buf + 1) = root;
                    *(root_buf + 2) = kfft_primei_root(p, root);

                    break;
                }
                root_buf +=3;
            }
        }
#else
        (void)root_buf; // fake operation for disable warning
#endif /* RADER and not MEMLESS */
        n /= p;
        *facbuf++ = p;
        *facbuf++ = n;
    } while (n > 1);
}

static inline void
kfft_kinit(kfft_kplan_t* st, int nfft, int inverse_fft, int level) {
#if defined(KFFT_RADER_ALGO) && !defined(KFFT_MEMLESS_MODE)
    kf_factor(nfft, st->factors, st->roots);
#else
    kf_factor(nfft, st->factors, NULL);
#endif
    st->nfft = nfft;
    st->inverse = inverse_fft;
    st->level = level;

#ifndef KFFT_MEMLESS_MODE
    for (int i = 0; i < nfft; ++i) {
        kfft_scalar phase = -2 * KFFT_CONST_PI * i / nfft;
        if (st->inverse)
            phase *= -1;
        kf_cexp(st->twiddles + i, phase);
    }
#endif /* memless */
}

/*
 *
 * User-callable function to allocate all necessary storage space for the fft.
 *
 * The return value is a contiguous block of memory, allocated with malloc.  As such,
 * It can be freed with free(), rather than a kfft-specific function.
 * */
kfft_kplan_t*
kfft_kconfig(int nfft, int inverse_fft, int level, void* mem, size_t* lenmem) {
    kfft_kplan_t* st = NULL;
#ifndef KFFT_MEMLESS_MODE
    size_t memneeded = sizeof(struct kfft_kstate) + sizeof(kfft_cpx) * (nfft); /* twiddle factors*/
#else
    size_t memneeded = sizeof(struct kfft_kstate); /* twiddle factors*/
#endif /* memless */

    if (lenmem == NULL) {
        st = (kfft_kplan_t*)KFFT_MALLOC(memneeded);
        kfft_trace("Alloc kernel plan: %p. Memneed: %zu\n", (void*)st, memneeded);
    } else {
        if (mem != NULL && *lenmem >= memneeded) {
            st = (kfft_kplan_t*)mem;
            kfft_trace("Redefine kernel plan: %p. Memneed: %zu\n", (void*)st, memneeded);
        }
        *lenmem = memneeded;
    }
    if (st) {
        kfft_kinit(st, nfft, inverse_fft, level);

        kfft_trace("Config: nfft - %d | inv - %d | lvl - %d", st->nfft, st->inverse, st->level);
        kfft_trace("        memory - %zu\n", memneeded);
    }
    return st;
}

static inline void
kfft_kstride(kfft_kplan_t* st, const kfft_cpx* fin, kfft_cpx* fout, int in_stride) {
    if (fin == fout) {
        // NOTE: this is not really an in-place FFT algorithm.
        // It just performs an out-of-place FFT into a temp buffer
        kfft_cpx* tmpbuf = (kfft_cpx*)KFFT_TMP_ALLOC(sizeof(kfft_cpx) * st->nfft);

        kfft_trace("ALLOC temp buffer: %p\n", (void*)tmpbuf);
        kf_work(tmpbuf, fin, 1, in_stride, st->factors, st);
        memcpy(fout, tmpbuf, sizeof(kfft_cpx) * st->nfft);
        kfft_trace("FREE temp buffer: %p\n", (void*)tmpbuf);

        KFFT_TMP_FREE(tmpbuf);
    } else {
        kf_work(fout, fin, 1, in_stride, st->factors, st);
    }
}

void
__kfft(kfft_kplan_t* cfg, const kfft_cpx* fin, kfft_cpx* fout) {
    kfft_kstride(cfg, fin, fout, 1);
}
