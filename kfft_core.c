#include "_kfft_guts.h"
#include "_kfft_bf.h"
#include "_kfft_gen.h"

static void
kf_work(kfft_cpx* Fout, const kfft_cpx* f, const size_t fstride, int in_stride, int* factors,
        const kfft_kplan_p st) {
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
kf_factor(int n, int* facbuf) {
    int p = 4;
    double floor_sqrt;
    floor_sqrt = floor(sqrt((double)n));

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
        n /= p;
        *facbuf++ = p;
        *facbuf++ = n;
    } while (n > 1);
}

static inline void
kfft_kinit(kfft_kplan_p st, int nfft, int inverse_fft, int level) {
    int i;
    kf_factor(nfft, st->factors);

    st->nfft = nfft;
    st->inverse = inverse_fft;
    st->level = level;

    for (i = 0; i < nfft; ++i) {
        const kfft_scalar pi = KFFT_CONST_PI;
        kfft_scalar phase = -2 * pi * i / nfft;
        if (st->inverse)
            phase *= -1;
        kf_cexp(st->twiddles + i, phase);
    }
}

/*
 *
 * User-callable function to allocate all necessary storage space for the fft.
 *
 * The return value is a contiguous block of memory, allocated with malloc.  As such,
 * It can be freed with free(), rather than a kfft-specific function.
 * */
kfft_kplan_p
kfft_kconfig(int nfft, int inverse_fft, int level, void* mem, size_t* lenmem) {
    kfft_kplan_p st = NULL;
    size_t memneeded = sizeof(struct kfft_state) + sizeof(kfft_cpx) * (nfft); /* twiddle factors*/

    if (lenmem == NULL) {
        st = (kfft_kplan_p)KFFT_MALLOC(memneeded);
    } else {
        if (mem != NULL && *lenmem >= memneeded)
            st = (kfft_kplan_p)mem;
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
kfft_kstride(kfft_kplan_p st, const kfft_cpx* fin, kfft_cpx* fout, int in_stride) {
    if (fin == fout) {
        // NOTE: this is not really an in-place FFT algorithm.
        // It just performs an out-of-place FFT into a temp buffer
        kfft_cpx* tmpbuf = (kfft_cpx*)KFFT_TMP_ALLOC(sizeof(kfft_cpx) * st->nfft);

        printf("ALLOC temp buffer: %p\n", tmpbuf);
        kf_work(tmpbuf, fin, 1, in_stride, st->factors, st);
        memcpy(fout, tmpbuf, sizeof(kfft_cpx) * st->nfft);
        printf("FREE temp buffer: %p\n", tmpbuf);

        KFFT_TMP_FREE(tmpbuf);
    } else {
        kf_work(fout, fin, 1, in_stride, st->factors, st);
    }
}

void
__kfft(kfft_kplan_p cfg, const kfft_cpx* fin, kfft_cpx* fout) {
    kfft_kstride(cfg, fin, fout, 1);
}
