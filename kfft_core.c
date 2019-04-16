#include "_kfft_guts.h"
#include "_kfft_bf.h"
#include "_kfft_gen.h"

#if defined (USE_RADER_ALGO)
static inline int
_kfr_ispower ( unsigned a,
               unsigned b,
               unsigned p )
{
    uint64_t res = 1;           // Initialize result
    uint64_t x = (uint64_t)a;   // Potential overflow fix
    uint64_t y = (uint64_t)b;   // Potential overflow fix

    x = x % p; // Update x if it is more than or equal to p

    while (y > 0)
    {
        // If y is odd, multiply x with result
        if (y & 1)
            res = (res * x) % p; // FIXME Maybe integer overflow
        // y must be even now
        y = y >> 1; // y = y/2
        x = (x * x) % p; // FIXME Maybe integer overflow
    }
    return (res != 1) ? 1 : 0;
}

// WARNING in kfft num - always prime. Don't check this
static unsigned
_kfr_find_root (unsigned num) {
    unsigned phi = num - 1;
    unsigned n = phi;

    unsigned primes [MAX_ROOTS];
    unsigned count = 0;

    for (unsigned i = 2; i*i <= n; i++) {
        if (n % i == 0){
            primes [count++] = i;
            while (n % i == 0)
                n /= i;
        }
    }
    if (n > 1)
        primes [ count++ ] = n;

    for (unsigned res=2; res<=num; ++res) {
		int ok = 1;
		for (unsigned i = 0; count > i && ok; ++i) {
			ok &= _kfr_ispower (res, phi / primes [i], num);
            kfft_trace ("%d\t%d\t%d\n", res, i, primes [i]);
        }
        if (ok)  return res;
	}
    return 0;
}
#endif /* USE_RADER_ALGO */

static
void kf_work(
        kiss_fft_cpx * Fout,
        const kiss_fft_cpx * f,
        const size_t fstride,
        int in_stride,
        int * factors,
        const __fft_cfg st
        )
{
    kiss_fft_cpx * Fout_beg=Fout;
    const int p=*factors++; /* the radix  */
    const int m=*factors++; /* stage's fft length/p */
    const kiss_fft_cpx * Fout_end = Fout + p*m;

#ifdef _OPENMP
    // use openmp extensions at the
    // top-level (not recursive)
    if (fstride==1 && p<=5)
    {
        int k;

        // execute the p different work units in different threads
#       pragma omp parallel for
        for (k=0;k<p;++k)
            kf_work( Fout +k*m, f+ fstride*in_stride*k,fstride*p,in_stride,factors,st);
        // all threads have joined by this point

        switch (p) {
            case 2: kf_bfly2(Fout,fstride,st,m); break;
            case 3: kf_bfly3(Fout,fstride,st,m); break;
            case 4: kf_bfly4(Fout,fstride,st,m); break;
            case 5: kf_bfly5(Fout,fstride,st,m); break;
            default: kf_bfly_generic(Fout,fstride,st,m,p);break;
        }
        return;
    }
#endif

    if (m==1) {
        do{
            *Fout = *f;
            f += fstride*in_stride;
        }while(++Fout != Fout_end );
    }else{
        do{
            // recursive call:
            // DFT of size m*p performed by doing
            // p instances of smaller DFTs of size m,
            // each one takes a decimated version of the input
            kf_work( Fout , f, fstride*p, in_stride, factors,st);
            f += fstride*in_stride;
        }while( (Fout += m) != Fout_end );
    }

    Fout=Fout_beg;

    // recombine the p smaller DFTs
    switch (p) {
        case 2: kf_bfly2(Fout,fstride,st,m); break;
        case 3: kf_bfly3(Fout,fstride,st,m); break;
        case 4: kf_bfly4(Fout,fstride,st,m); break;
        case 5: kf_bfly5(Fout,fstride,st,m); break;
        default: kf_bfly_generic(Fout,fstride,st,m,p);break;
    }
}

/*  facbuf is populated by p1,m1,p2,m2, ...
    where
    p[i] * m[i] = m[i-1]
    m0 = n                  */
static
void kf_factor(int n,int * facbuf)
{
    int p=4;
    double floor_sqrt;
    floor_sqrt = floor( sqrt((double)n) );

    /*factor out powers of 4, powers of 2, then any remaining primes */
    do {
        while (n % p) {
            switch (p) {
                case 4: p = 2; break;
                case 2: p = 3; break;
                default: p += 2; break;
            }
            if (p > floor_sqrt)
                p = n;          /* no more factors, skip to end */
        }
        n /= p;
        *facbuf++ = p;
        *facbuf++ = n;
    } while (n > 1);
}

static inline void
__kiss_fft_init (__fft_cfg   st,
                 int         nfft,
                 int         inverse_fft,
                 int         level)
{
        int i;
        st->nfft=nfft;
        st->inverse = inverse_fft;
        st->level   = level;

        kf_factor(nfft,st->factors);
#if defined (USE_RADER_ALGO)
        if (level > 0) {
            // TODO Rader initialization
            st->prime_root = _kfr_find_root (nfft);
        } else {
#endif /* USE_RADER_ALGO */
            for (i=0;i<nfft;++i) {
                const double pi=3.141592653589793238462643383279502884197169399375105820974944;
                double phase = -2*pi*i / nfft;
                if (st->inverse)
                    phase *= -1;
                kf_cexp(st->twiddles+i, phase );
            }
#if defined (USE_RADER_ALGO)
        }
#endif
}

/*
 *
 * User-callable function to allocate all necessary storage space for the fft.
 *
 * The return value is a contiguous block of memory, allocated with malloc.  As such,
 * It can be freed with free(), rather than a kiss_fft-specific function.
 * */
__fft_cfg
__kiss_fft_config (int nfft,
                   int inverse_fft,
                   int level,
                   void * mem,
                   size_t * lenmem )
{
    __fft_cfg st=NULL;
    size_t memneeded = sizeof(struct kiss_fft_state)
        + sizeof(kiss_fft_cpx)*(nfft-1); /* twiddle factors*/

    if ( lenmem==NULL) {
        st = ( __fft_cfg)KISS_FFT_MALLOC( memneeded );
    }else{
        if (mem != NULL && *lenmem >= memneeded)
            st = (__fft_cfg)mem;
        *lenmem = memneeded;
    }
    if (st) {
        __kiss_fft_init(st, nfft, inverse_fft, level);
    }
    return st;
}


static inline void
__kiss_fft_stride(__fft_cfg st,const kiss_fft_cpx *fin,kiss_fft_cpx *fout,int in_stride)
{
    if (fin == fout) {
        //NOTE: this is not really an in-place FFT algorithm.
        //It just performs an out-of-place FFT into a temp buffer
        kiss_fft_cpx * tmpbuf = (kiss_fft_cpx*)KISS_FFT_TMP_ALLOC( sizeof(kiss_fft_cpx)*st->nfft);
        kf_work(tmpbuf,fin,1,in_stride, st->factors,st);
        memcpy(fout,tmpbuf,sizeof(kiss_fft_cpx)*st->nfft);
        KISS_FFT_TMP_FREE(tmpbuf);
    }else{
        kf_work( fout, fin, 1,in_stride, st->factors,st );
    }
}

void
__kiss_fft(__fft_cfg cfg,const kiss_fft_cpx *fin,kiss_fft_cpx *fout)
{
    __kiss_fft_stride (cfg,fin,fout,1);
}
