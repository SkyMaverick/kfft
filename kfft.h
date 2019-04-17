#ifndef KFFT_H
#define KFFT_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 ATTENTION!
 If you would like a :
 -- a utility that will handle the caching of fft objects
 -- real-only (no imaginary time component ) FFT
 -- a multi-dimensional FFT
 -- a command-line utility to perform ffts
 -- a command-line utility to perform fast-convolution filtering

 Then see kfc.h kfftr.h kfftnd.h fftutil.c kiss_fastfir.c
  in the tools/ directory.
*/

#ifdef USE_SIMD
    #include <xmmintrin.h>
    #define kfft_scalar __m128
    #define KFFT_MALLOC(nbytes) _mm_malloc(nbytes,16)
    #define KFFT_FREE _mm_free
#else
    #define KFFT_MALLOC malloc
    #define KFFT_FREE free
#endif


#ifndef kfft_scalar
    #define kfft_scalar double
#endif

#ifndef KFFT_RADER_LEVEL
    #define KFFT_RADER_LEVEL 3
#endif
#ifndef KFFR_RADER_LIMIT
    #define KFFT_RADER_LIMIT 50
#endif

typedef struct {
    kfft_scalar r;
    kfft_scalar i;
}kfft_cpx;

typedef struct kfftr_state * kfft_cfg;

kfft_cfg
kfft_config  (int         nfft,
                  int         inverse_fft,
                  void *      mem,
                  size_t *    lenmem);
/*
 nfft must be even

 If you don't care to allocate space, use mem = lenmem = NULL
*/


void
kfft (kfft_cfg              cfg,
          const kfft_scalar *   timedata,
          kfft_cpx *            freqdata);
/*
 input timedata has nfft scalar points
 output freqdata has nfft/2+1 complex points
*/

void
kffti (kfft_cfg             cfg,
           const kfft_cpx *     freqdata,
           kfft_scalar *        timedata);
/*
 input freqdata has  nfft/2+1 complex points
 output timedata has nfft scalar points
*/

int
kfft_next_fast_size(int n);

void
kfft_free (kfft_cfg* cfg);

static inline size_t
kfft_get_size (const int n) {
    size_t memneeded = 0;
    kfft_config (n, 0, NULL, &memneeded);
    return memneeded;
}

#ifdef __cplusplus
}
#endif

#endif
