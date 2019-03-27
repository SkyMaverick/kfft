#ifndef KISS_FFT_H
#define KISS_FFT_H

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

 Then see kfc.h kiss_fftr.h kiss_fftnd.h fftutil.c kiss_fastfir.c
  in the tools/ directory.
*/

#ifdef USE_SIMD
    #include <xmmintrin.h>
    #define kiss_fft_scalar __m128
    #define KISS_FFT_MALLOC(nbytes) _mm_malloc(nbytes,16)
    #define KISS_FFT_FREE _mm_free
#else
    #define KISS_FFT_MALLOC malloc
    #define KISS_FFT_FREE free
#endif


#ifndef kiss_fft_scalar
    #define kiss_fft_scalar double
#endif

typedef struct {
    kiss_fft_scalar r;
    kiss_fft_scalar i;
}kiss_fft_cpx;

typedef struct kiss_fftr_state * kiss_fft_cfg;

kiss_fft_cfg
kiss_fft_config  (int         nfft,
                  int         inverse_fft,
                  int         delta,
                  int         step,
                  void *      mem,
                  size_t *    lenmem);
/*
 nfft must be even

 If you don't care to allocate space, use mem = lenmem = NULL
*/


void
kiss_fft (kiss_fft_cfg              cfg,
          const kiss_fft_scalar *   timedata,
          kiss_fft_cpx *            freqdata);
/*
 input timedata has nfft scalar points
 output freqdata has nfft/2+1 complex points
*/

void
kiss_ffti (kiss_fft_cfg             cfg,
           const kiss_fft_cpx *     freqdata,
           kiss_fft_scalar *        timedata);
/*
 input freqdata has  nfft/2+1 complex points
 output timedata has nfft scalar points
*/

int
kiss_fft_next_fast_size(int n);

void
kiss_fft_free (kiss_fft_cfg* cfg);

static inline size_t
kiss_fft_get_size (const int n) {
    size_t memneeded = 0;
    kiss_fft_config (n, 0, 0, 0, NULL, &memneeded);
    return memneeded;
}

#ifdef __cplusplus
}
#endif

#endif
