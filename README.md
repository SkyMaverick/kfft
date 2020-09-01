## ABOUT
LibKFFT is a free/open-source (ZLib licensed) C library providing fast fourier transform and other related operations for comlex and scalar double or float sequences.
## DEVELOPMENT
LibKFFT is deeply redesigned fork  [KissFFT project](https://github.com/mborgerding/kissfft) by Mark Borgerding.

## REALISED FEATURES
The features of LibKFFT include:
*   complex and scalar fourier transformation operations (forward and inverse)
    *   64-bit float type (double) as default
    *   32-bit float type (float) if set option KFFT_HALF_SCALAR
*   SIMD acceleration transform operations (now support SSE, SSE2, SSE3 (optional)) (build option)
*   Recursive  [Rader algorithm](https://en.wikipedia.org/wiki/Rader's_FFT_algorithm) for prime-number lenght sequenses (build option)
*   Internal math functions (aka sqrt, sincos etc.) without system libm functions (build option)
*   Parallelization with OpenMP framework (build option)
*   2D, sparse, convolution(+2D) operations as build-time extensions
*   Detail trace log in stdout for diagnostic (build option)
*   Tiny API support for simplify code in trivial usage

## BUILD AND TEST
__Need sofware__: meson >= 55.0, ninja >= 1.8.0
Test build on __OS with CC__:
*   _Linux_ (x86-64): GCC 9.3, Clang-10
*   _Windows_ (x86-64): MSVC 2015, 2017, Clang-10
*   _RPi3_ (armhf): GCC 8.3

__Build__: configure build in ___meson_options.txt___ if need and run `meson . ./<builddir> && ninja -C ./<builddir>`

__Test__ _(don't implemented on Windows)_: if change option `enable_tests` use test `ninja -C ./<builddir> <target>`
*  *test_bot* - run unit tests (need libcunit for build it)
*  *test_speed[_cpx,_scr]* - run speed tests and compare with popular implementations
*  *test_speed_kfft[_cpx,_scr]* - run speed tests for libkfft only (produce SVG and CSV)
*  *test_valide[_cpx,_scr]* - run validate tests (use FFTW3 for compare)

## BENCHMARKS
Compare with popular implementations
Test machine:
*   __CPU__ - Intel i5-2540M (Sandy)
*   __RAM__ - Kingstone DDR3 1333MHz 16Gb
*   __OS__  - Linux Mint 20

Implementation:
*   FFTW (3.3.8)
*   KissFFT (commit b2e0e60).
*   KFFT (0.8.0)

Build options:
*   __libfftw__: *build in Ubuntu 20.04 repo*
*   __libkiss__: `-DFIXED_POINT=32`
*   __libkfft__: `-DKFFT_USE_SIMD -DKFFT_USE_SYSMATH -DKFFT_USE_OPENMP -DKFFT_RADER_ALGO -DKFFT_RADER_LIMIT=25 -DKFFT_PLAN_LEVEL=7 -DKFFT_BFLY_LEVEL=5`
... libkfft and libkiss make with [Meson](http://https://mesonbuild.com) mode `Release build with PGO and LTO

Bench variants:
*   __Range__: 1 .. 1000000
*   __Step__: 101

This step may test all variants such as: normal lenght, prime-lenght and other difficult cases

## Bench result
__Complex test evaluation compare__
![complex test](https://raw.githubusercontent.com/SkyMaverick/kfft/master/docs/img/cmpcpx.svg "Complex 1D sequence transform")
__Scalar test evaluation compare__
![scalar test](https://raw.githubusercontent.com/SkyMaverick/kfft/master/docs/img/cmpscr.svg "Scalar 1D sequence transform")
