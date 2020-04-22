if not cc.has_header('immintrin.h')
    error ('Header <immintrin.h> not found. Disable enable_simd option. Maybe use fast-math optimization.')
endif
if (not (osid == 'windows')) and (not cc.has_header('cpuid.h'))
    error ('Header <cpuid.h> not found. Disable enable_simd option. Maybe use fast-math optimization.')
endif


if ccfam.startswith('x86')
    kfft_config += {'KFFT_ARCH_X86': true}
    
    # SSE2 check support
    sse2_example = cc.run('''
            #include <immintrin.h>
            int main () {
                double a[1],b[1],x;
                __m128d *av,*bv,xr;

                av = (__m128d*)a;
                bv = (__m128d*)b;
                xr = _mm_set1_pd(x);
                return 0;
            }''', name: 'SSE2 instructions')
    if sse2_example.compiled() and (sse2_example.returncode() == 0)
        kfft_config += {'KFFT_HAVE_SSE': true}
    endif

    # AVX check support
    avx_example = cc.run('''
            #include <immintrin.h>
            int main () {
                double a[1],b[1],x;
                __m256d *av,*bv,xr;

                av = (__m256d*)a;
                bv = (__m256d*)b;
                xr = _mm256_set1_pd(x);
                return 0;
            }''', args: '-mavx', name: 'AVX instructions')
    if avx_example.compiled() and (avx_example.returncode() == 0)
        kfft_config += {'KFFT_HAVE_AVX': true}
    endif

    # AVX2 check support
    avx2_example = cc.run('''
            #include <immintrin.h>
            #include <inttypes.h>

            int main () {
                unsigned a = 1, b = 1;
                __m256i av,bv,xr;

                av = _mm256_set1_epi32(a);
                bv = _mm256_set1_epi32(b);
                xr = _mm256_xor_si256(av, bv);
                return 0;
            }''', args: '-mavx2', name: 'AVX2 instructions')

    if avx2_example.compiled() and (avx2_example.returncode() == 0)
        kfft_config += {'KFFT_HAVE_AVX2': true}
    endif
endif # X86 architecture

if (ccfam == 'aarch64') or (ccid == 'arm')
    kfft_config += {'KFFT_ARCH_ARM': true}
    #TODO arm NEON extensions
endif

kfft_src += 'simd/kfft_simd.c'
kfft_headers += 'simd/kfft_simd.h'
kfft_incdirs += './simd'
