if not cc.has_header('immintrin.h')
    error ('Header <immintrin.h> not found. Disable enable_simd option. Maybe use fast-math optimization.')
endif
if (not (osid == 'windows')) and (not cc.has_header('cpuid.h'))
    error ('Header <cpuid.h> not found. Disable enable_simd option. Maybe use fast-math optimization.')
endif


if ccfam.startswith('x86')
    kfft_config += {'KFFT_ARCH_X86': true}
    
# SSE check support
    if cc.compiles('''
            #include <immintrin.h>
            int main () {
                float x; __m128 xr;
                xr = _mm_set1_ps(x);
                return 0;
            }''', name: 'SSE instructions')

        kfft_config += {'KFFT_HAVE_SSE': true}
    endif

# SSE2 check support
    if cc.compiles('''
            #include <immintrin.h>
            int main () {
                double x; __m128d xr;
                xr = _mm_set1_pd(x);
                return 0;
            }''', name: 'SSE2 instructions')
        kfft_config += {'KFFT_HAVE_SSE2': true}
    endif

# SSE3 check support
    if cc.compiles('''
            #include <immintrin.h>
            int main () {
                double x = 1; __m128d xr;
                xr = _mm_loaddup_pd(&x);
                return 0;
            }''', args: '-msse3', name: 'SSE3 instructions')
        kfft_config += {'KFFT_HAVE_SSE3': true}
    endif

# AVX check support
    if cc.compiles('''
            #include <immintrin.h>
            int main () {
                double x; __m256d xr;
                xr = _mm256_set1_pd(x);
                return 0;
            }''', args: '-mavx', name: 'AVX instructions')
        kfft_config += {'KFFT_HAVE_AVX': true}
    endif

# AVX2 check support
    if cc.compiles('''
            #include <immintrin.h>
            int main () {
                unsigned a = 1,b = 1;
                __m256i av,bv,xr;
                av = _mm256_set1_epi32(a);
                bv = _mm256_set1_epi32(b);
                xr = _mm256_xor_si256(av, bv);
                return 0;
            }''', args: '-mavx2', name: 'AVX2 instructions')
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
