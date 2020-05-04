if not cc.has_header('immintrin.h')
    error ('Header <immintrin.h> not found. Disable enable_simd option. Maybe use fast-math optimization.')
endif
if (not (osid == 'windows')) and (not cc.has_header('cpuid.h'))
    error ('Header <cpuid.h> not found. Disable enable_simd option. Maybe use fast-math optimization.')
endif


if ccfam.startswith('x86')
    kfft_config += {'KFFT_ARCH_X86': true}
    
    if get_option('list_vex86').contains('sse')
        subdir('sse')
    endif

    if get_option('list_vex86').contains('avx')
        subdir('avx')
    endif

    if get_option('list_vex86').contains('avx2')
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
    endif
endif # X86 architecture

if (ccfam == 'aarch64') or (ccid == 'arm')
    kfft_config += {'KFFT_ARCH_ARM': true}
    #TODO arm NEON extensions
endif

kfft_src += files('''
    kfft_simd_check.c
'''.split())

kfft_headers += files('''
    kfft_simd.h
    kfft_simd_check.h
'''.split())

kfft_incdirs += '''
    ./simd
    ./simd/sse
    ./simd/avx
    ./simd/avx2
'''.split()
