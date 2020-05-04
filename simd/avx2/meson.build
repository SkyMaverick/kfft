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

    avx2_args = (not (ccid == 'msvc')) ? '''
                    -mavx2
                '''.split() : []

    kfft_vex += {
        'name' : 'avx2',
        'src' : files ('''
                '''.split()),
        'cargs': avx2_args,
    }
endif
