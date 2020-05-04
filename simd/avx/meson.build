# AVX check support
if cc.compiles('''
        #include <immintrin.h>
        int main () {
            double x; __m256d xr;
            xr = _mm256_set1_pd(x);
            return 0;
        }''', args: '-mavx', name: 'AVX instructions')
    kfft_config += {'KFFT_HAVE_AVX': true}

    avx_args = (not (ccid == 'msvc')) ? '''
                    -mavx
                '''.split() : []

    kfft_vex += {
        'name' : 'avx',
        'src' : files ('''
                '''.split()),
        'cargs': avx_args,
    }
endif
