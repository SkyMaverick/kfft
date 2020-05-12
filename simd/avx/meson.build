# AVX check support
avx_args = (ccid == 'msvc') ? '''
                /arch:AVX
            '''.split() : '''
                -mavx
            '''.split()
if cc.compiles('''
        #include <immintrin.h>
        int main () {
            double x; __m256d xr;
            xr = _mm256_set1_pd(x);
            return 0;
        }''', args: avx_args, name: 'AVX instructions')
    kfft_config += {'KFFT_HAVE_AVX': true}

    kfft_vex += {
        'name' : 'avx',
        'src' : files ('''
                  kfft_bfly.c
                  kfft_generic.c
                  kfft_math.c
                '''.split()),
        'cargs': avx_args,
    }
endif
