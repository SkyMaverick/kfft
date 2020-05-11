# SSE2 check support
if cc.compiles('''
        #include <immintrin.h>
        int main () {
            double x; __m128d xr;
            xr = _mm_set1_pd(x);
            return 0;
        }''', name: 'SSE2 instructions')
    kfft_config += {'KFFT_HAVE_SSE2': true}

    sse_args = (ccid == 'msvc') ? '''
                    /ARCH:SSE2
                '''.split() : '''
                    -msse
                    -msse2
                '''.split()

    # SSE3 check support
    if (get_option('enable_sse3'))
        if cc.compiles('''
                #include <immintrin.h>
                int main () {
                    double x = 1; __m128d xr;
                    xr = _mm_loaddup_pd(&x);
                    return 0;
                }''', args: '-msse3', name: 'SSE3 instructions')
            kfft_config += {'KFFT_HAVE_SSE3': true}
        endif

    sse_args += (ccid == 'msvc') ? '''
                    /ARCH:SSE3
                '''.split() : '''
                    -msse3
                '''.split()
    endif

    kfft_vex += {
        'name' : 'sse',
        'src' : files ('''
                  kfft_bfly.c
                  kfft_generic.c
                  kfft_math.c
                '''.split()),
        'cargs': sse_args,
    }
endif
