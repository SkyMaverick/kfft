# SSE2 check support
sse_args = (ccid == 'msvc') ? '''
                /arch:SSE2
            '''.split() : '''
                -msse
                -msse2
            '''.split()

if cc.compiles('''
        #include <immintrin.h>
        int main () {
            double x; __m128d xr;
            xr = _mm_set1_pd(x);
            return 0;
        }''', args: sse_args, name: 'SSE2 instructions')
    kfft_config += {'KFFT_HAVE_SSE2': true}

    # SSE3 check support
    if (get_option('enable_sse3'))
        sse_args = (ccid == 'msvc') ? '''
                        /arch:AVX
                    '''.split() : '''
                        -msse3
                    '''.split()
        if cc.compiles('''
                #include <immintrin.h>
                int main () {
                    double x = 1; __m128d xr;
                    xr = _mm_loaddup_pd(&x);
                    return 0;
                }''', args: sse_args, name: 'SSE3 instructions')
            kfft_config += {'KFFT_HAVE_SSE3': true}
        endif
    endif
else
    error ('Target machine don\'t support SSE2 instructions')
endif
