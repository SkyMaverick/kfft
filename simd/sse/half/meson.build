# SSE2 check support
sse_args = (ccid == 'msvc') ? '''
                /arch:SSE
            '''.split() : '''
                -msse
            '''.split()

if cc.compiles('''
        #include <immintrin.h>
        int main () {
            float x; __m128 xr;
            xr = _mm_set1_ps(x);
            return 0;
        }''', args: sse_args, name: 'SSE instructions')
    kfft_config += {'KFFT_HAVE_SSE': true}
    
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
endif
