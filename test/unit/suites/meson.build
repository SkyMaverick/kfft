# --- simple test suites ---
suites = {}

if (not get_option('enable_sysmath'))
    suites = {
        'api_complex_1d': files ('api_cpx1.c'),
        'api_complex_2d': files ('api_cpx2.c'),
        'api_scalar_1d': files ('api_sclr1.c'),
        'api_scalar_2d': files ('api_sclr2.c'),
        'test_sqrt': files ('math_sqrt.c'),
        'test_sincos': files ('math_sincos.c'),
        'sys_pool': files ('sys_pool.c')
    }
    test_deps += cc.find_library('m', required: true)
endif

foreach targets, src : suites
    lib_suite = shared_library (targets, src + test_utils_src,
                                include_directories: [kfft_inc, test_inc],
                                link_with: libkfft_static,
                                dependencies: test_deps,
                                install: false)
endforeach
