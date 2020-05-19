# --- simple test suites ---
suites = {}

if (not get_option('enable_sysmath'))
    suites = {
        'test_sqrt': files ('math_sqrt.c'),
        'test_sincos': files ('math_sincos.c')
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
