# --- simple test suites ---
suites = {
    'test_sqrt': files ('''
                math_sqrt.c
            '''.split())
}

foreach targets, src : suites
    lib_suite = shared_library (targets, src + test_utils_src,
                                include_directories: [kfft_inc, test_inc],
                                link_with: libkfft,
                                dependencies: test_deps,
                                install: false)
endforeach
