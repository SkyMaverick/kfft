fftw_dep = dependency ('fftw3', required: false)
if not fftw_dep.found()
    message ('FFTW3 not found. Try with compilator.')
    fftw_dep = cc.found_library('fftw3', required: false)
endif

test_link = [libkfft]

app_demo =      executable ('demo', 'demo.c',
                            link_with: test_link,
                            include_directories: kfft_inc,
                            install: false)

test_src = '''
    test.c
'''.split()

app_usdemo =    executable ('usdemo', test_src,
                            c_args: '-DSIMPLE_APP',
                            link_with: test_link,
                            include_directories: kfft_inc,
                            install: false)

app_ustest =    executable ('ustest', test_src,
                            link_with: test_link,
                            include_directories: kfft_inc,
                            install: false)

if fftw_dep.found()
    app_ukfcmp =    executable ('ukfcmp', test_src,
                            c_args: ['-DFFTW_COMPARE'],
                            link_with: test_link,
                            dependencies: fftw_dep,
                            include_directories: kfft_inc,
                            install: false)
    app_ukfcmpwp =  executable ('ukfcmpwp', test_src,
                            c_args: ['-DFFTW_COMPARE', '-DCHECK_WITHOUT_PLAN'],
                            link_with: test_link,
                            dependencies: fftw_dep,
                            include_directories: kfft_inc,
                            install: false)
    app_aukfcmp =   executable ('aukfcmp', test_src,
                            c_args: ['-DFFTW_COMPARE', '-DLOG_AUTO'],
                            link_with: test_link,
                            dependencies: fftw_dep,
                            include_directories: kfft_inc,
                            install: false)
    app_aukfcmpwp = executable ('aukfcmpwp', test_src,
                            c_args: ['-DFFTW_COMPARE', '-DCHEC_WITHOUT_PLAN', '-DLOG_AUTO'],
                            link_with: test_link,
                            dependencies: fftw_dep,
                            include_directories: kfft_inc,
                            install: false)
endif
