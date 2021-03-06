valide_bench = files('valid_kfft_bench.py')

test_link = [libkfft]
c_args = []

if get_option('enable_half_scalar')
    if ccid == 'msvc'
        c_args += '/DTEST_HALF_SCALAR'
    else
        c_args += '-DTEST_HALF_SCALAR'
    endif
    
    fftw_dep = dependency ('fftw3f', required: false)
    if not fftw_dep.found()
        message ('FFTW3 not found. Try with compilator.')
        fftw_dep = cc.find_library('fftw3f', required: false)
    endif
else
    fftw_dep = dependency ('fftw3', required: false)
    if not fftw_dep.found()
        message ('FFTW3 not found. Try with compilator.')
        fftw_dep = cc.find_library('fftw3', required: false)
    endif
endif

if fftw_dep.found()
    fft_cpx_app = executable ('run_kf_cpx', files('kfft_fftw_cpx.c'),
                               c_args : c_args,
                               dependencies: fftw_dep,
                               link_with: libkfft,
                               include_directories: kfft_inc,
                               install: false)
    fft_scr_app = executable ('run_kf_scr', files('kfft_fftw_scr.c'),
                               c_args : c_args,
                               dependencies: fftw_dep,
                               link_with: libkfft,
                               include_directories: kfft_inc,
                               install: false)

    test_speed_cpx = run_target ('test_valide_cpx',
                              depends: [libkfft,fft_cpx_app],
                              command: [rt_python,
                                        valide_bench,
                                        fft_cpx_app,
                                        get_option('valide_test_limit').to_string(),
                                        get_option('valide_test_round').to_string()])
    test_speed_scr = run_target ('test_valide_scr',
                              depends: [libkfft,fft_scr_app],
                              command: [rt_python,
                                        valide_bench,
                                        fft_scr_app,
                                        get_option('valide_test_limit').to_string(),
                                        get_option('valide_test_round').to_string()])

else
    error('Validate test framework must have libfftw as dependency')
endif
