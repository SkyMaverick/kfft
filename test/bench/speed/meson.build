speed_bench = files('speed_bench.py')
speed_kfft_bench = files('speed_kfft_bench.py')

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
    fftw_cpx_app = executable ('run_fftw_cpx', files('tb_cpx_fftw.c'),
                               c_args : c_args,
                               dependencies: fftw_dep,
                               include_directories: kfft_inc,
                               install: false)
    fftw_scr_app = executable ('run_fftw_scr', files('tb_scr_fftw.c'),
                               c_args : c_args,
                               dependencies: fftw_dep,
                               include_directories: kfft_inc,
                               install: false)
endif

kfft_cpx_app = executable ('run_kfft_cpx', files('tb_cpx_kfft.c'),
                           c_args : c_args,
                           link_with: test_link,
                           include_directories: kfft_inc,
                           install: false)
kfft_scr_app = executable ('run_kfft_scr', files('tb_scr_kfft.c'),
                           c_args : c_args,
                           link_with: test_link,
                           include_directories: kfft_inc,
                           install: false)
subdir('kiss_fft')
kiss_cpx_app = executable ('run_kiss_cpx', files('tb_cpx_kiss.c'),
                           c_args : c_args,
                           link_with: kiss_lib,
                           include_directories: kiss_inc,
                           install: false)
kiss_scr_app = executable ('run_kiss_scr', files('tb_scr_kiss.c'),
                           c_args : c_args,
                           link_with: kiss_lib,
                           include_directories: kiss_inc,
                           install: false)


test_speed_cpx = run_target ('test_speed_cpx',
                          command: [rt_python,
                                    speed_bench,
                                    kfft_cpx_app,
                                    fftw_cpx_app,
                                    kiss_cpx_app,
                                    get_option('speed_test_sequence_limit').to_string(),
                                    get_option('speed_test_step').to_string(),
                                    meson.current_build_dir()])
test_speed_kfft_cpx = run_target ('test_speed_kfft_cpx',
                              command: [rt_python,
                                        speed_kfft_bench,
                                        kfft_cpx_app,
                                        get_option('speed_test_sequence_limit').to_string(),
                                        get_option('speed_test_step').to_string(),
                                        meson.current_build_dir()])
test_speed_scr = run_target ('test_speed_scr',
                          command: [rt_python,
                                    speed_bench,
                                    kfft_scr_app,
                                    fftw_scr_app,
                                    kiss_scr_app,
                                    get_option('speed_test_sequence_limit').to_string(),
                                    get_option('speed_test_step').to_string(),
                                    meson.current_build_dir()])
speed_kfft_test_scr = run_target ('test_speed_kfft_scr',
                              command: [rt_python,
                                        speed_kfft_bench,
                                        kfft_scr_app,
                                        get_option('speed_test_sequence_limit').to_string(),
                                        get_option('speed_test_step').to_string(),
                                        meson.current_build_dir()])
