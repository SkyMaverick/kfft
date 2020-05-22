speed_bench = files('speed_bench.py')

fftw_dep = dependency ('fftw3', required: false)
if not fftw_dep.found()
    message ('FFTW3 not found. Try with compilator.')
    fftw_dep = cc.find_library('fftw3', required: false)
endif

test_link = [libkfft]

kfft_cpx_app = executable ('run_kfft', files('tb_cpx_kfft.c'),
                           link_with: test_link,
                           include_directories: kfft_inc,
                           install: false)
fftw_cpx_app = executable ('run_fftw', files('tb_cpx_fftw.c'),
                           dependencies: fftw_dep,
                           include_directories: kfft_inc,
                           install: false)
subdir('kiss_fft')
kiss_cpx_app = executable ('run_kiss', files('tb_cpx_kiss.c'),
                           link_with: kiss_lib,
                           include_directories: kiss_inc,
                           install: false)


speed_test = run_target ('speed_test',
                          command: [rt_python,
                                    speed_bench,
                                    kfft_cpx_app,
                                    fftw_cpx_app,
                                    kiss_cpx_app,
                                    get_option('speed_test_sequence_limit').to_string(),
                                    get_option('speed_test_step').to_string(),
                                    meson.current_build_dir()])
