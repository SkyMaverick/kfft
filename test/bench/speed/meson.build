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
