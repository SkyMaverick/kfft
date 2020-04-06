fftw_dep = dependency ('fftw3', required: false)
if not fftw_dep.found()
    message ('FFTW3 not found. Try with compilator.')
    fftw_dep = cc.find_library('fftw3', required: false)
endif

test_link = [libkfft]

test_src = '''
    app.c
'''.split()

kfft_apps = {
    'usdemo': ['-DSIMPLE_APP'],
    'ustest': []
}

foreach test_app, app_args : kfft_apps
    A = executable (test_app, test_src,
                    c_args: app_args,
                    link_with: test_link,
                    include_directories: kfft_inc,
                    install: false)
endforeach

if fftw_dep.found()
    fftw_apps = {
        'ukfcmp' : ['-DFFTW_COMPARE'],
        'ukfcmpwp' : ['-DFFTW_COMPARE', '-DCHECK_WITHOUT_PLAN'],
        'aukfcmp' : ['-DFFTW_COMPARE', '-DLOG_AUTO'],
        'aukfcmpwp' : ['-DFFTW_COMPARE', '-DCHECK_WITHOUT_PLAN', '-DLOG_AUTO']
    }

    foreach test_app, app_args : fftw_apps
        A = executable (test_app, test_src,
                        c_args: app_args,
                        link_with: test_link,
                        dependencies: fftw_dep,
                        include_directories: kfft_inc,
                        install: false)
    endforeach
endif
