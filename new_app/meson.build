app_name = 'kfft'

test_link = [libkfft]

configure_file (
    output: 'config.h',
    configuration : {
        'VER_MAJOR': meson.project_version().split('.')[0],
        'VER_MINOR': meson.project_version().split('.')[1],
        'VER_PATCH': meson.project_version().split('.')[2],
        
        'APP_NAME': '"' + app_name + '"',
    
    }
)

app_dep = []
if not (osid == 'windows')
    app_dep += cc.find_library('dl', required: true)
endif
app_rpath = ':'.join ([
    '$ORIGIN'/'.',
    '$ORIGIN'/'..'
])

app_fft = executable (app_name, 'kfft_main.c',
                       link_with: test_link,
                       include_directories: kfft_inc,
                       dependencies: app_dep,
                       build_rpath: app_rpath,
                       install: true)
