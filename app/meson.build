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

app_fft = executable (app_name, 'fft_app.c',
                       link_with: test_link,
                       include_directories: kfft_inc,
                       install: true)
