project ('kfft', 'c',
         version: '0.1.3',
         license: 'Zlib',
         default_options:[
            'c_std=gnu11',
            'warning_level=3',
         ],
         meson_version: '>=0.52')

message ('Build configuration - @0@'.format(get_option('buildtype')))

kfft_src = '''
    kfft.c
    kfft.h
    kfft_alloc.c
    kfft_alloc.h
    kfft_core.c
    kfft_guts.h
    kfft_rader.c
'''.split()

cc = meson.get_compiler('c')
ccid = cc.get_id()
osid = target_machine.system()

kfft_cflags  = '''
'''.split()

kfft_ldflags = '''
'''.split()

kfft_defs    = '''
'''.split()

if get_option('enable_traces')
    kfft_defs += 'KFFT_TRACE'
endif
if get_option('enable_simd')
    if (cc.has_header('xmmintrin.h'))
        kfft_defs += 'KFFT_USE_SIMD'
    else
        warning ('SIMD header not found. SIMD disabled')
    endif
endif
if get_option('enable_alloca')
    if (cc.has_header('alloca.h'))
        kfft_defs += 'KFFT_USE_ALLOCA'
    else
        warning ('Alloca not found on this system. Disabled.')
    endif
endif
if get_option('enable_extended')
    kfft_defs += 'KFFT_RADER_ALGO'
endif

foreach item : kfft_cflags
    if cc.has_argument (item)
        add_project_arguments (item, language: 'c')
    endif
endforeach
foreach item : kfft_ldflags
    if cc.has_link_argument (item)
        add_project_link_arguments (item, language: 'c')
    endif
endforeach
foreach item : kfft_defs
    item = ((ccid == 'msvc') ? '/D': '-D') + item
    add_project_arguments (item, language: 'c')
endforeach

kfft_deps = []
if not (osid == 'windows')
    kfft_deps += cc.find_library('m', required: true)
endif

kfft_inc= include_directories('.')

libkfft = shared_library (meson.project_name(), kfft_src,
                          dependencies: kfft_deps,
                          install: true)

if get_option('enable_tests')
    subdir('test')
endif
