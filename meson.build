project ('kfft', 'c',
         version: '0.1.5',
         license: 'Zlib',
         default_options:[
            'c_std=gnu11',
            'warning_level=3',
         ],
         meson_version: '>=0.52')

message ('Build configuration - @0@'.format(get_option('buildtype')))

kfft_ver = meson.project_version().split('.')

kfft_src = '''
    kfft_guts.h
    kfft.h
    kfft.c
    kfft_alloc.c
    kfft_alloc.h
    kfft_core.c
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

kfft_cfg = configuration_data()
kfft_cfg.set('KFFT_MAJOR_VERSION', kfft_ver[0])
kfft_cfg.set('KFFT_MINOR_VERSION', kfft_ver[1])
kfft_cfg.set('KFFT_PATCH_VERSION', kfft_ver[2])

kfft_cfg.set('KFFT_TRACE', get_option('enable_traces'))
kfft_cfg.set('KFFT_USE_SIMD', get_option('enable_simd'))
kfft_cfg.set('KFFT_USE_ALLOCA', get_option('enable_alloca'))
kfft_cfg.set('KFFT_USE_SYSMATH', get_option('enable_sysmath'))
kfft_cfg.set('KFFT_RADER_ALGO', get_option('enable_extended'))
kfft_cfg.set('KFFT_MEMLESS_MODE', get_option('enable_memless'))

configure_file (input: 'config.h.in',
                output: 'config.h',
                configuration: kfft_cfg)

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
