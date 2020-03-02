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
    kfft.h
    kfft_alloc.c
    kfft_cpx.c
    kfft_conv.c
    kfft_ext.c
    kfft_real.c
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

kfft_deps = []

# Validate system depends

if get_option ('enable_sysmath')
    if not (osid == 'windows')
        kfft_deps += cc.find_library('m', required: true)
    endif
endif
if get_option ('enable_simd')
    if not cc.has_header('xmmintrin.h')
        error ('Header <xmmintrin.h> not found. Disable enable_simd option. Maybe use fast-math optimization.')
    endif
endif
if get_option ('enable_fastmath')
    if not (ccid == 'msvc')
        message ('Use CC fast-math optimization. This maybe dangerouse optimization. Disable this if application crash.')
        kfft_cflags += '''
            -ffast-math
            -funroll-loops
        '''.split()
    endif
endif
if get_option ('enable_alloca')
    if not cc.has_header('alloca.h')
        error ('Header <alloca.h> not found. Disable enable_alloca option')
    endif
endif
if get_option ('build_with_native')
    if (ccid == 'gcc') or (ccid == 'clang')
        kfft_cflags += '''
            -march=native
        '''.split()
    endif
endif
build_type = get_option('build_library_api')
if build_type == 'export'
    kfft_defs += ['KFFT_EXPORTS']
elif build_type == 'import'
    kfft_defs += ['KFFT_IMPORTS']
else
    message ('No defines API')
endif

# Define configuration file

kfft_cfg_data = configuration_data({
    'KFFT_MAJOR_VERSION': kfft_ver[0],
    'KFFT_MINOR_VERSION': kfft_ver[1],
    'KFFT_PATCH_VERSION': kfft_ver[2],
    
    'KFFT_TRACE'        : get_option('enable_traces')   ,
    'KFFT_USE_SIMD'     : get_option('enable_simd')     ,
    'KFFT_USE_ALLOCA'   : get_option('enable_alloca')   ,
    'KFFT_USE_SYSMATH'  : get_option('enable_sysmath')  ,
    'KFFT_RADER_ALGO'   : get_option('enable_extended') ,
    'KFFT_MEMLESS_MODE' : get_option('enable_memless')  ,
    'KFFT_RADER_LIMIT'  : get_option('limit_rader')     ,
    'KFFT_PLAN_LEVEL'   : get_option('level_plan_max')  ,
    'KFFT_BFLY_LEVEL'   : get_option('level_bfly_max')  ,
})

kfft_cfg = configure_file (input: 'config.h.in',
                           output: 'kfft_config.h',
                           configuration: kfft_cfg_data)

# Apply all argments

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

kfft_inc= include_directories('''
        .
        ./incs
    '''.split())

optimize_opts = []
if not (get_option('buildtype') == 'debug')
    optimize_opts += [
            get_option('build_with_lto') ? 'b_lto=true' : 'b_lto=false',
            get_option('build_with_pgo') ? 'b_pgo=use'  : 'b_pg=none',
        ]
endif

# build library
libkfft = shared_library (meson.project_name(), kfft_src,
                          dependencies: kfft_deps,
                          gnu_symbol_visibility: 'hidden',
                          override_options : optimize_opts,
                          include_directories: kfft_inc,
                          install: true)

install_headers('kfft.h', subdir: 'kfft')
install_headers('''
        incs/kfft_alloc.h
        incs/kfft_cpx.h
        incs/kfft_conv.h
        incs/kfft_macro.h
        incs/kfft_real.h
        incs/kfft_system.h
    '''.split(), subdir: 'kfft/incs')
install_headers(kfft_cfg, subdir: 'kfft/incs')

if get_option('enable_tests')
    subdir('test')
endif
