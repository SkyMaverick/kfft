project ('kfft', 'c',
         version: '0.7.0',
         license: 'Zlib',
         default_options:[
            'c_std=gnu11',
            'warning_level=3',
            'strip=true',
            'b_pie=true'
         ],
         meson_version: '>=0.52')

message ('Build configuration - @0@'.format(get_option('buildtype')))

kfft_ver = meson.project_version().split('.')

kfft_src = files('''
    kfft.h
    kfft_pool.c
    kfft_cpx.c
    kfft_scalar.c
    kfft_ext.c
    kfft_math.c
    kfft_shift.c
'''.split())

kfft_headers = files('''
    incs/kfft_pool.h
    incs/kfft_algo.h
    incs/kfft_cpx.h
    incs/kfft_ext.h
    incs/kfft_macro.h
    incs/kfft_math.h
    incs/kfft_scalar.h
    incs/kfft_shift.h
    incs/kfft_system.h
    incs/kfft_trace.h
    incs/kfft_types.h
'''.split())

kfft_inc=[]

kfft_inc += include_directories('''
    .
    incs
    ext
'''.split())

cc = meson.get_compiler('c')
ccid = cc.get_id()
ccfam = target_machine.cpu_family()
osid = target_machine.system()

kfft_cflags  = '''
'''.split()

kfft_ldflags = '''
'''.split()

kfft_defs    = '''
'''.split()

kfft_deps = []

kfft_config = {
    'KFFT_MAJOR_VERSION': kfft_ver[0],
    'KFFT_MINOR_VERSION': kfft_ver[1],
    'KFFT_PATCH_VERSION': kfft_ver[2],

    'KFFT_DYNAPI_ENABLE': get_option('enable_dynamical_api'),
    
    'KFFT_TRACE'        : get_option('enable_traces')       ,
    'KFFT_USE_SIMD'     : get_option('enable_simd') and
                          (not get_option('enable_half_scalar')),
    'KFFT_USE_ALLOCA'   : get_option('enable_alloca')       ,
    'KFFT_USE_SYSMATH'  : get_option('enable_sysmath')      ,
    'KFFT_RADER_ALGO'   : get_option('enable_extended')     ,
    'KFFT_MEMLESS_MODE' : get_option('enable_memless')      ,
    'KFFT_HALF_SCALAR'  : get_option('enable_half_scalar')  ,

    'KFFT_RADER_LIMIT'  : get_option('limit_rader')         ,
    'KFFT_PLAN_LEVEL'   : get_option('level_plan_max')      ,
    'KFFT_BFLY_LEVEL'   : get_option('level_bfly_max')      ,
}

if get_option ('enable_sysmath')
    if not (osid == 'windows')
        kfft_deps += cc.find_library('m', required: true)
    endif
else
    subdir('math')
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
if get_option ('enable_openmp')
    fomp = dependency ('openmp', language: 'c')
    if fomp.found()
        kfft_config += {'KFFT_USE_OPENMP': get_option('enable_openmp')}
        kfft_deps += fomp
    endif
endif
if get_option ('build_with_native')
    if (ccid == 'gcc') or (ccid == 'clang')
        kfft_cflags += '''
            -march=native
        '''.split()
    endif
else
    if (ccid == 'gcc') or (ccid == 'clang')
        if (ccfam == 'x86')
            kfft_cflags += '''
                -march=i686
                -mtune=generic
            '''.split()
        endif
        if (ccfam == 'x86_64')
            kfft_cflags += '''
                -march=x86-64
                -mtune=generic
            '''.split()
        endif
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

# Validate system depends
kfft_vex = []
if get_option('enable_simd')
    subdir('simd')
endif # Enable SIMD

# Add all checked extensions
subdir('ext')

# Define configuration file
kfft_cfg = configure_file (input: 'config.h.in',
                           output: 'kfft_config.h',
                           configuration: configuration_data(kfft_config))

optimize_opts = []
if not (get_option('buildtype') == 'debug')
    optimize_opts += [
            get_option('build_with_lto') ? 'b_lto=true' : 'b_lto=false',
            get_option('build_with_pgo') ? 'b_pgo=use'  : 'b_pg=none',
        ]
endif

#kfft_inc = include_directories (kfft_incdirs)

kfft_link = []
foreach vex : kfft_vex
    kfft_link += static_library (vex.get('name'), vex.get('src'),
                    include_directories: kfft_inc,
                    c_args: vex.get('cargs'),
                    gnu_symbol_visibility: 'hidden',
                    install: false)
endforeach
# build library
libkfft = shared_library (meson.project_name(), kfft_src,
                          dependencies: kfft_deps,
                          gnu_symbol_visibility: 'hidden',
                          override_options : optimize_opts,
                          include_directories: kfft_inc,
                          link_with: kfft_link,
                          install: true)

install_headers('kfft.h', subdir: 'kfft')
install_headers('kfft_tiny.h', subdir: 'kfft')
install_headers(kfft_headers, subdir: 'kfft/incs')
install_headers(kfft_cfg, subdir: 'kfft/incs')

if get_option('enable_tests')
    libkfft_static = static_library (meson.project_name(), kfft_src,
                              dependencies: kfft_deps,
                              override_options : optimize_opts,
                              include_directories: kfft_inc,
                              link_with: kfft_link,
                              install: false)

    subdir('test')
endif

if meson.is_subproject()
    kfft_lib = declare_dependency (
                link_with: libkfft,
                include_directories: kfft_inc
            )
endif
