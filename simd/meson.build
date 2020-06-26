if not cc.has_header('immintrin.h')
    error ('Header <immintrin.h> not found. Disable enable_simd option. Maybe use fast-math optimization.')
endif
if (not (osid == 'windows')) and (not cc.has_header('cpuid.h'))
    error ('Header <cpuid.h> not found. Disable enable_simd option. Maybe use fast-math optimization.')
endif

kfft_config += {'KFFT_USE_SIMD': true}

if ccfam.startswith('x86')
    kfft_config += {'KFFT_ARCH_X86': true}
    
    if not get_option('enable_half_scalar')
        if get_option('list_vex86').contains('avx')
            subdir('avx')
        endif
    endif
    if get_option('list_vex86').contains('sse')
        subdir('sse')
    endif
endif # X86 architecture

if (ccfam == 'aarch64') or (ccid == 'arm')
    kfft_config += {'KFFT_ARCH_ARM': true}
    #TODO arm NEON extensions
endif

kfft_src += files('''
    kfft_simd_check.c
'''.split())

kfft_headers += files('''
    kfft_simd.h
    kfft_simd_check.h
'''.split())

kfft_inc += include_directories('''
    .
    sse
    avx
'''.split())
