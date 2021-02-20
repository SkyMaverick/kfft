if get_option('enable_half_scalar')
    subdir('half')
else
    subdir('norm')
endif

kfft_inc += include_directories('''
    .
'''.split())

kfft_vex += {
    'name' : 'sse',
    'src' : files ('''
              kfft_sse.c
            '''.split()),
    'cargs': sse_args,
}
