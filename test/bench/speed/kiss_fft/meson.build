kiss_src = files ('''
    kiss_fft.c
    kiss_fftr.c
'''.split())

c_args = '''
    -mtune=native
    -ffast-math
    -fomit-frame-pointer
    -unroll-loops
    -dA
    -fverbose-asm
'''.split()

kiss_inc = include_directories('.')
kiss_dep = cc.find_library('m', required: false)

kiss_lib = shared_library('kiss', kiss_src,
                           c_args: c_args,
                           dependencies: kiss_dep,
                           install: false)
