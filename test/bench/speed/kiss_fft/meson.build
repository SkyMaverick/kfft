kiss_src = files ('''
    kiss_fft.c
    kiss_fftr.c
'''.split())

c_args=[]
if ccid == 'msvc'
    c_args = '''
        /DFIXED_POINT=32
    '''.split()
else
    c_args = '''
        -mtune=native
        -ffast-math
        -fomit-frame-pointer
        -unroll-loops
        -dA
        -fverbose-asm
        -DFIXED_POINT=32
    '''.split()
endif
kiss_dep = dependency ('openmp', language: 'c')

kiss_inc = include_directories('.')
kiss_dep = cc.find_library('m', required: false)

kiss_lib = shared_library('kiss', kiss_src,
                           c_args: c_args,
                           dependencies: kiss_dep,
                           install: false)
kiss_lib_static = static_library('kiss', kiss_src,
                           c_args: c_args,
                           dependencies: kiss_dep,
                           install: false)
