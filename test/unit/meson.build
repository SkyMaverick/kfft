test_name = 'test_bot'
test_src = files('''
    run.c
'''.split())
test_utils_src = files('''
    util.c
'''.split())

test_deps = [
    dependency('cunit', required: true),
    cc.find_library('dl', required: true),
]

test_inc = include_directories('.')
test_app = executable(test_name,
                      test_src + test_utils_src,
                      dependencies: test_deps,
                      include_directories: test_inc,
                      install: false)
subdir('suites')

testbot = run_target ('test_bot',
                      depends: test_app,
                      command: test_app)
test ('test_bot framewrok', test_app)
