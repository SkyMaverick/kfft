test_name = 'test_bot'
test_src = files('''
    run.c
    util.c
    util.h
'''.split())

test_deps = [
    dependency('cunit', required: true),
    cc.find_library('dl', required: true),
]

test_app = executable(test_name,
                      test_src,
                      dependencies: test_deps,
                      install: false)
subdir('suites')

testbot = run_target ('test_bot',
                      depends: test_app,
                      command: test_app)
test ('test_bot framewrok', test_app)
