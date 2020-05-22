tests = '''
    speed
'''.split()

rt_python = import('python').find_installation()
if rt_python.found()
    rt_python_ver = run_command (rt_python, '--version').stdout().strip()
    message ('Python version found: @0@'.format(rt_python_ver))
    
    foreach test : tests
        subdir(test)
    endforeach
else
    
endif
