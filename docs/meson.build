doxy_app = find_program('doxygen')

if doxy_app.found()
    doxy_file = files('''
        Doxyfile
    '''.split())

    doc_target = run_target ( 'docs', command:[doxy_app, doxy_file] )
else
    warning ('Documentation files not create. Please install doxygen and run this target again.')
endif
