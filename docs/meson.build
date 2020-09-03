app_doxy = find_program('doxygen', required: false)

if app_doxy.found()
    doxy_file = files('''
        Doxyfile
    '''.split())

    doc_target = run_target ( 'docs',
                              command:[app_doxy, doxy_file])
else
    warning ('''Documentation files not create.
              Please install doxygen and run this target again.''')
endif
