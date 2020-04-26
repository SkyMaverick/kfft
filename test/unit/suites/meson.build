# --- simple test suites ---
suites = {
    #TODO 'suite' : 'main_app depends files'
}

foreach item_file, item_adds : suites
    lib_name = item_file.split('.')[0]
    lib_suite = shared_library (lib_name, item_file + item_adds,
                                include_directories: kfft_inc,
                                link_with: libkfft,
                                install: false)
endforeach
