static int
prepare_2d(app_mode_t* M) {
    if ((M->len > 0) && (M->x > 0) && (M->len % M->x > 0))
        return 1;

    M->y = M->len / M->x;
    return 0;
}

static inline char*
cmd_line_parse(int argc, char* argv[], app_mode_t* mode) {
    char* ret = NULL;

    int opt = 0;
    while ((opt = getopt(argc, argv, FMT_OPTSTRING)) != -1) {
        switch (opt) {
        case 'g':
            mode->flags |= KFFT_FLAG_GENERIC;
            break;
        case 'G':
            mode->flags |= KFFT_FLAG_GENERIC_ONLY;
            break;
        case 'i':
            mode->flags |= KFFT_FLAG_INVERSE;
            break;
        case 's':
            mode->is_shift = true;
            break;
        case 'S':
            mode->is_cpx = false;
            break;
        case 'x':
            mode->x = atol(optarg);
            if (mode->x == 0)
                exit(1);
            mode->is_2d = true;
            break;
        case '?':
            display_help();
            exit(0);
        case 'v':
            fprintf(stdout, "%d.%d.%d\n", VER_MAJOR, VER_MINOR, VER_PATCH);
            exit(0);
        case 'V':
            display_info();
            exit(0);
        }
    }

    if (optind < argc) {
        size_t tmp_size = 1;
        char* buf = calloc(1, tmp_size);
        if (buf) {
            do {
                char* old = buf;
                tmp_size += strlen(argv[optind]) + 1;

                buf = realloc(buf, tmp_size);
                if (buf == NULL) {
                    free(old);
                    goto bailout;
                }

                strcat(buf, argv[optind]);
                strcat(buf, " ");
            } while (argc > ++optind);

            buf[strlen(buf) - 1] = '\0';
        } /* buf allocated */
        ret = buf;
    }

    if (ret == NULL) {
        ret = read_stdin();
        mode->is_stdin = true;
    }
bailout:
    return ret;
}
