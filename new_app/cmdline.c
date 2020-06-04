#include "info.c"

#if defined(KFFT_2D_ENABLE)
static int
prepare_2d(state_t* M) {
    //    if ((M->len > 0) && (M->x > 0) && (M->len % M->x > 0))
    //        return 1;
    //
    //    M->y = M->len / M->x;
    //    return 0;
}
#endif /* KFFT_2D_ENABLE */

static int
parse_sparse_arg(char* arg, state_t* M) {
    //    size_t olen = strlen(arg);
    //    char* buf = calloc(olen + 1, sizeof(char));
    //    if (buf) {
    //        strncpy(buf, optarg, olen);
    //
    //        char* s = strchr(buf, ':');
    //        if (s)
    //            *s = '\0';
    //
    //        M->dim = atoi(buf);
    //        if (strlen(buf) < olen)
    //            M->step = atoi(s + 1);
    //
    //        M->is_sparse = true;
    //        free(buf);
    //
    //        return 0;
    //    } else {
    //        return 1;
    //    }
}

static inline char*
pipe_read_stdin(state_t* st) {
    // TODO
    return NULL;
}

static inline char*
manual_read_stdin(int argc, char* argv[], int idx) {
    return NULL;
}

static char*
cmd_line_parse(int argc, char* argv[], state_t* st) {
    char* ret = NULL;

    int opt = 0;
    while ((opt = getopt(argc, argv, FMT_OPTSTRING)) != -1) {
        switch (opt) {
        case 'b':
            st->mode |= KFA_MODE_BINARY;
            break;
        case 'f':
            st->mode |= KFA_MODE_STDIN;
            // TODO
        case 'g':
            st->mode |= KFA_MODE_GENERIC;
            break;
        case 'G':
            st->mode |= KFA_MODE_GENONLY;
            break;
        case 'i':
            st->mode |= KFA_MODE_INVERSE;
            break;
        case 's':
            st->mode |= KFA_MODE_SHIFT;
            break;
        case 'S':
            st->mode |= KFA_MODE_SCALAR;
            break;
        case 'd': {
            st->mode |= KFA_MODE_SPARSE;
            //            parse_sparse_arg(optarg, mode);
            break;
        }
        case 'x':
            st->mode |= KFA_MODE_2D;
            break;
        case 'v':
            fprintf(stdout, "%d.%d.%d\n", VER_MAJOR, VER_MINOR, VER_PATCH);
            exit(0);
        case 'V':
            display_info(st);
            exit(0);
        case '?':
            display_help();
            exit(0);
        }
    }
// buffer_work:
//     if (optind < argc) {
//         size_t tmp_size = 1;
//         char* buf = calloc(1, tmp_size);
//         if (buf) {
//             do {
//                 char* old = buf;
//                 tmp_size += strlen(argv[optind]) + 1;
//
//                 buf = realloc(buf, tmp_size);
//                 if (buf == NULL) {
//                     free(old);
//                     goto bailout;
//                 }
//
//                 strcat(buf, argv[optind]);
//                 strcat(buf, " ");
//             } while (argc > ++optind);
//
//             buf[strlen(buf) - 1] = '\0';
//         } /* buf allocated */
//         ret = buf;
//     }
//
//     if (ret == NULL) {
//         ret = read_stdin();
//         mode->is_stdin = true;
//     }
bailout:
    return ret;
}
