#include "info.c"

static int
parse_2d_arg(char* arg, state_t* st) {
    int ret = atoi(arg);
    if (ret > 0) {
        st->dims.x = ret;
        return 0;
    }
    return 1;
}

static int
parse_sparse_arg(char* arg, state_t* st) {
    size_t olen = strlen(arg);
    char* buf = calloc(olen + 1, sizeof(char));
    if (buf == NULL)
        return 1;

    strncpy(buf, optarg, olen);

    char* s = strchr(buf, ':');
    if (s)
        *s = '\0';

    st->sparse.dx = atoi(buf);
    if (strlen(buf) < olen)
        st->sparse.sx = atoi(s + 1);

    free(buf);
    return 0;
}

#if !defined(KFFT_OS_WINDOWS)
int
stdin_check(void) {
    fd_set rd;
    struct timeval tv = {1, 0};
    int ret;

    FD_ZERO(&rd);
    FD_SET(STDIN_FILENO, &rd);
    ret = select(1, &rd, NULL, NULL, &tv);

    return (ret > 0);
}
#else
int
stdin_check(void) {
    // FIXME
    return 1;
}
#endif /* not KFFT_OS_WINDOWS */

void*
realloc_align(void* mem, size_t old_sz, size_t new_sz, state_t* st) {
    void* ret = KRNL_FUNCS(st).cb_malloc(new_sz);
    if (ret) {
        memset(ret, 0, new_sz);
        memcpy(ret, mem, old_sz);
    }

    free(mem);
    return ret;
}

static inline char*
pipe_read_stdin(state_t* st) {
    char buf[STDIN_BUF_SIZE];
    size_t ret_size = 1;
    char* ret = NULL;

    if (stdin_check() > 0) {
        ret = calloc(STDIN_BUF_SIZE, sizeof(kfft_scalar));
        if (ret) {
            size_t n = 0;
            while ((n = fread(buf, 1, STDIN_BUF_SIZE, stdin)) > 0) {
                size_t new_size = (st->mode & KFA_MODE_BINARY)
                                      ? (ret_size + n) + (ret_size + n) % sizeof(kfft_scalar)
                                      : ret_size + n;
                if ((ret = realloc_align(ret, ret_size, new_size + 1, st)) != NULL) {
                    strcat(ret, buf);
                    ret_size += new_size;
                } else
                    return NULL;
            }
            ret[ret_size] = '\0';
        } /* ret allocated */
    }
    st->buf.lenght = ret_size;

    return ret;
}

static inline char*
manual_read_stdin(int argc, char* argv[], int idx) {
    char* buf = NULL;
    if (idx < argc) {
        size_t tmp_size = 1;
        buf = calloc(1, tmp_size);
        if (buf) {
            do {
                char* old = buf;
                tmp_size += strlen(argv[idx]) + 1;

                buf = realloc(buf, tmp_size);
                if (buf == NULL) {
                    free(old);
                    return NULL;
                }

                strcat(buf, argv[idx]);
                strcat(buf, " ");
            } while (argc > ++idx);

            buf[strlen(buf) - 1] = '\0';
        } /* buf allocated */
    }
    return buf;
}

static inline unsigned
post_buffer_process(state_t* st) {
    if (st->mode & KFA_MODE_SCALAR) {
        if (st->mode & KFA_MODE_INVERSE) {
            st->out_lenght = (st->lenght % 2) ? (st->lenght + 1) / 2 : st->lenght / 2;
        } else {
            st->out_lenght = st->lenght * 2;
        }
    } else {
        st->out_lenght = (st->lenght % 2) ? st->lenght + 1 : st->lenght;
    }

    if (st->mode & KFA_MODE_2D) {
        if (st->lenght % st->dims.x)
            return KFA_RET_FAIL_ARGS;
    }
    if (st->mode & KFA_MODE_SPARSE) {
        if (st->lenght % (st->sparse.dx + st->sparse.sx))
            return KFA_RET_FAIL_ARGS;
    }
    printf("%zu\n", st->lenght);
    printf("%zu\n", st->out_lenght);
    return KFA_RET_SUCCESS;
}

static inline kfft_scalar*
b2s_binary(char* buffer, state_t* st) {
    return (kfft_scalar*)buffer;
}

static inline kfft_scalar*
b2s_manual(char* buffer, state_t* st) {
    size_t len = 1;

    if (buffer == NULL)
        return NULL;

    char* args = buffer;
    // Analize
    while ((args = strchr(args, ' ')) != NULL)
        len++, *args = '\0', args++;
    args = buffer;

    kfft_scalar* tmp = KRNL_FUNCS(st).cb_malloc((len + 1) * sizeof(kfft_scalar));
    if (tmp) {
        for (size_t i = 0; i < len; i++) {
            tmp[i] = (kfft_scalar)atof(args);
            args += strlen(args) + 1;

            st->lenght++;
        }
    }
    free(buffer);

    if (len % 2)
        st->lenght += 1;

    if (post_buffer_process(st) != KFA_RET_SUCCESS) {
        free(tmp);
        tmp = NULL;
    }

    return tmp;
}

static kfft_scalar*
buf2scalar(char* buffer, state_t* st) {
    return (st->mode & KFA_MODE_BINARY) ? b2s_binary(buffer, st) : b2s_manual(buffer, st);
}

static kfft_scalar*
cmd_line_parse(int argc, char* argv[], state_t* st) {
    int opt = 0;
    while ((opt = getopt(argc, argv, FMT_OPTSTRING)) != -1) {
        switch (opt) {
        case 'b':
            st->mode |= KFA_MODE_BINARY;
            break;
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
            if (parse_sparse_arg(optarg, st))
                return NULL;
            st->mode |= KFA_MODE_SPARSE;
            break;
        }
        case 'x':
            if (parse_2d_arg(optarg, st))
                return NULL;
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
        case 'f':
            st->mode |= KFA_MODE_STDIN;
            return buf2scalar(manual_read_stdin(argc, argv, optind), st);
        }
    }
    return buf2scalar(pipe_read_stdin(st), st);
}
