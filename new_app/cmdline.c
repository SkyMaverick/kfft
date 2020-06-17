#include "info.c"

static inline bool
in_iscpx(state_t* st) {
    if (KFA_CHECK(st, SCALAR)) {
        if (KFA_CHECK(st, INVERSE))
            return true;
    } else {
        return true;
    }
    return false;
}

static inline bool
out_iscpx(state_t* st) {
    if (KFA_CHECK(st, SCALAR))
        if (KFA_CHECK(st, INVERSE))
            return false;
    return true;
}

static inline void
calculate_io(state_t* st, size_t in_size) {
    if (in_iscpx(st)) {
        if (in_size % 2)
            in_size += 1;
    }
    st->in_len = in_size;
    st->in_count = (in_iscpx(st)) ? in_size / 2 : in_size;

    st->out_count = st->in_count;
    if (out_iscpx(st)) {
        st->out_len = st->in_count * 2;
    }
}

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

static inline kfft_scalar*
parse_stdin_buffer(char* buf, state_t* st) {
    size_t count = 1;

    char* tmp = buf;
    while (tmp = strchr(tmp, ' '))
        count++, tmp++;

    calculate_io(st, count);
    kfft_scalar* ret = KRNL_FUNCS(st).cb_malloc(st->in_len * st->scsz);
    if (ret) {
        char* old;
        tmp = old = buf;

        count = 0;
        while (tmp = strchr(tmp, ' ')) {
            *tmp++ = '\0';
            ret[count] = (kfft_scalar)atof(old);

            old = tmp;
            count++;
        }
        ret[count] = (kfft_scalar)atof(old);
    }
    free(buf);
    return ret;
}
#define KFA_BUF_SIZE 0x000F
static inline kfft_scalar*
read_stdin_pipe(state_t* st) {
    char buf[KFA_BUF_SIZE];
    char* tmp = NULL;

    size_t n = 0;
    while ((n = fread(buf, 1, KFA_BUF_SIZE, stdin)) > 0) {
        tmp = realloc(tmp, st->buf.lenght + n + 1);
        if (tmp == NULL)
            return NULL;
        memcpy(&tmp[st->buf.lenght], buf, n);
        st->buf.lenght += n;
    }
    tmp[st->buf.lenght - 1] = '\0';

    return parse_stdin_buffer(tmp, st);
}

static kfft_scalar*
read_stdin_args(int argc, char* argv[], int idx, state_t* st) {
    kfft_scalar* tmp = NULL;
    if (idx < argc) {
        size_t needed = argc - idx;

        tmp = KRNL_FUNCS(st).cb_malloc(st->in_len * st->scsz);
        if (tmp) {
            memset(tmp, 0, needed * st->scsz);
            for (size_t i = 0; idx < argc; idx++, i++)
                tmp[i] = (kfft_scalar)atof(argv[idx]);
        }
        calculate_io(st, needed);
    }
    return tmp;
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
            return read_stdin_args(argc, argv, optind, st);
        }
    }
    return read_stdin_pipe(st);
}
