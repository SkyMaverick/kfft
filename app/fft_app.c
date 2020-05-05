#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kfft.h"
#if defined(KFFT_OS_WINDOWS)
    #include "getopt_win.h"
#else
    #include <unistd.h>
    #include <getopt.h>
#endif

#include "const.h"
#include "config.h"

typedef struct {
    void* buf;

    bool is_shift;
    bool is_cpx;
    bool is_2d;
    bool is_stdin;

    size_t len;
    size_t x;
    size_t y;

    uint32_t flags;
} app_mode_t;

static size_t
parse_buffer(void** out, char* buf, bool as_cpx) {
    char* args = buf;

    size_t len = 1;
    // Analize
    while ((args = strchr(args, ' ')) != NULL)
        len++, *args = '\0', args++;

    args = buf;

    size_t alloc_size = (as_cpx) ? ((len % 2) ? ((len + 1) / 2) : (len / 2)) * sizeof(kfft_cpx)
                                 : len * sizeof(kfft_scalar);

    kfft_scalar* tmp = calloc(len, alloc_size);
    if (tmp) {
        for (size_t i = 0; i < len; i++) {
            tmp[i] = (kfft_scalar)atof(args);
            args += strlen(args) + 1;
        }
    }
    *out = tmp;
    return (as_cpx) ? ((len % 2) ? ((len + 1) / 2) : (len / 2)) : len;
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

static char*
read_stdin(void) {
    char buf[STDIN_BUF_SIZE];
    size_t ret_size = 1;
    char* ret = NULL;

    if (stdin_check() > 0) {
        ret = malloc(STDIN_BUF_SIZE * sizeof(char));
        if (ret) {
            ret[0] = '\0';
            size_t n = 0;
            while ((n = fread(buf, 1, STDIN_BUF_SIZE, stdin)) > 0) {
                char* old = ret;
                ret_size += n;
                ret = realloc(ret, ret_size + 1);
                if (ret == NULL) {
                    free(old);
                    return NULL;
                };
                strcat(ret, buf);
            }
            ret[ret_size] = '\0';
        } /* ret allocated */
    }
    return ret;
}

static void
write_stdout(kfft_scalar* in, size_t sz) {
    char buf[80];
    size_t out_size = 1;

    char* out = calloc(1, sizeof(STDOUT_BUF_SIZE));
    if (out) {
        out[0] = '\0';
        for (size_t i = 0; i < sz; i++) {
            sprintf(buf, "%.3f ", in[i]);
            out_size += strlen(buf);
            char* old = out;
            out = realloc(out, out_size);
            if (out == NULL) {
                free(old);
            }
            strcat(out, buf);
        }
        out[out_size - 2] = '\0'; // erase last space

        fprintf(stdout, "%s", out);
        fflush(stdout);

        free(out);
    } /* out allocated */
}
#else
// TODO Windows
static char*
read_stdin(void) {
    return 0;
}
static void
write_stdout(kfft_scalar* in, size_t sz) {}
#endif /* not KFFT_OS_WINDOWS */

#include "cmd.c"
#include "cwork.c"
#include "swork.c"

int
main(int argc, char* argv[]) {
    // clang-format off
    app_mode_t mode = {
            .buf = NULL,
            .is_shift = false,
            .is_cpx = true,
            .is_2d = false,
            .is_stdin = true,
            
            .len = 0,
            .x = 0,
            .y = 0,
            
            .flags = 0};
    // clang-format on

    kfft_return_t ret = KFFT_RET_SUCCESS;
    char* buffer = cmd_line_parse(argc, argv, &mode);

    if (buffer) {
        ret = (mode.is_cpx) ? work_cpx(buffer, &mode) : work_scalar(buffer, &mode);
        free(buffer);
    } else {
        display_help();
    }
    fflush(stdout);
    fclose(stdout);

    if (ret != KFFT_RET_SUCCESS)
        fprintf(stderr, "%s\n", kfft_strerr(ret));

    return ret;
}
