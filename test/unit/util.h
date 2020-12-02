/* Modules for CUnit test implementation
 (c) a_dobkin (habrahabr.ru)
 */

#ifndef _UTIL_H_
#define _UTIL_H_

#include <CUnit/Basic.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define _TEST(name) static void test_##name()

#define _ADD_TEST(suite, name)                                                                     \
    if (NULL == CU_add_test(suite, #name, (CU_TestFunc)test_##name)) {                             \
        CU_cleanup_registry();                                                                     \
        return;                                                                                    \
    }

#define _MSG_OFF(out, err)                                                                         \
    do {                                                                                           \
        dup2(STDOUT_FILENO, (out));                                                                \
        dup2(STDERR_FILENO, (err));                                                                \
        close(STDOUT_FILENO);                                                                      \
        close(STDERR_FILENO);                                                                      \
    } while (0)

#define _MSG_ON(out, err)                                                                          \
    do {                                                                                           \
        dup2((out), STDOUT_FILENO);                                                                \
        dup2((err), STDERR_FILENO);                                                                \
        close((out));                                                                              \
        close((err));                                                                              \
    } while (0)

#define STDIO_OFF(id)                                                                              \
    int _mfdout_##id = 0, _mfderr_##id = 0;                                                        \
    _MSG_OFF(_mfdout_##id, _mfderr_##id);

#define STDIO_ON(id) _MSG_ON(_mfdout_##id, _mfderr_##id);

CU_pSuite
CUnitCreateSuite(const char* title);
void
CUnitInitialize(void);
void
CUnitUInitialize(void);

#endif
