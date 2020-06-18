#include "config.h"

#define FMT_OPTSTRING "gGisSd:x:vV?"

#define STDIN_BUF_SIZE 0xFF
#define STDOUT_BUF_SIZE 0xFF

#define __HELP_SRING                                                                               \
    "Fast Fourier Transform for input sequence. Use kfft library.\n"                               \
    "Use as part of an operations pipeline:\n"                                                     \
    "\t... <sequence> | " APP_NAME " [-" FMT_OPTSTRING "] | " APP_NAME " [-" FMT_OPTSTRING         \
    "] | <sequence>...\n"                                                                          \
    "\nInput/output buffer format - float / double numbers separated by spaces \n"                 \
    "Input format: \n"                                                                             \
    "\t- " APP_NAME " <args> x0 x1 x2 .. xN (for scalar sequence) \n"                              \
    "\t- " APP_NAME " <args> x0.r x0.i x1.r x1.i .. xN.r xN.i (for complex sequence) \n"           \
    "Output format: \n"                                                                            \
    "\t- x0 x1 x2 .. xN (for scalar sequence) \n"                                                  \
    "\t- x0.r x0.i x1.r x1.i .. xN.r xN.i (for complex sequence) \n"                               \
    "\nArguments: \n"                                                                              \
    "\t-g - use generic algorithm for prime size chuncks evalation \n"                             \
    "\t-G - use ONLY generic algorithm for all sequence evalation \n"                              \
    "\t-i - use inverse transformation \n"                                                         \
    "\t-S - use scalar functions for evalate \n"                                                   \
    "\t-d - use sparse functionality (arguments dims:step)\n"                                      \
    "\t-s - use shift functionality \n"                                                            \
    "\t-x - translate linear buffer as 2d buffer X * (lenght / X) \n"                              \
    "\t-v - version info \n"                                                                       \
    "\t-V - extended version info \n"

#define HELP_CONST __HELP_SRING

const char* const help_msg = HELP_CONST;
