#define STDIN_BUF_SIZE 0xFF
#define STDOUT_BUF_SIZE 0xFF

const char* const help_msg = " \
Fast Fourier Transform for input sequence. Use kfft library.\n\
Use program:\n\
\tkfft [gGidsSx:f:vV?] <sequence>\n\
Input/output buffer format - float / double numbers separated by spaces \n\
Input format: \n\
\t- kfft <args> x0 x1 x2 .. xN (for scalar sequence) \n\
\t- kfft <args> x0.r x0.i x1.r x1.i .. xN.r xN.i (for complex sequence) \n\
Output format: \n\
\t- x0 x1 x2 .. xN (for scalar sequence) \n\
\t- x0.r x0.i x1.r x1.i .. xN.r xN.i (for complex sequence) \n\
";
