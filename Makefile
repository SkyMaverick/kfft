KFVER=130

ifeq ($(shell uname -s),Darwin)
	SHARED := -Wl,-install_name,libkfft.dylib -o libkfft.dylib
else
	SHARED := -Wl,-soname,libkfft.so -o libkfft.so
endif

CC=gcc

DemoApp=demo
UTestSimple=ustest
UTestCmp=ukfcmp
UTestCmpWP=ukfcmpwp
UAutoTestCmp=aukfcmp
UAutoTestCmpWP=aukfcmpwp

AppAll=${DemoApp} ${UTestSimple} ${UTestCmp} ${UTestCmpWP} ${UAutoTestCmp} ${UAutoTestCmpWP}

all: lib ${AppAll}

# ===================================

demo: lib
	${CC} ${CFLAGS} -L. -I. test/demo.c -lkfft -lm -Wl,-rpath,. -o ${DemoApp}
ustest: lib
	${CC} ${CFLAGS} -L. -I. test/test.c -lkfft -lm -Wl,-rpath,. -o ${UTestSimple}
ukfcmp: lib
	${CC} ${CFLAGS} -L. -I. test/test.c -DFFTW_COMPARE -lkfft -lfftw3 -lm -Wl,-rpath,. -o ${UTestCmp}
ukfcmpwp: lib
	${CC} ${CFLAGS} -L. -I. test/test.c -DFFTW_COMPARE -DCHECK_WITHOUT_PLAN -lkfft -lfftw3 -lm -Wl,-rpath,. -o ${UTestCmpWP}
aukfcmp: lib
	${CC} ${CFLAGS} -L. -I. test/test.c -DFFTW_COMPARE -DLOG_AUTO -lkfft -lfftw3 -lm -Wl,-rpath,. -o ${UAutoTestCmp}
aukfcmpwp: lib
	${CC} ${CFLAGS} -L. -I. test/test.c -DFFTW_COMPARE -DCHECK_WITHOUT_PLAN -DLOG_AUTO -lkfft -lfftw3 -lm -Wl,-rpath,. -o ${UAutoTestCmpWP}

lib:
	gcc ${CFLAGS} -fPIC -c *.c  -o kfft.o
	ar crus libkfft.a kfft.o
	gcc -shared ${CFLAGS} $(SHARED) kfft.o

clean:
	rm -f kfft*.tar.gz *~ *.pyc kfft*.zip *.a *.o *.so ${AppAll}

asm: kfft.s

kfft.s: kfft.c kfft.h _kfft_guts.h
	[ -e kfft.s ] && mv kfft.s kfft.s~ || true
	gcc -S kfft.c -O3 -mtune=native -ffast-math -fomit-frame-pointer -unroll-loops -dA -fverbose-asm
	gcc -o kfft_short.s -S kfft.c -O3 -mtune=native -ffast-math -fomit-frame-pointer -dA -fverbose-asm -DFIXED_POINT
	[ -e kfft.s~ ] && diff kfft.s~ kfft.s || true
