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
KfftSrc=kfft.c kfft.h _kfft_guts.h _kfft_bf.h
KfftObjs=kfft.o kfft_core.o

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

%.o: %.c
	 ${CC} ${CFLAGS} ${LIBFLAGS} -fPIC -c $< -o $@
    
lib: kfft.o kfft_core.o
	ar crus libkfft.a ${KfftObjs}
	gcc -shared ${CFLAGS} $(SHARED) ${KfftObjs}

clean:
	rm -f kfft*.tar.gz *~ *.pyc kfft*.zip *.a *.o *.so *.s ${AppAll}
