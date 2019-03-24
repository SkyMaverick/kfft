KFVER=130

ifeq ($(shell uname -s),Darwin)
	SHARED := -Wl,-install_name,libkissfft.dylib -o libkfft.dylib
else
	SHARED := -Wl,-soname,libkissfft.so -o libkfft.so
endif

AppName=app

all:
	gcc -Wall -g -fPIC -c *.c -Dkiss_fft_scalar=double -o kfft.o
	ar crus libkfft.a kfft.o
	gcc -shared $(SHARED) kfft.o
# ===================================
	gcc -Wall -g -Dkiss_fft_scalar=double -L. -I. test/testapp.c -lkfft -lm -Wl,-rpath,. -o ${AppName}
doc:
	@echo "Start by reading the README file.  If you want to build and test lots of stuff, do a 'make testall'"
	@echo "but be aware that 'make testall' has dependencies that the basic kissfft software does not."
	@echo "It is generally unneeded to run these tests yourself, unless you plan on changing the inner workings"
	@echo "of kissfft and would like to make use of its regression tests."

clean:
	rm -f kfft*.tar.gz *~ *.pyc kfft*.zip *.a *.o *.so ${AppName}

asm: kfft.s

kfft.s: kfft.c kfft.h _kfft_guts.h
	[ -e kfft.s ] && mv kfft.s kfft.s~ || true
	gcc -S kfft.c -O3 -mtune=native -ffast-math -fomit-frame-pointer -unroll-loops -dA -fverbose-asm 
	gcc -o kfft_short.s -S kfft.c -O3 -mtune=native -ffast-math -fomit-frame-pointer -dA -fverbose-asm -DFIXED_POINT
	[ -e kfft.s~ ] && diff kfft.s~ kfft.s || true
