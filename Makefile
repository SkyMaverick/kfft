KFVER=130

ifeq ($(shell uname -s),Darwin)
	SHARED := -Wl,-install_name,libkissfft.dylib -o libkissfft.dylib
else
	SHARED := -Wl,-soname,libkissfft.so -o libkissfft.so
endif

all:
	gcc -Wall -g -fPIC -c *.c -Dkiss_fft_scalar=double -o kiss_fft.o
	ar crus libkissfft.a kiss_fft.o
	gcc -shared $(SHARED) kiss_fft.o
# ===================================
	gcc -Wall -g -Dkiss_fft_scalar=double -L. -I. test/testapp.c -lkissfft -lm -Wl,-rpath,. -o app
doc:
	@echo "Start by reading the README file.  If you want to build and test lots of stuff, do a 'make testall'"
	@echo "but be aware that 'make testall' has dependencies that the basic kissfft software does not."
	@echo "It is generally unneeded to run these tests yourself, unless you plan on changing the inner workings"
	@echo "of kissfft and would like to make use of its regression tests."

clean:
	rm -f kiss_fft*.tar.gz *~ *.pyc kiss_fft*.zip 

asm: kiss_fft.s

kiss_fft.s: kiss_fft.c kiss_fft.h _kiss_fft_guts.h
	[ -e kiss_fft.s ] && mv kiss_fft.s kiss_fft.s~ || true
	gcc -S kiss_fft.c -O3 -mtune=native -ffast-math -fomit-frame-pointer -unroll-loops -dA -fverbose-asm 
	gcc -o kiss_fft_short.s -S kiss_fft.c -O3 -mtune=native -ffast-math -fomit-frame-pointer -dA -fverbose-asm -DFIXED_POINT
	[ -e kiss_fft.s~ ] && diff kiss_fft.s~ kiss_fft.s || true
