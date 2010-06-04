# Uncomment one pair of the following for SGI Irix 6.x complib.sgimath
# IBM ESSL, or HP veclib (Convex MLIB). The default will use
# a relatively slow built-in FFT from Numerical Recipes
#
FFTFLAGS=-DSGI_COMPLIB
FFTLIBS=-lcomplib.sgimath -lftn
#
#FFTFLAGS=-DIBM_ESSL
#FFTLIBS=-lessl
#
#FFTFLAGS=-DHP_VECLIB
#FFTLIBS=-lveclib
#
#FFTFLAGS=
#FFTLIBS=

CFLAGS= -g
CC=mpicc

tests:: benchfft.exe testfft.exe

.c.o:
	$(CC) -c $(CFLAGS) $?

OBJECTS = mpi3dfft.o parautil.o fftwrap.o prf.o

fftwrap.o : fftwrap.c
	$(CC) -c $(CFLAGS) $(FFTFLAGS) fftwrap.c


benchfft.exe : benchfft.o $(OBJECTS)
	$(CC) $(CFLAGS) benchfft.o $(OBJECTS) -o benchfft.exe $(FFTLIBS) -lm

testfft.exe : testfft.o $(OBJECTS)
	$(CC) $(CFLAGS) testfft.o $(OBJECTS) -o testfft.exe $(FFTLIBS) -lm



