CC = mpicc
CFLAGS = -Wall -Wextra -I..
LDLIBS = -lm

CFLAGS += $(shell pkg-config --cflags scalapack-openmpi)
LDLIBS += $(shell pkg-config --libs scalapack-openmpi)

all: bench

bench: speed_lmdif1 speed_lmdif1_cross
	mpiexec ./speed_lmdif1 ./speed_lmdif1_cross
	
.phony: bench clean

speed_lmdif1: speed_lmdif1.o ls.o             ../libminpack.a
speed_lmdif1_cross: speed_lmdif1_cross.o eq.o ../libminpack.a

clean:
	rm -f *.o
	rm -f speed_lmdif1 speed_lmdif1_cross
