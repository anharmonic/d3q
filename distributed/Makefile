CC = mpicc
CFLAGS= -O3 -Wall -Wextra -I. -g
LDFLAGS= -g

all: default
default: libminpack.a

.phony: clean test bench

OBJS = utils.o extrablacs.o enorm.o qrsolv.o lmpar.o plmdif1.o plmdif.o

$(OBJS): pminpack.h

libminpack.a:  $(OBJS)
	ar r $@ $(OBJS); ranlib $@

check: libminpack.a
	cd tests ; $(MAKE) check

bench: libminpack.a
	cd benchmarks ; $(MAKE) bench

clean:
	rm -f *.o libminpack.a
	(cd tests ; $(MAKE) clean)
	(cd benchmarks ; $(MAKE) clean)
