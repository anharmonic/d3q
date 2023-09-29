include ../../../make.inc
#CFLAGS= -O3 -Wall -Wextra -I. -g #-traceback  -qopenmp -std=c11
#LDFLAGS= -g 

all: libminpack.a
default: libminpack.a

.phony: clean test bench

OBJS = enorm.o qrsolv.o fdjac2.o fdjac1.o lmpar.o dogleg.o r1updt.o r1mpyq.o \
       lmbase.o lmdif1.o lmdif.o lmder.o lmder1.o \
       hybrbase.o hybrd.o hybrj.o hybrd1.o hybrj1.o \
       chkder.o

$(OBJS): minpack.h

libminpack.a:  $(OBJS)
	ar r $@ $(OBJS); ranlib $@

check: libminpack.a
	cd tests ; $(MAKE) check

bench: libminpack.a
	cd benchmarks ; $(MAKE) bench

clean:
	rm -f *.o libminpack.a
