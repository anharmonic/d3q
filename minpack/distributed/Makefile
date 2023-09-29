include ../../../make.inc
MPICC = mpicc
#CFLAGS= -O3 -Wall -Wextra -I. -g
#LDFLAGS= -g

.SUFFIXES :
.SUFFIXES : .o .c

# most fortran compilers can directly preprocess c-like directives: use
#       $(MPIF90) $(F90FLAGS) -c $<
# if explicit preprocessing by the C preprocessor is needed, use:
#       $(CPP) $(CPPFLAGS) $< -o $*.F90
#       $(MPIF90) $(F90FLAGS) -c $*.F90 -o $*.o
# remember the tabulator in the first column !!!


.c.o:
	$(MPICC) $(CFLAGS)  -c $<


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
