CC = mpicc
CFLAGS += -Wall -Wextra -I.. -g
LDLIBS += -lm 

CFLAGS += $(shell pkg-config --cflags scalapack-openmpi)
LDLIBS += $(shell pkg-config --libs scalapack-openmpi)

all: check

.phony: check clean
	
test_plmdif1.o ls.o: ls.h

test_plmdif1: test_plmdif1.o ls.o ../libminpack.a 
test_extrablacs: test_extrablacs.o ../libminpack.a 

check: test_extrablacs  test_plmdif1 
	prove -f --exec mpiexec ./test_extrablacs ./test_plmdif1 

clean:
	rm -f *.o test_plmdif1 test_extrablacs
