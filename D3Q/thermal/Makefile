#.SUFFIXES : .o .f90
F90 = $(MPIF90)

.f90.o:
	$(F90) $(F90FLAGS) -c $<

include ../../make.sys

all : thermalk.x
# for intel ifort:
#F90FLAGS = -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow
#F90 = mpif90
#LIBS= -mkl=parallel  -lrt 

# for g95 with ATLAS
#F90FLAGS = -g -O3 -Wunused
#F90 = mpif90
#LIBS=-L/opt/g95/lib -llapack  -lf77blas -lcblas -latlas  -lrt


OBJS=thermalk.o linpack.o libs.o nanosec.o mat_input.o modules.o
thermalk.x: $(OBJS)
	$(F90) -o thermalk.x  $(OBJS)   $(LIBS)  $(F90FLAGS)

clean:
	rm -f *.o *.mod *.x

#dependencies
thermalk.o: modules.o
thermalk.o: mat_input.o
mat_input.o: modules.o

