.SUFFIXES : .o .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

all : mat3q2r.x

include ../../make.sys
F90 = $(MPIF90)

#F90FLAGS = -O3 -ftrace=full -fbounds-check
#F90 = g95
#LIBS= -L/opt/g95/lib -llapack -lf77blas -lcblas -latlas -lrt

mat3q2r.x: mat3q2r.o libs.o
	$(F90) -o mat3q2r.x  mat3q2r.o libs.o  $(LIBS)  $(F90FLAGS)

clean:
	rm -rf *.o *.mod *.x
