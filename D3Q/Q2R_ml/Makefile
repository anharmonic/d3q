.SUFFIXES : .o .f90

include ../../make.sys
F90 = $(MPIF90)

.f90.o:
	$(F90) $(F90FLAGS) -c $<

all : mat3q2r.x

#F90FLAGS = -O3 -ftrace=full -fbounds-check -g
#F90 = g95
#LIBS= -L/opt/g95/lib -llapack -lf77blas -lcblas -latlas -lrt

mat3q2r.x: mat2q2r.o libs.o
	$(F90) $(F90FLAGS) -o mat2q2r.x mat2q2r.o libs.o  $(LIBS)

clean:
	rm -rf *.o *.mod *.x
