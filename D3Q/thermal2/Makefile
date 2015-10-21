# Makefile for thermal2

include ../../make.sys

STRIP = true #strip

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
	  $(MOD_FLAG)../src $(MOD_FLAG)../../PW/src  $(MOD_FLAG).

IFLAGS=-I../../include -I../../PHonon/PH -I/opt/intel/composer_xe_2013_sp1.4.211/mkl/include/intel64/lp64/

EXTOBJS = ../src/libd3q.a \
../../PHonon/PH/libph.a \
../../PHonon/PH/libphaux.a \
../../PW/src/libpw.a

QEMODS = ../../Modules/libqemod.a

TLDEPS= bindir mods libs pw

#IBOBJS        = ../../flib/pthermal2.a 
LIBOBJS = ../../flib/flib.a ../../flib/ptools.a ../../clib/clib.a ../../iotk/src/libiotk.a
# line for xlf90 on bluegene:
#LIBOBJS = ../../flib/flib.a ../../flib/ptools.a ../../clib/clib.a ../../iotk/src/libiotk.a ../../lapack-3.2/lapack.a


# lines for xlf90 on bluegene
#F90FLAGS += -qxlf2003=polymorphic -qsmp=omp
#F90FLAGS += -qsmp=omp
# intel fortran, use -fopenmp on gfortran:
#F90FLAGS += -check all -traceback -g
#F90FLAGS += -openmp 
#LDFLAGS  += -openmp #-pg -g

D3QTOOLSLIBS = \
asr2.o \
casimir.o \
code_input.o \
dynbubble.o \
final_state.o \
functions.o \
fc2_interp.o \
fc3_interp.o \
input_fc.o \
isotopes.o \
linewidth.o \
more_constants.o \
mpi_thermal.o \
nanoclock.o \
nanosec.o \
nist_isotopes_db.o \
ph_velocity.o \
posix_signal.o \
q_grids.o \
quter.o \
timers.o \
variational_tk.o
#vector.o

EXECUTABLES = lw.x tk.x db.x asr3.x sqom.x sparse.x q2r.x r2q.x qq2rr.x

all : tldeps $(EXECUTABLES)

mpitest.x : PROGRAM_mpitest.o $(D3QTOOLSLIBS)
	$(LD) $(LDFLAGS) -o $@ PROGRAM_mpitest.o $(D3QTOOLSLIBS) \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/thermal2/$@ d3_$@ )


tk.x : PROGRAM_tk.o $(D3QTOOLSLIBS)
	$(LD) $(LDFLAGS) -o $@ PROGRAM_tk.o $(D3QTOOLSLIBS) \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/thermal2/$@ d3_$@ )

lw.x : PROGRAM_lw.o $(D3QTOOLSLIBS)
	$(LD) $(LDFLAGS) -o $@ PROGRAM_lw.o $(D3QTOOLSLIBS) \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/thermal2/$@ d3_$@ )

db.x : PROGRAM_db.o $(D3QTOOLSLIBS)
	$(LD) $(LDFLAGS) -o $@ PROGRAM_db.o $(D3QTOOLSLIBS) \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/thermal2/$@ d3_$@ )

r2q.x : PROGRAM_r2q.o $(D3QTOOLSLIBS)
	$(LD) $(LDFLAGS) -o $@ PROGRAM_r2q.o $(D3QTOOLSLIBS) \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/thermal2/$@ d3_$@ )

sparse.x : PROGRAM_sparse.o $(D3QTOOLSLIBS)
	$(LD) $(LDFLAGS) -o $@ PROGRAM_sparse.o $(D3QTOOLSLIBS) \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/thermal2/$@ d3_$@ )

sqom.x : PROGRAM_sqom.o $(D3QTOOLSLIBS)
	$(LD) $(LDFLAGS) -o $@ PROGRAM_sqom.o $(D3QTOOLSLIBS) \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/thermal2/$@ d3_$@ )

q2r.x : PROGRAM_q2r.o $(D3QTOOLSLIBS)
	$(LD) $(LDFLAGS) -o $@ PROGRAM_q2r.o $(D3QTOOLSLIBS) \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/thermal2/$@ d3_$@ )

qq2rr.x : PROGRAM_qq2rr.o $(D3QTOOLSLIBS)
	$(LD) $(LDFLAGS) -o $@ PROGRAM_qq2rr.o $(D3QTOOLSLIBS) \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/thermal2/$@ d3_$@ )

asr3.x : PROGRAM_asr3.o $(D3QTOOLSLIBS)
	$(LD) $(LDFLAGS) -o $@ PROGRAM_asr3.o $(D3QTOOLSLIBS) \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/thermal2/$@ d3_$@ )


d3q :
	( cd ../; make src)
tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.mod *.d *.i *.L

include make.depend

