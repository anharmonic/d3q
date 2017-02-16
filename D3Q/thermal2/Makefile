# Makefile for thermal2

include ../../make.inc

STRIP = true #strip

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
	  $(MOD_FLAG)../src $(MOD_FLAG)../../PW/src  $(MOD_FLAG).

IFLAGS+=-I../../FFTXlib -I../../include -I../../PHonon/PH 
#F90FLAGS += -check all -traceback -g
# For profiling
#F90FLAGS += -p
#LDFLAGS += -p

# Remove the following line if you want to use intel compiler but not MKL 
# vectorized vdCos and vdSin
DFLAGS+= -D__HASVTRIG

# Uncomment the following two lines to enable openmp in thermal2
# use "-fopenmp" for gfortran and xlf compiler, "-openmp" for intel compiler
#F90FLAGS+=-fopenmp
#LDFLAGS+=-fopenmp

EXTOBJS = ../src/libd3q.a \
../../PHonon/PH/libph.a \
../../PHonon/PH/libphaux.a \
../../PW/src/libpw.a

QEMODS = ../../Modules/libqemod.a

TLDEPS= bindir mods libs pw-lib lrmods ph

LIBOBJS = ../../LR_Modules/liblrmod.a ../../clib/clib.a ../../iotk/src/libiotk.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a


# lines for xlf90 on bluegene
#F90FLAGS += -qxlf2003=polymorphic -qsmp=omp
#F90FLAGS += -qsmp=omp
# intel fortran, use -fopenmp on gfortran:
#F90FLAGS += -check all -traceback -g
#F90FLAGS += -openmp 
#LDFLAGS  += -openmp #-pg -g

THERMALK2_OBJS = \
mpi_thermal.o \
asr2.o \
casimir.o \
code_input.o \
dynbubble.o \
final_state.o \
functions.o \
fc2_interp.o \
f3_bwfft.o \
fc3_interp.o \
input_fc.o \
isotopes.o \
linewidth.o \
merge_degenerate.o \
more_constants.o \
nanoclock.o \
nanosec.o \
nist_isotopes_db.o \
ph_system.o \
ph_velocity.o \
posix_signal.o \
q_grids.o \
quter.o \
timers.o \
variational_tk.o

EXECUTABLES = lw.x tk.x asr3.x sqom.x sparse.x q2r.x r2q.x qq2rr.x import3py.x recenter.x

all : tldeps d3q-libs thermal2-libs $(EXECUTABLES)

tk.x : PROGRAM_tk.o tldeps d3q-libs thermal2-libs
	$(LD) $(LDFLAGS) -o $@ PROGRAM_tk.o libthermal2.a \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/thermal2/$@ d3_$@ )

lw.x : PROGRAM_lw.o tldeps d3q-libs thermal2-libs
	$(LD) $(LDFLAGS) -o $@ PROGRAM_lw.o libthermal2.a \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/thermal2/$@ d3_$@ )

db.x : PROGRAM_db.o tldeps d3q-libs thermal2-libs
	$(LD) $(LDFLAGS) -o $@ PROGRAM_db.o libthermal2.a \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/thermal2/$@ d3_$@ )

r2q.x : PROGRAM_r2q.o tldeps d3q-libs thermal2-libs
	$(LD) $(LDFLAGS) -o $@ PROGRAM_r2q.o libthermal2.a \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/thermal2/$@ d3_$@ )

sparse.x : PROGRAM_sparse.o tldeps d3q-libs thermal2-libs
	$(LD) $(LDFLAGS) -o $@ PROGRAM_sparse.o libthermal2.a \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/thermal2/$@ d3_$@ )

sqom.x : PROGRAM_sqom.o tldeps d3q-libs thermal2-libs
	$(LD) $(LDFLAGS) -o $@ PROGRAM_sqom.o libthermal2.a \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/thermal2/$@ d3_$@ )

q2r.x : PROGRAM_q2r.o tldeps d3q-libs thermal2-libs
	$(LD) $(LDFLAGS) -o $@ PROGRAM_q2r.o libthermal2.a \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/thermal2/$@ d3_$@ )

qq2rr.x : PROGRAM_qq2rr.o tldeps d3q-libs thermal2-libs
	$(LD) $(LDFLAGS) -o $@ PROGRAM_qq2rr.o libthermal2.a \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/thermal2/$@ d3_$@ )

asr3.x : PROGRAM_asr3.o tldeps d3q-libs thermal2-libs
	$(LD) $(LDFLAGS) -o $@ PROGRAM_asr3.o libthermal2.a \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/thermal2/$@ d3_$@ )

import3py.x : PROGRAM_import3py.o tldeps d3q-libs thermal2-libs
	$(LD) $(LDFLAGS) -o $@ PROGRAM_import3py.o libthermal2.a \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/thermal2/$@ d3_$@ )

recenter.x : PROGRAM_recenter.o tldeps d3q-libs thermal2-libs
	$(LD) $(LDFLAGS) -o $@ PROGRAM_recenter.o libthermal2.a \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/thermal2/$@ d3_$@ )

d3q-libs :
	( cd ../src/; make libd3q.a)

thermal2-libs : libthermal2.a

libthermal2.a : $(THERMALK2_OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

#d3q :
#	( cd ../; make src)

tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.mod *.d *.i *.L

include make.depend

