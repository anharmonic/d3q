# Makefile for thermal2

include ../../make.inc

STRIP = true #strip

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
	  $(MOD_FLAG)../src $(MOD_FLAG)../../PW/src  $(MOD_FLAG).

IFLAGS+=-I../../FFTXlib -I../../include -I../../PHonon/PH 

# Remove the following line if you want to use intel compiler but not MKL 
# vectorized vdCos and vdSin
DFLAGS+= -D__HASVTRIG

# lines for xlf90 on bluegene
#F90FLAGS += -qxlf2003=polymorphic -qsmp=omp
#F90FLAGS += -qsmp=omp
# intel fortran, use -fopenmp on gfortran:
#F90FLAGS += -check all -traceback -g
#F90FLAGS += -openmp 
#LDFLAGS  += -openmp #-pg -g

# For profiling
#F90FLAGS += -p
#LDFLAGS += -p

#F90FLAGS += -g77libs
#F90FLAGS += -Mbounds -Mchkptr -Mchkstk -traceback


EXTOBJS = \
../src/libd3q.a \
../../PHonon/PH/libph.a \
../../PHonon/PH/libphaux.a \
../../LR_Modules/liblrmod.a \
../../PW/src/libpw.a

QEMODS = ../../Modules/libqemod.a

TLDEPS= bindir mods libs pw-lib lrmods ph

LIBOBJS = ../../clib/clib.a ../../iotk/src/libiotk.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a

THERMALK2_OBJS = \
mpi_thermal.o \
asr2.o \
casimir.o \
code_input.o \
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

EXECUTABLES = d3_lw.x d3_tk.x d3_asr3.x d3_sqom.x d3_sparse.x d3_q2r.x d3_r2q.x d3_qq2rr.x d3_import3py.x d3_recenter.x

MYNAME       = $(shell basename $(CURDIR))
FATHERNAME   = $(shell cd ..; basename `pwd -P`)
THELINKPLACE = $(FATHERNAME)/$(MYNAME)

all : tldeps d3q-libs thermal2-libs $(EXECUTABLES)

d3_tk.x : PROGRAM_tk.o tldeps d3q-libs thermal2-libs
	$(LD) $(LDFLAGS) -o $@ PROGRAM_tk.o libthermal2.a \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../$(THELINKPLACE)/$@ $@ )

d3_lw.x : PROGRAM_lw.o tldeps d3q-libs thermal2-libs
	$(LD) $(LDFLAGS) -o $@ PROGRAM_lw.o libthermal2.a \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../$(THELINKPLACE)/$@ $@ )

#db.x : PROGRAM_db.o tldeps d3q-libs thermal2-libs
#	$(LD) $(LDFLAGS) -o $@ PROGRAM_db.o libthermal2.a \
#	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
#	$(STRIP)  $@
#	- ( cd ../../bin ; ln -fs ../$(THELINKPLACE)/$@ d3_$@ )

d3_r2q.x : PROGRAM_r2q.o tldeps d3q-libs thermal2-libs
	$(LD) $(LDFLAGS) -o $@ PROGRAM_r2q.o libthermal2.a \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../$(THELINKPLACE)/$@ $@ )

d3_sparse.x : PROGRAM_sparse.o tldeps d3q-libs thermal2-libs
	$(LD) $(LDFLAGS) -o $@ PROGRAM_sparse.o libthermal2.a \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../$(THELINKPLACE)/$@ $@ )

d3_sqom.x : PROGRAM_sqom.o tldeps d3q-libs thermal2-libs
	$(LD) $(LDFLAGS) -o $@ PROGRAM_sqom.o libthermal2.a \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../$(THELINKPLACE)/$@ $@ )

d3_q2r.x : PROGRAM_q2r.o tldeps d3q-libs thermal2-libs
	$(LD) $(LDFLAGS) -o $@ PROGRAM_q2r.o libthermal2.a \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../$(THELINKPLACE)/$@ $@ )

d3_qq2rr.x : PROGRAM_qq2rr.o tldeps d3q-libs thermal2-libs
	$(LD) $(LDFLAGS) -o $@ PROGRAM_qq2rr.o libthermal2.a \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../$(THELINKPLACE)/$@ $@ )

d3_asr3.x : PROGRAM_asr3.o tldeps d3q-libs thermal2-libs
	$(LD) $(LDFLAGS) -o $@ PROGRAM_asr3.o libthermal2.a \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../$(THELINKPLACE)/$@ $@ )

d3_import3py.x : PROGRAM_import3py.o tldeps d3q-libs thermal2-libs
	$(LD) $(LDFLAGS) -o $@ PROGRAM_import3py.o libthermal2.a \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../$(THELINKPLACE)/$@ $@ )

d3_recenter.x : PROGRAM_recenter.o tldeps d3q-libs thermal2-libs
	$(LD) $(LDFLAGS) -o $@ PROGRAM_recenter.o libthermal2.a \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../$(THELINKPLACE)/$@ $@ )

d3q-libs :
	( cd ../src/; make libd3q.a)

thermal2-libs : libthermal2.a

libthermal2.a : $(THERMALK2_OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.mod *.d *.i *.L

include make.depend

