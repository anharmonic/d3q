# Makefile for tools

include ../../make.sys

STRIP = true #strip

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
	  $(MOD_FLAG)../src $(MOD_FLAG)../../PW/src  $(MOD_FLAG).  

IFLAGS=-I../../include

EXTOBJS = ../src/libd3q.a \
../../PHonon/PH/libph.a \
../../PW/src/libpw.a 

QEMODS = ../../Modules/libqemod.a

TLDEPS= bindir mods libs pw

LIBOBJS        = ../../flib/ptools.a ../../flib/flib.a ../../clib/clib.a ../../iotk/src/libiotk.a

#D3QTOOLSOBJS = \
#lw.o \

#MPIF90 = ifort
F90FLAGS += -fopenmp -ffpe-trap=invalid,zero,overflow -fcheck=all -fbacktrace -g
LDFLAGS  += -fopenmp

D3QTOOLSLIBS = \
asr2.o \
input_fc.o \
interp_fc.o \
nanoclock.o \
nanosec.o \
ph_velocity.o \
q_grid.o \
functions.o \

EXECUTABLES = lw.x asr3.x

all : tldeps $(EXECUTABLES)

lw.x : lw.o $(D3QTOOLSLIBS) d3q
	$(LD) $(LDFLAGS) -o $@ lw.o $(D3QTOOLSLIBS) \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS) 
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/tools/$@ . )

asr3.x : asr3.o $(D3QTOOLSLIBS) d3q
	$(LD) $(LDFLAGS) -o $@ asr3.o $(D3QTOOLSLIBS) \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS) 
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/tools/$@ . )


d3q :
	( cd ../; make src)
tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.mod *.d *.i *.L

include make.depend

