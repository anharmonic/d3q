# Makefile for tools

include ../../make.sys

STRIP = true #strip

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
	  $(MOD_FLAG)../src $(MOD_FLAG)../../PW/src  $(MOD_FLAG).  

IFLAGS=-I../../include -I../../PHonon/PH

EXTOBJS = ../src/libd3q.a \
../../PHonon/PH/libph.a \
../../PW/src/libpw.a 

QEMODS = ../../Modules/libqemod.a

TLDEPS= bindir mods libs pw

LIBOBJS        = ../../flib/ptools.a ../../flib/flib.a ../../clib/clib.a ../../iotk/src/libiotk.a

#D3QTOOLSOBJS = \
#lw.o \

F90FLAGS += -check all -traceback -g
#F90FLAGS += -openmp #-g -axAVX -vec-report=1 -fp-model fast
#-xAVX 
LDFLAGS  += -openmp #-g
#F90FLAGS +=  -g -profile-loops=all
#LDFLAGS  +=  -g

D3QTOOLSLIBS = \
asr2.o \
final_state.o \
functions.o \
input_fc.o \
interp_fc.o \
linewidth.o \
more_constants.o \
nanoclock.o \
nanosec.o \
ph_velocity.o \
q_grid.o \
quter.o \
quter_standalone.o \
sparse_fc.o \
vector.o 

EXECUTABLES = lw.x asr3.x sqom.x gen_sparse.x q2r.x

all : tldeps $(EXECUTABLES)

lw.x : lw.o $(D3QTOOLSLIBS) d3q
	$(LD) $(LDFLAGS) -o $@ lw.o $(D3QTOOLSLIBS) \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS) 
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/tools/$@ . )

gen_sparse.x : gen_sparse.o $(D3QTOOLSLIBS) d3q
	$(LD) $(LDFLAGS) -o $@ gen_sparse.o $(D3QTOOLSLIBS) \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS) 
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/tools/$@ . )


sqom.x : sqom.o $(D3QTOOLSLIBS) d3q
	$(LD) $(LDFLAGS) -o $@ sqom.o $(D3QTOOLSLIBS) \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS) 
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/tools/$@ . )

q2r.x : q2r.o $(D3QTOOLSLIBS) d3q
	$(LD) $(LDFLAGS) -o $@ q2r.o $(D3QTOOLSLIBS) \
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

