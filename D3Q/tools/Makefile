# Makefile for tools

include ../../make.sys

STRIP = strip

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
	  $(MOD_FLAG)../src $(MOD_FLAG).  

EXTOBJS = ../src/libd3q.a \
../../PHonon/PH/libph.a \
../../PW/src/libpw.a 

QEMODS = ../../Modules/libqemod.a

TLDEPS= bindir mods libs pw 

LIBOBJS        = ../../flib/flib.a ../../clib/clib.a ../../iotk/src/libiotk.a ../../FFTXlib/libqefft.a

D3QTOOLSOBJS = \
diffd3.o \
d3_sc2c.o \
xml2giorgia.o

D3QTOOLSLIBS = \
read_d3dyn.o \
sc2c_tools.o 

#sc2c.o 
#sc3c.o 
EXECUTABLES = xml2giorgia.x d3_sc2c.x diffd3.x

all : tldeps $(EXECUTABLES)

xml2giorgia.x : xml2giorgia.o $(D3QTOOLSLIBS)
	$(LD) $(LDFLAGS) -o $@ xml2giorgia.o $(D3QTOOLSLIBS) \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS) 
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/tools/$@ . )
d3_sc2c.x : d3_sc2c.o $(D3QTOOLSLIBS)
	$(LD) $(LDFLAGS) -o $@ d3_sc2c.o $(D3QTOOLSLIBS) \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS) 
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/tools/$@ . )
diffd3.x : diffd3.o $(D3QTOOLSLIBS)
	$(LD) $(LDFLAGS) -o $@ diffd3.o $(D3QTOOLSLIBS) \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(LIBS) 
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/tools/$@ . )



tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f pwi2xsf pwi2xsf_old *.x *.o *~ *.F90 *.mod *.d *.i *.L

include make.depend

