# Makefile for tools

include ../../make.inc

STRIP = true

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../Modules $(MOD_FLAG)../../LAXlib $(MOD_FLAG)../../UtilXlib \
          $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../PW/src $(MOD_FLAG)../../LR_Modules $(MOD_FLAG)../../PHonon/PH \
          $(MOD_FLAG)../src \
          $(MOD_FLAG).

IFLAGS+=-I../../FFTXlib -I../../UtilXlib -I../../include -I../../PHonon/PH

PHOBJS = ../../PHonon/PH/libph.a ../../PHonon/PH/libphaux.a
LRMODS = ../../LR_Modules/liblrmod.a
PWOBJS = ../../PW/src/libpw.a
QEMODS = ../../Modules/libqemod.a ../../KS_Solvers/Davidson/libdavid.a ../../KS_Solvers/CG/libcg.a \
         ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a ../../UtilXlib/libutil.a ../../upflib/libupf.a


EXTOBJS = ../src/libd3q.a \
          $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS)

TLDEPS= bindir mods libs pw-lib lrmods ph

LIBOBJS = ../../clib/clib.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a


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
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/tools/$@ . )
d3_sc2c.x : d3_sc2c.o $(D3QTOOLSLIBS)
	$(LD) $(LDFLAGS) -o $@ d3_sc2c.o $(D3QTOOLSLIBS) \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/tools/$@ . )
diffd3.x : diffd3.o $(D3QTOOLSLIBS)
	$(LD) $(LDFLAGS) -o $@ diffd3.o $(D3QTOOLSLIBS) \
	$(EXTOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)
	$(STRIP)  $@
	- ( cd ../../bin ; ln -fs ../D3Q/tools/$@ . )



tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f pwi2xsf pwi2xsf_old *.x *.o *~ *.F90 *.mod *.d *.i *.L

include make.depend

