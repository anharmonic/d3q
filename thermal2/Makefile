# Makefile for D3Q
# Adapted from TDDFPT main Makefile

SUBDIRS = src minpack/lapackified thermal2
PARADIRS = minpack/distributed  
OPTDIRS = tools

default: all

ifdef SCALAPACK_LIBS
all: $(SUBDIRS)
else
all: $(SUBDIRS) $(PARADIRS)
endif

more: all $(OPTDIRS)

clean: clean_subdirs 

veryclean: clean clean_examples

clean_subdirs :
	for D in $(SUBDIRS) $(OPTDIRS); do \
	   test -d $$D && (cd $$D ; $(MAKE) $(MFLAGS) clean) \
	done

clean_examples :
	if test -d Examples ; then \
	( cd Examples ; ./clean_all) ; fi 

$(SUBDIRS) $(PARADIRS) $(OPTDIRS): FORCE
	test -d $@ && (cd $@ && $(MAKE) $(MFLAGS) all) 

MYNAME = $(shell basename $(CURDIR))

distrib: veryclean
	( cd ..; \
	rm -f d3q-latest D3Q || exit 1; \
	ln -s $(MYNAME) d3q-latest; \
	ln -s d3q-latest D3Q; \
	tar --exclude .git  --exclude \*.odt --exclude protect -cJvf $(MYNAME).tar.xz d3q-latest D3Q $(MYNAME) ) 

FORCE:
