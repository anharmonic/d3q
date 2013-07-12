# Makefile for D3Q
# Adapted from TDDFPT main Makefile

SUBDIRS = src tools
OPTDIRS = Q2R3 Q2R_ml thermal thermal2  

default: all

all: $(SUBDIRS)
more: all $(OPTDIRS)

clean: clean_subdirs clean_examples

clean_subdirs :
	for D in $(SUBDIRS) $(OPTDIRS); do \
	   test -d $$D && (cd $$D ; $(MAKE) $(MFLAGS) clean) \
	done

clean_examples :
	if test -d examples ; then \
	( cd examples ; ./clean_all) ; fi 

$(SUBDIRS) $(OPTDIRS): FORCE
	test -d $@ && (cd $@ ; $(MAKE) $(MFLAGS) all) 

FORCE:
