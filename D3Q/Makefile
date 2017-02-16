# Makefile for D3Q
# Adapted from TDDFPT main Makefile

SUBDIRS = src thermal2  
OPTDIRS = 

default: all

all: $(SUBDIRS)
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

$(SUBDIRS) $(OPTDIRS): FORCE
	test -d $@ && (cd $@ ; $(MAKE) $(MFLAGS) all) 

distrib: veryclean
	( cd ..; \
	tar --exclude .svn -czvf d3q.tgz D3Q; \
	tar --exclude .svn -cJvf d3q.tar.xz D3Q )

FORCE:
