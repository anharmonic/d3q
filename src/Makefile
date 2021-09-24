# Makefile for D3Q

include ../../make.inc
include make.libs

#
# use recursive definitions in Makefile to modify the PATH for modules
#

MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../PHonon/PH \
          $(MOD_FLAG)../../LR_Modules \
          $(MOD_FLAG)../../dft-d3

#
# The above is in alphabetical order except when order matter during compilation
#
##### io_scattering.o \

D3QOBJS = \
allocate_d3.o \
cmdline_param.o \
d2mxc.o \
d3_basis.o \
d3com.o \
d3_debug.o \
d3_exc_gc.o \
d3_exc.o \
d3_grid.o \
d3_h_psi.o \
d3_init.o \
d3_iofiles.o \
d3ionq.o \
d3matrix_io2.o \
d3matrix.o \
d3_nlcc.o \
d3_open.o \
d3q.o \
d3_readin.o \
d3_recover.o \
d3_reset.o \
d3_restart.o \
d3_setup.o \
d3_shuffle.o \
d3_symmetry.o \
d3_valence.o \
dpsi1dpsi2dv3.o \
dpsi1dv2dpsi3.o \
dq1rhodq23v.o \
dq_vscf.o \
drho_add_phase.o \
dvdpsi.o \
efermi_shift.o \
gen_dpdvp.o \
h_psiq.o \
incdrhoscf2.o \
print_clock_d3.o \
qstar_d3.o \
rhodq123v.o \
rotate_d3.o \
run_nscf_d3.o \
set_kplus3q.o \
setup_d3gc.o \
solve_linter_d3q.o \
stop_d3.o \
w_1gauss.o \
write_d3dyn.o
#default : d3q


PHOBJS = ../../PHonon/PH/libph.a
LRMODS = ../../LR_Modules/liblrmod.a
PWOBJS = ../../PW/src/libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a
QEMODS = $(BASEMODS)

TLDEPS= pw ph pp

all :   pw ph pp libd3q.a d3q.x

libd3q.a : $(D3QOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

d3q.x : d3q.o libd3q.a $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	d3q.o  libd3q.a	$(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../D3Q/src/d3q.x . )

pw :
	cd ../../ ; make pw

ph :
	cd ../../ ; make ph

pp :
	cd ../../ ; make pp

ld1 :
	cd ../../ ; make ld1

clean :
	- /bin/rm -f  *.o *~ *.d *.mod *.i libd3q.a liblr.a

include make.depend
