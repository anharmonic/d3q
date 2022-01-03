set(src_d3q
    src/allocate_d3.f90
    src/cmdline_param.f90
    src/d2mxc.f90
    src/d3_basis.f90
    src/d3com.f90
    src/d3_debug.f90
    src/d3_exc.f90
    src/d3_exc_gc.f90
    src/d3_grid.f90
    src/d3_h_psi.f90
    src/d3_init.f90
    src/d3_iofiles.f90
    src/d3ionq.f90
    src/d3matrix.f90
    src/d3matrix_io2.f90
    src/d3_nlcc.f90
    src/d3_open.f90
    src/d3_readin.f90
    src/d3_recover.f90
    src/d3_reset.f90
    src/d3_restart.f90
    src/d3_setup.f90
    src/d3_shuffle.f90
    src/d3_symmetry.f90
    src/d3_valence.f90
    src/dpsi1dpsi2dv3.f90
    src/dpsi1dv2dpsi3.f90
    src/dq1rhodq23v.f90
    src/dq_vscf.f90
    src/drho_add_phase.f90
    src/dvdpsi.f90
    src/efermi_shift.f90
    src/gen_dpdvp.f90
    src/h_psiq.f90
    src/incdrhoscf2.f90
    src/print_clock_d3.f90
    src/qstar_d3.f90
    src/rhodq123v.f90
    src/rotate_d3.f90
    src/run_nscf_d3.f90
    src/set_kplus3q.f90
    src/setup_d3gc.f90
    src/solve_linter_d3q.f90
    src/stop_d3.f90
    src/w_1gauss.f90
    src/write_d3dyn.f90)

qe_add_library(qe_d3q ${src_d3q})
target_link_libraries(
    qe_d3q
    PRIVATE qe_modules
            qe_lr_modules
            qe_phonon_ph
            qe_pw
            qe_mpi_fortran
            qe_upflib
	    qe_xclib
            qe_fftx)

###########################################################
# d3q.x
###########################################################
set(src_d3q_x src/d3q.f90)
qe_add_executable(qe_d3q_exe ${src_d3q_x})
set_target_properties(qe_d3q_exe PROPERTIES OUTPUT_NAME d3q.x)
target_link_libraries(qe_d3q_exe 
	PRIVATE 
	qe_modules
        qe_lr_modules
	qe_phonon_ph
	qe_pw
	qe_lax
        qe_xclib
	qe_fftx
	qe_d3q)

qe_install_targets(
    # Libraries
    qe_d3q
    # Executables
    qe_d3q_exe)



set (src_thermal2
     thermal2/ACRS.f90
     thermal2/asr2.f90
     thermal2/casimir.f90
     thermal2/code_input.f90
     thermal2/decompose_d2.f90
     thermal2/dynbubble.f90
     thermal2/eos.f90
     thermal2/f3_bwfft.f90
     thermal2/fc2_interp.f90
     thermal2/fc3_interp.f90
     thermal2/final_state.f90
     thermal2/functions.f90
     thermal2/harmonic.f90
     thermal2/input_fc.f90
     thermal2/isotopes.f90
     thermal2/linewidth.f90
     thermal2/mc_grids.f90
     thermal2/merge_degenerate.f90
     thermal2/more_constants.f90
     thermal2/mpi_thermal.f90
     thermal2/nanoclock.f90
     thermal2/nanosec.c
     thermal2/neutrons.f90
     thermal2/nist_isotopes_db.f90
     thermal2/overlap.f90
     thermal2/ph_dos.f90
     thermal2/ph_system.f90
     thermal2/ph_velocity.f90
     thermal2/posix_signal.f90
     thermal2/q_grids.f90
     thermal2/quter.f90
     thermal2/read_md.f90
     thermal2/timers.f90
     thermal2/variational_tk.f90)

qe_add_library(qe_thermal2 ${src_thermal2})
target_link_libraries(
    qe_thermal2
    PRIVATE qe_modules
            qe_lr_modules
            qe_phonon_ph
            qe_pw
            qe_mpi_fortran
            qe_upflib
            qe_xclib
            qe_fftx
	    qe_d3q)

#set_target_properties(qe_thermal2)

#library
qe_install_targets( qe_thermal2 )

##########################################################
# thermal2 list of binaries:
##########################################################

set(qe_thermal2_exe
	lw
        tk
        qha
        asr3
        sqom
        sparse
        q2r
        r2q
        qq2rr
        import_shengbte
        recenter
        interpolate2
        tdph
	db
        rmzeu)

add_custom_target(qe_d3_thermal2
    COMMENT
        "phonon-phonon interactions")

foreach(exe IN ITEMS ${qe_thermal2_exe})
  add_dependencies(qe_d3_thermal2 qe_d3_${exe}_exe)

  set(d3_${exe}_x thermal2/PROGRAM_${exe}.f90)
  qe_add_executable(qe_d3_${exe}_exe ${d3_${exe}_x})
  set_target_properties(qe_d3_${exe}_exe PROPERTIES OUTPUT_NAME d3_${exe}.x)
  target_link_libraries(qe_d3_${exe}_exe  
	PRIVATE
       	qe_thermal2
	qe_modules
        qe_lr_modules
	qe_phonon_ph
	qe_pw
	qe_lax
	qe_fftx
	qe_d3q)

  #library
  qe_install_targets(qe_d3_${exe}_exe)
endforeach()

##########################################################
add_custom_target(d3q
    DEPENDS
        qe_d3q_exe
        qe_d3_thermal2
    COMMENT
        "D3Q code to compute third derivative of total energy, and tools to compute phonon linewidth and thermal conductivity")
