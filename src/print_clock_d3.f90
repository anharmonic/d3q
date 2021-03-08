!
! Copyright (C) 2011 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

SUBROUTINE print_clock_d3
  USE io_global,  ONLY : stdout
  IMPLICIT NONE

!   CALL print_clock ('')

  WRITE( stdout, * )
  CALL print_clock ('D3TOTEN')
  !
  WRITE( stdout, * )
  WRITE( stdout,  * ) '     NSCF calculation of ground-state wfcs'
  CALL print_clock('PWSCF')
  CALL print_clock('init_run')
  CALL print_clock('electrons')
  CALL print_clock('c_bands')
  CALL print_clock('cegterg')
  !
  WRITE( stdout, * )
  WRITE( stdout,  * ) '     NSCF calculation of perturbed wfcs'
  CALL print_clock ('generate_dwfc')
  CALL print_clock ('solve_linter')
  CALL print_clock ('ortho')
  CALL print_clock ('cgsolve')
  CALL print_clock ('incdrhoscf')
  CALL print_clock ('set_efsh')
  !
  WRITE( stdout, * )
  CALL print_clock ('ch_psi')
  CALL print_clock ('h_psiq')
  CALL print_clock ('last')

  WRITE( stdout, * )
  CALL print_clock ('firstfft')
  CALL print_clock ('product')
  CALL print_clock ('secondfft')

  WRITE( stdout, * )
  WRITE( stdout,  * ) '     D3 setup and pre-calculations'
  CALL print_clock ('d3_setup')
  CALL print_clock ('d3_setup_noq')
  CALL print_clock ('d3_init')
  CALL print_clock ('gen_dpsi1dv2psi')
  CALL print_clock ('dpsi1dv2psi')

  WRITE( stdout, * )
  WRITE( stdout,  * ) '     D3 matrix elements'
  CALL print_clock('dpsi1dv2dpsi3')
  CALL print_clock('dpsi1dpsi2dv3')
  CALL print_clock('dq1rhodq23v')
  CALL print_clock('rhodq123v')
  CALL print_clock('d3_nlcc_0')
  CALL print_clock('d3_nlcc_123')
  CALL print_clock('d3ionq')
  CALL print_clock('d3_smr_ijk')
  CALL print_clock('d3_smr_ij')
  CALL print_clock('d3_smr_gamma')
  CALL print_clock('d3_exc')
  CALL print_clock('d3_exc_gc')
  CALL print_clock('d3matrix')


  WRITE( stdout, * )
  WRITE( stdout,  * ) '     Input/output'
  CALL print_clock('d3_add_rho_core')
  CALL print_clock('drho_add_phase')
  CALL print_clock('drho_change_q')
  CALL print_clock('read_drho')
  CALL print_clock('davcio')

  WRITE( stdout, * )
  WRITE( stdout,  * ) '     General routines'
  CALL print_clock ('calbec')
  CALL print_clock ('cft3')
  CALL print_clock ('cft3s')
  CALL print_clock ('cinterpolate')
  CALL print_clock ('dq_vscf')

  WRITE( stdout, * )
#ifdef __MPI
  WRITE( stdout,  * ) '     Parallel routines'
  CALL print_clock ('fft_scatter')
#endif
  !
  RETURN
END SUBROUTINE print_clock_d3


SUBROUTINE print_clock_d3_short
  USE io_global,  ONLY : stdout
  IMPLICIT NONE

!   CALL print_clock ('')

  WRITE( stdout, * )
  CALL print_clock('D3TOTEN')
  WRITE( stdout,  '(5x,"--> ",a)' ) 'wfc and dwfc'
  CALL print_clock('nscf')
  CALL print_clock('generate_dwfc')
  WRITE( stdout,  '(5x,"--> ",a)' ) 'matrix elements'
  CALL print_clock('dpsi1dv2dpsi3')
  CALL print_clock('dpsi1dpsi2dv3')
  CALL print_clock('dq1rhodq23v')
  CALL print_clock('rhodq123v')
  CALL print_clock('d3_nlcc_0')
  CALL print_clock('d3_nlcc_123')
  CALL print_clock('d3ionq')
  CALL print_clock('d3_smr_ijk')
  CALL print_clock('d3_smr_ij')
  CALL print_clock('d3_smr_gamma')
  CALL print_clock('d3_exc')
  CALL print_clock('d3matrix')
  WRITE( stdout,  '(5x,"--> ",a)' ) 'total I/O'
  CALL print_clock('davcio')
  !
  RETURN
END SUBROUTINE print_clock_d3_short
