!
! Copyright (C) 20020-2022 Lorenzo Paulatto & Ibrahim G Garba
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE tdph_module
  USE kinds,      ONLY : DP
  USE input_fc,   ONLY : ph_system_info
  USE decompose_d2, ONLY : dynmat_basis, sym_and_star_q
  USE input_fc,           ONLY : forceconst2_grid
  USE iso_c_binding
#include "mpi_thermal.h"
  !
  !INTEGER :: nfar=0

  TYPE tdph_input_type
    !
    CHARACTER(len=256) :: ai = 'md'
    CHARACTER(len=256) :: fmd = 'md.out'
    CHARACTER(len=256) :: ftau, fforce, ftoten
    CHARACTER(len=256) :: file_mat2 = 'mat2R.periodic'
    CHARACTER(len=8) :: fit_type = "force"
    CHARACTER(len=8) :: minimization = "ph+zstar"
    CHARACTER(len=9) :: basis = "mu"
    INTEGER :: nfirst, nskip, nmax, nprint
    REAL(DP) :: e0, thr, T, randomization, alpha_rigid
    !
  END TYPE tdph_input_type
  !
  ! Global variables used by CHI2_LMDIF
  INTEGER :: nq_wedge, nq1, nq2, nq3, nqmax
  TYPE(ph_system_info) :: Si
  TYPE(dynmat_basis),ALLOCATABLE :: dmb(:)
  TYPE(sym_and_star_q),ALLOCATABLE :: symq(:)
  TYPE(tdph_input_type)  :: input
  TYPE(forceconst2_grid) :: fcout
  INTEGER                :: n_steps, nat_sc
  REAL(DP),ALLOCATABLE   :: u(:,:,:), force_harm(:,:,:),  h_energy(:), force_md(:,:,:)
  INTEGER,ALLOCATABLE :: rank(:)
  INTEGER :: zrank
  INTEGER :: ulog
  INTEGER :: mdata
  REAL(dp),ALLOCATABLE :: zstar(:,:,:), zstar_sc(:,:,:), zbasis(:,:,:,:), &
                          force_rgd(:,:,:), tau_sc_alat(:,:)
  COMPLEX(DP),ALLOCATABLE :: rbdyn(:,:,:,:)
  REAL(DP) :: omega_sc, mat_ij(3,3), at_sc(3,3), bg_sc(3,3)
  REAL(DP), PARAMETER :: gamma(3) = (/0._dp,0._dp,0._dp/)

  ! Not necessarily global, but put here for debugging:
  INTEGER ::  mdata_tot, nph
  REAL(kind=C_DOUBLE),ALLOCATABLE,TARGET :: ph_coefficients(:), ph_coefficients0(:), diff_tot(:)
  !
  CONTAINS
  !
  ! read everything from files mat2R and mat3R
  SUBROUTINE READ_INPUT_TDPH(input)
    USE code_input,           ONLY : parse_command_line
    USE cmdline_param_module, ONLY : cmdline_to_namelist
    !
    IMPLICIT NONE
    !
    TYPE(tdph_input_type),INTENT(out) :: input
    CHARACTER(len=256)  :: input_file
    INTEGER :: ios
    !
    ! Input variable, and defaul values:
    CHARACTER(len=256) :: ai  = 'md'
    CHARACTER(len=256) :: fmd = 'md.out'
    CHARACTER(len=256) :: ftau   = 'positions.dat'
    CHARACTER(len=256) :: fforce = 'forces.dat'
    CHARACTER(len=256) :: ftoten = 'pioud.dat'
    CHARACTER(len=256) :: file_mat2 = 'mat2R.periodic'
    CHARACTER(len=8) :: fit_type = "force"
    CHARACTER(len=8) :: minimization = "ph+zstar"
    CHARACTER(len=9) :: basis = "mu"
    INTEGER :: nfirst=1000, nskip=100, nmax=-1, nprint=1000, nread=-1

    INTEGER :: input_unit, aux_unit, err1, err2
    !CHARACTER(len=6), EXTERNAL :: int_to_char
    INTEGER,EXTERNAL :: find_free_unit
    REAL(DP) :: e0 = 0._dp, thr = 1.d-8, T=-1._dp, randomization=0._dp, alpha_rigid = 1._dp
    !
    NAMELIST  / tdphinput / &
        fmd, fforce, ftau, ftoten, &
        ai, file_mat2, fit_type, &
        nfirst, nskip, nmax, nprint, nread, &
        e0, minimization, thr, T, basis, &
        randomization, alpha_rigid

    ioWRITE(*,*) "Waiting for input"
    !
    IF(ionode)THEN
      input_file="input.TDPH"
      CALL parse_command_line(input_file)
      IF(TRIM(input_file)=="-")THEN
        ioWRITE(stdout,'(2x,3a)') "Warning! Reading standard input will probably not work with MPI"
        input_unit = 5
        err1=0
      ELSE
        ioWRITE(stdout,'(2x,3a)') "Reading input file '", TRIM(input_file), "'"
        !input_unit = find_free_unit()
        OPEN(newunit=input_unit, file=input_file, status="OLD", action="READ",iostat=err1)
      ENDIF
      !
      IF(err1==0)THEN
        aux_unit = find_free_unit()
        READ(input_unit, tdphinput, iostat=err1)
        ioWRITE(stdout,'(2x,3a)') "merging with command line arguments"
      ELSE
        !ioWRITE(stdout,'(2x,3a)') "no input file, trying with command line arguments"
        CALL errore("tdph", "bad input file", 1)
      ENDIF
      OPEN(unit=aux_unit, file=TRIM(input_file)//".tmp~", status="UNKNOWN", action="READWRITE")
      CALL cmdline_to_namelist("tdphinput", aux_unit)
      REWIND(aux_unit)
      READ(aux_unit, tdphinput, iostat=err2)
      IF(err2==0)CLOSE(aux_unit, status="DELETE")

      IF(err1/=0 .and. err2/=0 .and. ionode)  WRITE(stdout,'(2x,3a)') "Warning: no input file or parameters!"

      ioWRITE(stdout, tdphinput)
    ENDIF
    !
    CALL bcast_namelist_variables()
    !
    !IF(ANY(<1))  CALL errore("READ_INPUT_TDPH","Invalid nk",1)

    IF(e0==0._dp .and. fit_type=='energy') &
        CALL errore("tdph","need zero energy (e0 = total energy of unit cell) to fit energy difference", 1)

    IF(nmax>0 .and. nread>0) CALL errore("tdph", "You cannot specify both nread and nmax",1)
    IF(nmax<0 .and. nread>0) nmax = nfirst+nskip*nread
    IF(nmax<0 .and. nread<0) nmax = 50000 ! default value when nothing specified
    IF(ANY((/nfirst,nmax,nskip/)<1)) &
        CALL errore("tdph","wrong parameters", 1)
    !
    input%ai            = ai
    input%fmd           = fmd
    input%ftau          = ftau
    input%fforce        = fforce
    input%ftoten        = ftoten
    input%file_mat2     = file_mat2
    input%fit_type      = fit_type
    input%minimization  = minimization
    input%nfirst        = nfirst
    input%nskip         = nskip
    input%nmax          = nmax
    input%nprint        = nprint
    input%e0            = e0
    input%thr           = thr
    input%T             = T
    input%basis         = basis
    input%randomization = randomization
    input%alpha_rigid   = alpha_rigid
    !
    CONTAINS
    SUBROUTINE bcast_namelist_variables()
        USE mpi_thermal, ONLY : mpi_broadcast
        IMPLICIT NONE
        CALL mpi_broadcast(ai)
        CALL mpi_broadcast(fmd)
        CALL mpi_broadcast(ftau)
        CALL mpi_broadcast(fforce)
        CALL mpi_broadcast(ftoten)
        CALL mpi_broadcast(file_mat2)
        CALL mpi_broadcast(fit_type)
        CALL mpi_broadcast(nfirst)
        CALL mpi_broadcast(nskip)
        CALL mpi_broadcast(nmax)
        CALL mpi_broadcast(nprint)
        CALL mpi_broadcast(nread)
        CALL mpi_broadcast(e0)
        CALL mpi_broadcast(minimization)
        CALL mpi_broadcast(thr)
        CALL mpi_broadcast(T)
        CALL mpi_broadcast(basis)
        CALL mpi_broadcast(randomization)
        CALL mpi_broadcast(alpha_rigid)
    END SUBROUTINE
  END SUBROUTINE READ_INPUT_TDPH
  !
  SUBROUTINE set_qe_global_geometry(Si)
    USE input_fc,           ONLY : ph_system_info
    USE cell_base,          ONLY : at, bg, celldm, ibrav, omega
    USE ions_base,          ONLY : nat, ityp, ntyp => nsp, atm, tau, amass
    USE noncollin_module,   ONLY : m_loc, nspin_mag
    !
    IMPLICIT NONE
    TYPE(ph_system_info),INTENT(in) :: Si
    !
    ! Quantum-ESPRESSO symmetry subroutines use global variables
    ! we copy the system data from structure Si
    ntyp   = Si%ntyp
    nat    = Si%nat
    IF(allocated(tau)) DEALLOCATE(tau)
    ALLOCATE(tau(3,nat))
    IF(allocated(ityp)) DEALLOCATE(ityp)
    ALLOCATE(ityp(nat))
    celldm = Si%celldm
    at     = Si%at
    bg     = Si%bg
    omega  = Si%omega
    atm(1:ntyp)    = Si%atm(1:ntyp)
    amass(1:ntyp)  = Si%amass(1:ntyp)
    tau(:,1:nat)   = Si%tau(:,1:nat)
    ityp(1:nat)    = Si%ityp(1:nat)

  END SUBROUTINE
  !
!-----------------------------------------------------------------------
  SUBROUTINE chi_lmdif_c(mdata_tot_, npars_, pars_, diff_tot_, iswitch) BIND(c, name="chi_lmdif_c")
    !-----------------------------------------------------------------------
    ! Calculates the square difference, fdiff2, btw harmonic and ab-initio
    ! forces for n_steps molecur dyanmics simulation
    !
    !USE tdph_module,  ONLY : nfar
    USE iso_c_binding
    USE mpi_thermal,     ONLY : mpi_bsum, allgather_vec, my_id
    USE decompose_d2,    ONLY : recompose_fc
    USE timers,          ONLY : t_chi2, t_comm, t_zstar, t_rigid
    USE harmonic_module, ONLY : harmonic_force_md
    USE input_fc,        ONLY : write_fc2
    ! rigid block model:
    USE decompose_zstar, ONLY : recompose_zstar, zstar_to_supercell
    USE rigid_d3,           ONLY : rgd_blk_d3
    USE asr2_module,        ONLY : impose_asr2

    IMPLICIT NONE
    INTEGER(kind=C_INT),INTENT(in)    :: mdata_tot_, npars_
    REAL(kind=C_DOUBLE),INTENT(in)    :: pars_(npars_)
    REAL(kind=C_DOUBLE),INTENT(inout) :: diff_tot_(mdata_tot_)
    INTEGER(kind=C_INT),INTENT(inout) :: iswitch

    REAL(kind=c_DOUBLE) :: diff(mdata)
    INTEGER,SAVE :: iter = 0
    CHARACTER (LEN=6),  EXTERNAL :: int_to_char
    REAL(DP) :: chi2
    INTEGER :: istep,i,j
    !REAL(DP),ALLOCATABLE :: ph_aux(:), zstar_aux(:)
    REAL(DP),ALLOCATABLE,SAVE :: ph_aux(:), zstar_aux(:), zstar_old(:)
    LOGICAL,SAVE :: done_zstar=.false.

    CALL t_chi2%start()

    IF(.not.allocated(ph_aux)) &
          ALLOCATE(ph_aux(nph),zstar_aux(zrank),zstar_old(zrank))
    IF(npars_ == nph)THEN
      ph_aux    = pars_(1:nph)
      zstar_aux = ph_coefficients(nph+1:nph+zrank)
    ELSEIF (npars_ == zrank) THEN
      done_zstar=.false.
      ph_aux    = ph_coefficients(1:nph)
      zstar_aux = pars_(1:zrank)
    ELSEIF (npars_ == nph+zrank) THEN
      done_zstar=.false.
      ph_aux    = pars_(1:nph)
      zstar_aux = pars_(nph+1:nph+zrank)
    ELSE
      CALL errore('minim','what',1)
    ENDIF

    CALL recompose_fc(Si, nq_wedge, symq, dmb, rank, nph, ph_aux,&
                      nq1, nq2, nq3, nqmax, 0, fcout)

    CALL harmonic_force_md(n_steps, nat_sc, Si, fcout, u, force_harm, h_energy)

    ! Computing long-range forces is quite slow, we try to do it only when it is really necessary
    IF(Si%lrigid)THEN ! I put this IF first because the others are always false when this one is
    IF(.not.done_zstar.AND.ANY(zstar_old /= zstar_aux)) THEN
      !
      !print*, "recom. zstar", zstar_aux
      zstar_old = zstar_aux
      !
      CALL t_zstar%start()
      done_zstar = .true.
      CALL recompose_zstar(Si%nat, zrank, zbasis, zstar_aux, zstar)
      CALL zstar_to_supercell(Si%nat, nat_sc, zstar, zstar_sc)
      rbdyn = 0._dp
      CALL t_rigid%start()
      CALL rgd_blk_d3(2,2,2, nat_sc, rbdyn, gamma, tau_sc_alat, Si%epsil, zstar_sc, bg_sc, &
                   omega_sc, Si%alat, .false., +1._dp) !, alpha=input%alpha_rigid)
      CALL t_rigid%stop()
!$OMP PARALLELDO DEFAULT(shared) PRIVATE(istep,i,j,mat_ij)
      DO istep = 1, n_steps
        DO i = 1, nat_sc
          force_rgd(:,i,istep) = 0._dp
          DO j = 1, nat_sc
            ! Dynamical matrix at Gamma should be real, let's remove any noise
            mat_ij = DBLE(rbdyn(:,:,i,j))
            force_rgd(:,i,istep) =  force_rgd(:,i,istep) - MATMUL(mat_ij,u(:,j,istep))
          END DO
        END DO
      END DO
!$OMP END PARALLELDO
      CALL t_zstar%stop()
    ENDIF
    !
    ! Always include rigid contribution, even if it was not computed on this call
    force_harm =  force_harm + force_rgd
    ENDIF ! Si%lrigid
    !
    SELECT CASE(input%fit_type)
    CASE('force', 'forces')
      force_harm(1:3,1:nat_sc,1:n_steps) = force_harm(1:3,1:nat_sc,1:n_steps)&
                                         - force_md(1:3,1:nat_sc,1:n_steps)
      diff = RESHAPE( force_harm(1:3,1:nat_sc,1:n_steps), (/ mdata /) )
    CASE('energy')
      ! not implemented with MPI_GATHER
      STOP 100
      !DO i = 1, n_steps
      !diff = diff + (h_energy(i)-toten_md(i))**2
      !ENDDO
    CASE('thforce')
      ! not implemented with MPI_GATHER
      STOP 100
      !diff = diff + RESHAPE( (ABS(force_harm(:,:,i) - force_md(:,:,i))**2 *EXP(-toten_md(i)) ), (/ mfc /) )
    CASE DEFAULT
      CALL errore("tdph", 'unknown chi2 method', 1)
    END SELECT

    CALL t_comm%start()
    CALL allgather_vec(mdata, diff, diff_tot_)
    CALL t_comm%stop()

    ! lmdif_c uses 0 to signal new steps, 1 at the first evaluation
    IF(iswitch==0 .or. iswitch==1)THEN
      iter = iter+1
      chi2 = SUM(diff**2)
      CALL mpi_bsum(chi2)
      chi2=SQRT(chi2)/mdata_tot_
      ioWRITE(*,'(i10,e15.7)') iter, chi2
      IF(isnan(chi2))THEN
        WRITE(10000+my_id, *)  diff
        STOP 777
      ENDIF

      ioWRITE(ulog, "(2i10,e17.6,9999e15.6)") iter, iswitch, chi2, pars_
      ! Every input%nprint steps have a look
      IF(MODULO(iter,input%nprint)==0) &
        CALL write_fc2("matOUT.iter_"//TRIM(int_to_char(iter)), Si, fcout)
    ENDIF

    CALL t_chi2%stop()
    !
  END SUBROUTINE chi_lmdif_c
    !
!  void fcn(void * farg, int m, int n, const double *x, double *fvec)
!-----------------------------------------------------------------------
  SUBROUTINE chi_lmdif_c_para(farg, mdata_tot_, npars_, pars_, diff_tot_, iswitch_) &
    BIND(c, name="chi_lmdif_c_para")
    !-----------------------------------------------------------------------
    ! Calculates the square difference, fdiff2, btw harmonic and ab-initio
    ! forces for n_steps molecur dyanmics simulation
    USE iso_c_binding
    IMPLICIT NONE
    TYPE(c_ptr),INTENT(in) :: farg
    INTEGER(kind=C_INT),INTENT(in),VALUE    :: mdata_tot_, npars_, iswitch_
    REAL(kind=C_DOUBLE),INTENT(in)    :: pars_(npars_)
    REAL(kind=C_DOUBLE),INTENT(inout) :: diff_tot_(mdata_tot_)
    !
    INTEGER(kind=C_INT) :: iswitch_aux
    !
    !print*, "wrap>", mdata_tot_, npars_
    !print*, "wrap2>", pars_
    iswitch_aux = iswitch_
    CALL chi_lmdif_c(mdata_tot_, npars_, pars_, diff_tot_, iswitch_aux)
    !print*, ">>>>>", SUM(diff_tot_**2)
    !
  END SUBROUTINE
END MODULE

!
!----------------------------------------------------------------------------
PROGRAM tdph
  !----------------------------------------------------------------------------
  !
  USE kinds,              ONLY : DP
  USE constants,          ONLY : amu_ry, K_BOLTZMANN_SI, K_BOLTZMANN_RY, RYTOEV !ry_to_kelvin
  USE parameters,         ONLY : ntypx
  !USE mp,                 ONLY : mp_bcast
  !USE mp_global,          ONLY : mp_startup, mp_global_end
  !USE mp_world,           ONLY : world_comm
  !USE environment,        ONLY : environment_start, environment_end
  ! symmetry
  USE symm_base,          ONLY : s, invs, nsym, find_sym, set_sym_bl, &
                                  irt, copy_sym, nrot, inverse_s, t_rev
  USE noncollin_module,   ONLY : m_loc, nspin_mag
  USE lr_symm_base,       ONLY : rtau, nsymq, minus_q, irotmq, gi, gimq, invsymq
  USE control_lr,         ONLY : lgamma
  USE decompose_d2,       ONLY : smallg_q_fullmq, find_d2_symm_base, sym_and_star_q, dotprodmat, &
                                 make_qstar_d2, allocate_sym_and_star_q, tr_star_q, dynmat_basis
  USE cmdline_param_module
  USE input_fc,           ONLY : forceconst2_grid, ph_system_info, read_system, aux_system, read_fc2, &
                                 div_mass_fc2, multiply_mass_dyn, write_fc2
  USE fc2_interpolate,    ONLY : fftinterp_mat2, mat2_diag, dyn_cart2pat
  USE asr2_module,        ONLY : impose_asr2
  USE quter_module,       ONLY : quter
  USE mpi_thermal,        ONLY : start_mpi, stop_mpi, num_procs, mpi_broadcast, my_id
  USE tdph_module
  ! harmonic
  USE read_md_module,     ONLY : read_md, read_pioud, fc_to_supercell, read_max_steps_para
  USE harmonic_module,    ONLY : harmonic_force_md
  ! minipack
  USE lmdif_module,       ONLY : lmdif0
#if defined  (__SCALAPACK)
  USE lmdif_p_module,     ONLY : lmdif_p0, lmdif_c0, plmdif_c0
#else
  USE lmdif_p_module,     ONLY : lmdif_p0, lmdif_c0
#endif
  USE timers
  USE random_numbers,     ONLY : randy
  USE decompose_d2,       ONLY : recompose_fc
  USE rigid_d3,           ONLY : rgd_blk_d3
  USE iso_c_binding
  USE decompose_zstar,    ONLY : find_zstar_symm_base, dotprodzstar, recompose_zstar
  !
  IMPLICIT NONE
  !
  CHARACTER(len=7),PARAMETER :: CODE="TDPH"
  CHARACTER(len=256) :: fildyn, filout
  INTEGER :: ierr, nargs
  !
  INTEGER       :: nqq, nq_done

  REAL(DP),ALLOCATABLE      :: x_q(:,:), w_q(:)
  ! for harmonic_module
  REAL(DP),ALLOCATABLE      :: tau_md(:,:,:), tau_sc(:,:), &
                               force_diff(:), toten_md(:)
  ! for lmdf1
  INTEGER                   :: n_steps_tot, j_steps
  !
  REAL(DP) :: syq(3,48), force_ratio(3), aux
  LOGICAL :: sym(48), skip_equivalence, time_reversal !, lrigid_save
  !
  COMPLEX(DP),ALLOCATABLE :: phi(:,:,:,:), d2(:,:), w2(:,:), &
                             star_wdyn(:,:,:,:, :), star_dyn(:,:,:)
  REAL(DP),ALLOCATABLE :: decomposition(:), xqmax(:,:), &
                          metric(:)
  INTEGER :: i,j, istep, icar,jcar, na,nb, iq, iph, iswitch, first_step, n_skip
#if defined  (__SCALAPACK)
  INTEGER,POINTER :: ptrdummy => null()
#endif

  TYPE(forceconst2_grid) :: fc
  !
  CALL start_mpi()
  !
  ! Read namelist tdphinput
  CALL READ_INPUT_TDPH(input)

  CALL read_fc2(input%file_mat2, Si, fc)
  CALL impose_asr2("simple", Si%nat, fc, Si%zeu)
  !lrigid_save = Si%lrigid
  !Si%lrigid = .false.
  CALL aux_system(Si)
  CALL div_mass_fc2(Si, fc)
  CALL fc_to_supercell(Si, fc, at_sc, bg_sc, omega_sc, nat_sc, tau_sc, zstar_sc)
  !
  CALL set_qe_global_geometry(Si)
  !
  ! ######################### symmetry setup #########################
  ! Symmetry setup uses global variable, which is not ideal, but
  ! not a huge problem since we only need to do this part once
  ! at the beginning.
  ! ~~~~~~~~ setup bravais lattice symmetry ~~~~~~~~
  CALL set_sym_bl ( )
  ioWRITE(stdout, '(5x,a,i3)') "Symmetries of bravais lattice: ", nrot
  !
  ! ~~~~~~~~ setup crystal symmetry ~~~~~~~~
  IF(.not.allocated(m_loc))  THEN
    ALLOCATE(m_loc(3,Si%nat))
    m_loc = 0._dp
  ENDIF

  CALL find_sym ( Si%nat, Si%tau, Si%ityp, .false., m_loc )
  ioWRITE(stdout, '(5x,a,i3)') "Symmetries of crystal:         ", nsym
  !
  ! Find the reduced grid of q-points:
  skip_equivalence = .FALSE.
  time_reversal    = .TRUE.
  nq1 = fc%nq(1)
  nq2 = fc%nq(2)
  nq3 = fc%nq(3)
  nqmax = nq1*nq2*nq3
  ALLOCATE(x_q(3,nqmax), w_q(nqmax))
  CALL kpoint_grid( nsym, time_reversal, skip_equivalence, s, t_rev, Si%bg, nqmax,&
                    0,0,0, nq1,nq2,nq3, nq_wedge, x_q, w_q )
  !
  ioWRITE(stdout, *) "Generated ", nq_wedge, "points"

  ALLOCATE(rtau( 3, 48, Si%nat), d2(3*Si%nat,3*Si%nat))

  ! Variable to hold the dyn matrix and q-points of the entire grid
  nq_done = 0
  ALLOCATE(star_wdyn(3,3,Si%nat,Si%nat, nqmax))
  ALLOCATE(xqmax(3,nqmax))

  ! For every q-point in the irreducible wedge, we find its symmetry
  ! and the basis of the space of symmetry-constrained dynfactoramical matrices
  ! Again, this part uses gloabl variables which is a annoying, but once
  ! we have the set of basis matrices, we don't have to touch it
  ! anymore.
  ALLOCATE(dmb(nq_wedge))
  ALLOCATE(rank(nq_wedge))
  ALLOCATE(symq(nq_wedge))

  CALL t_init%start()

  ioWRITE(stdout, '("=====================")')
  ! First loop, find symmetry of every q-point
  Q_POINTS_LOOP_a : &
  DO iq = 1, nq_wedge
    ! ioWRITE(stdout, *) "____[[[[[[[", iq, "]]]]]]]]____"
    ! ioWRITE(stdout, '(i6, 3f12.4)') iq, x_q(:,iq)
    !
    ! ~~~~~~~~ setup small group of q symmetry ~~~~~~~~
    ! part 1: call smallg_q and the copy_sym,
    symq(iq)%xq = x_q(:,iq)
    minus_q = .true.

    sym = .false.
    sym(1:nsym) = .true.
    CALL smallg_q_fullmq(symq(iq)%xq, 0, Si%at, Si%bg, nsym, s, sym, minus_q)
    nsymq = copy_sym(nsym, sym)
    ! recompute the inverses as the order of sym.ops. has changed
    CALL inverse_s ( )

    ! part 2: this computes gi, gimq
    call set_giq (symq(iq)%xq,s,nsymq,nsym,irotmq,minus_q,gi,gimq)
!    WRITE(stdout, '(5x,a,i3)') "Symmetries of small group of q:", nsymq
!    IF(minus_q) WRITE(stdout, '(10x,a)') "in addition sym. q -> -q+G"
    !
    ! finally this does some of the above again and also computes rtau...
    CALL sgam_lr(Si%at, Si%bg, nsym, s, irt, Si%tau, rtau, Si%nat)
    !
    ! Now, I copy all the symmetry definitions to a derived type
    ! which I'm going to use EXCLUSIVELY from here on
    CALL allocate_sym_and_star_q(Si%nat, symq(iq))
    symq(iq)%nrot = nrot
    symq(iq)%nsym  = nsym
    symq(iq)%nsymq = nsymq
    symq(iq)%minus_q = minus_q
    symq(iq)%irotmq = irotmq
    symq(iq)%s = s
    symq(iq)%invs = invs
    symq(iq)%rtau = rtau
    symq(iq)%irt= irt
  ENDDO Q_POINTS_LOOP_a
  ioWRITE(stdout, '(5x,a,999i5)') "Number of sym.ops. of each q-point:", symq(:)%nsymq

  ioWRITE(stdout, '("=====================")')
  ! Find symmetric basis for each point (MPI parallelisation)
  Q_POINTS_LOOP_b1 : &
  DO iq = 1, nq_wedge
    ! Simple parallelization of basis search
    IF(MOD(iq,num_procs)==my_id)THEN
        ! the next subroutine uses symmetry from global variables to find the basis of crystal-symmetric
        ! matrices at this q point
        !IF(input%basis=='mu')THEN
          CALL fftinterp_mat2(symq(iq)%xq, Si, fc, d2, gamma)
          d2 = multiply_mass_dyn(Si,d2)
        !ENDIF
        CALL find_d2_symm_base(symq(iq)%xq, rank(iq), dmb(iq)%basis, &
         Si%nat, Si%at, Si%bg, symq(iq)%nsymq, symq(iq)%minus_q, &
         symq(iq)%irotmq, symq(iq)%rtau, symq(iq)%irt, symq(iq)%s, symq(iq)%invs, &
         d2, input%basis)
    ENDIF
  ENDDO Q_POINTS_LOOP_b1
  !
  ! Distribute basis via MPI
  CALL t_comm%start()
  Q_POINTS_LOOP_b2 : &
  DO iq = 1, nq_wedge
    CALL mpi_broadcast(rank(iq), root=MOD(iq,num_procs))
    IF(MOD(iq,num_procs)/=my_id) ALLOCATE(dmb(iq)%basis(Si%nat3, Si%nat3, rank(iq)))
    CALL mpi_broadcast(Si%nat3, Si%nat3, rank(iq), dmb(iq)%basis, root=MOD(iq,num_procs))
  ENDDO Q_POINTS_LOOP_b2
  CALL t_comm%stop()
  !
  ioWRITE(stdout, '("=====================")')
  ! Find star of q points
  Q_POINTS_LOOP_c : &
  DO iq = 1, nq_wedge
      !
    ! Calculate the list of points making up the star of q and of -q
    CALL tr_star_q(symq(iq)%xq, Si%at, Si%bg, symq(iq)%nsym, symq(iq)%s, symq(iq)%invs, &
                   symq(iq)%nq_star, symq(iq)%nq_trstar, symq(iq)%sxq, &
                   symq(iq)%isq, symq(iq)%imq, .false. )

    !ioWRITE(stdout, '(5x,a,2i5)') "Found star of q and -q", symq(iq)%nq_star, symq(iq)%nq_trstar
    syq = symq(iq)%sxq
    call cryst_to_cart(symq(iq)%nq_trstar, syq, Si%at, -1)
    DO i = 1, symq(iq)%nq_trstar
       syq(1,i) = MODULO(syq(1,i), 1._dp)
       syq(2,i) = MODULO(syq(2,i), 1._dp)
       syq(3,i) = MODULO(syq(3,i), 1._dp)
       syq(:,i) = syq(:,i) * (/nq1, nq2, nq3/)
       !ioWRITE(stdout,'(i4,3i3,l2)') i, NINT(syq(:,i)), (i>symq(iq)%nq_star)
    ENDDO

  ENDDO Q_POINTS_LOOP_c
  ioWRITE(stdout, '(5x,a,999i5)') "Points in the star of each q-point:", symq(:)%nq_trstar
  !
  ! Find symmetric basis for effetcive charges
  IF(Si%lrigid) THEN
    ! I use the symmetry of Gamma as it is the same as the one of the crystal
    IF(ANY(symq(1)%xq/=0._dp)) CALL errore("star","gamma should be the first point",1)
    CALL find_zstar_symm_base(zrank, zbasis, Si%nat, Si%at, Si%bg, symq(1)%nsym, symq(1)%irt, &
                              symq(1)%s, Si%zeu, "simple" )
  ELSE
    zrank = 0
  ENDIF
  !
  ! Number of degrees of freedom for the entire grid:
  nph = SUM(rank)!+zrank
  ioWRITE(stdout, '(5x,a,2i5)') "TOTAL number of degrees of freedom", nph

  ! Allocate a vector to hold the decomposed phonons over the entire grid
  ! I need single vector in order to do minimization, otherwise a derived
  ! type would be more handy
  ALLOCATE(ph_coefficients(nph+zrank), ph_coefficients0(nph+zrank))
  iph = 0
  Q_POINTS_LOOP2 : &
  DO iq = 1, nq_wedge
    !
    ! Interpolate the system dynamical matrix at this q
    CALL fftinterp_mat2(symq(iq)%xq, Si, fc, d2, gamma)
    ! Remove the mass factor, I cannot remove it before because the effective
    ! charges/long range interaction code assumes it is there
    d2 = multiply_mass_dyn(Si,d2)
    !
    ! Decompose the dynamical matrix over the symmetric basis at this q-point
    ioWRITE(stdout,'(2x,a,3f12.4)') "Phonon parameters xq=",symq(iq)%xq
    DO i = 1,rank(iq)
      iph = iph +1
      ph_coefficients(iph) = dotprodmat(3*Si%nat,d2, dmb(iq)%basis(:,:,i))
      ioWRITE(stdout,"(i3,1f12.6)") i, ph_coefficients(iph)
    ENDDO
    !
  ENDDO Q_POINTS_LOOP2
  !
  ! Decompose effective charges
  IF(Si%lrigid) THEN
    ioWRITE(stdout,'(2x,a)') "Effective charges parameters"
    DO i = 1,zrank
      iph = iph +1
      ph_coefficients(iph) = dotprodzstar(Si%nat,Si%zeu, zbasis(:,:,:,i))
      ioWRITE(stdout,"(i3,1f12.6)") i, ph_coefficients(iph)
    ENDDO
  ENDIF
  !
  ph_coefficients0 = ph_coefficients
!-----------------------------------------------------------------------
  ! Variables that can be adjusted according to need ...
  !
  !n_steps    = input%nmax   ! total molecular dynamics steps TO READ
  n_steps    = read_max_steps_para(input%nmax)   ! total molecular dynamics steps TO READ on this CPU
  first_step = input%nfirst ! start reading from this step
  n_skip     = input%nskip  ! number of steps to skip

  CALL t_init%stop()
!###################  end of initialization ####################################################

  ALLOCATE(tau_md(3,nat_sc,n_steps))
  ALLOCATE(force_md(3,nat_sc,n_steps))
  ALLOCATE(toten_md(n_steps))
  CALL t_read%start()
  IF(input%ai=="md") THEN
  CALL read_md(input%fmd, input%e0, nat_sc, Si%alat, at_sc, first_step, n_skip, n_steps, &
               n_steps_tot, tau_md, force_md, toten_md, tau_sc, u)
  ELSE IF(input%ai=="pioud")THEN
    CALL read_pioud(input%ftau, input%fforce, input%ftoten, input%e0, nat_sc, Si%alat, &
                    first_step, n_skip, n_steps, n_steps_tot, tau_md, force_md, toten_md, tau_sc, u)
  ELSE
    CALL errore("tdph","unknown input format", 1)
  ENDIF
  IF(n_steps_tot > input%nmax) CALL errore('tdph','read more steps than expected',1)
  CALL t_read%stop()

!###################  end of data input ####################################################
  ! Compute force from rigid block model (i.e. from long range effective charges interaction)
  ! and remove them from MD forces to exclude long range effect from fit
  IF(Si%lrigid)THEN
    !
    ALLOCATE(force_rgd(3,nat_sc,n_steps))
    ALLOCATE(rbdyn(3,3,nat_sc,nat_sc))
    ALLOCATE(tau_sc_alat(3,nat_sc))
    ALLOCATE(zstar(3,3,Si%nat)) ! will be used during minimization
    tau_sc_alat = tau_sc/Si%alat
  END IF ! lrigid

!###################  end of rigid block ####################################################

  mdata = 3*nat_sc*n_steps
  mdata_tot = 3*nat_sc*n_steps_tot
  ALLOCATE(force_diff(3*nat_sc))
  ALLOCATE(diff_tot(mdata_tot))
  ALLOCATE(force_harm(3,nat_sc,n_steps))

  IF(input%randomization/=0._dp)THEN
    IF(ionode)THEN
      !
      aux = 2*SUM(ABS(ph_coefficients))/nph*ABS(input%randomization)
      !
      IF(input%randomization>0._dp)THEN
        WRITE(*,"(x,a,f12.6)") "Adding random noise to initial phonons up to ", aux
        DO i = 1, nph
          ph_coefficients(i) = ph_coefficients(i) + 2*aux*(randy()-.5_dp)
        ENDDO
      ELSE IF(input%randomization<0._dp)THEN
        WRITE(*,"(x,a,f12.6)") "Randomizing initial phonons up to ", aux
        DO i = 1, nph
          ph_coefficients(i) = aux*2*(randy()-.5_dp)
        ENDDO
      ENDIF
      !
    ENDIF
    ! be sure that every CPU has the same coefficients
    CALL mpi_broadcast(nph, ph_coefficients)
  ENDIF

  OPEN(newunit=ulog, file="tdph.log", form="formatted", status="unknown")

  CALL t_minim%start()
  ioWRITE(*,*) "Starting minimization: ", TRIM(input%minimization)
  SELECT CASE (input%minimization)
  CASE("lmdif")
    CALL lmdif_p0(chi_lmdif_c, mdata_tot, nph, ph_coefficients(1:nph), diff_tot, input%thr, iswitch)
  CASE("ph+zstar")
    CALL lmdif_c0(chi_lmdif_c, mdata_tot, nph, ph_coefficients(1:nph), diff_tot, input%thr, iswitch)
     IF(Si%lrigid)THEN
       ioWRITE(*,*) "Minimizing z^star "
       ph_coefficients(nph+1:nph+zrank) = ph_coefficients(nph+1:nph+zrank)*0.9
       CALL lmdif_c0(chi_lmdif_c, mdata_tot, zrank, ph_coefficients(nph+1:nph+zrank), diff_tot, input%thr, iswitch)
     ENDIF
#if defined  (__SCALAPACK)
  CASE("para")
      CALL plmdif_c0(chi_lmdif_c_para, ptrdummy, mdata_tot, nph, ph_coefficients(1:nph), diff_tot, input%thr, iswitch)
#endif
  CASE("ph")
    CALL lmdif_c0(chi_lmdif_c, mdata_tot, nph, ph_coefficients(1:nph), diff_tot, input%thr, iswitch)
  CASE("global")
    CALL lmdif_c0(chi_lmdif_c, mdata_tot, nph+zrank, ph_coefficients, diff_tot, input%thr, iswitch)
  CASE("none")
    ! Do one estimation of the forces to have something to print out
    iswitch = 1
    CALL chi_lmdif_c(mdata_tot, nph+zrank, ph_coefficients, diff_tot, iswitch)
  CASE DEFAULT
    CALL errore("tdph","unknown minimization engine requested",1)
  END SELECT
  CALL t_minim%stop()
  !
!###################  end of minimization ####################################################
! Write to file the final matrices in "periodic" form, the final fitted forces, etc
  CALL recompose_fc(Si, nq_wedge, symq, dmb, rank, nph, ph_coefficients(1:nph),&
                    nq1, nq2, nq3, nqmax, 0, fcout)
  !Si%lrigid = lrigid_save
  IF(Si%lrigid) CALL recompose_zstar(Si%nat, zrank, zbasis, ph_coefficients(nph+1:nph+zrank), Si%zeu)
  !CALL impose_asr2("simple", Si%nat, fcout, Si%zeu)
  !Si%lrigid = lrigid_save
  CALL write_fc2("matOUT.periodic", Si, fcout)
  !Si%lrigid = .false.

  ! Force ratio
  OPEN(116,file="force_ratio.dat",status="unknown")
  !
  DO i = 1, n_steps
   WRITE(116,*) "i_step = ",i
   DO j = 1, nat_sc
     WHERE(force_md(:,j,i)/=0._dp)
       force_ratio = (force_harm(:,j,i))/force_md(:,j,i)
     ELSEWHERE
       force_ratio = 0._dp
     END WHERE

     ioWRITE(116,'(1f14.6,5x,3(3f14.6,5x))') DSQRT(norm2(force_ratio)),&
                             force_ratio, force_harm(:,j,i), force_md(:,j,i)
     END DO
  END DO
  CLOSE(116)
  !
  ! Harmonic energy
  !
!  OPEN(117,file="h_enr.dat",status="unknown")
!  DO i = 1, n_steps
!  !WRITE(117,*) "i_step = ",i
!    ioWRITE(117,'(i5,3(E13.6, 3x))') i, h_energy(i), toten_md(i), &
!                                         EXP(-toten_md(i)/(K_BOLTZMANN_RY*300.0_DP))
!  END DO
  !
  ! Write to file the matrices in "centered" form for later Fourier interpolation
  CALL recompose_fc(Si, nq_wedge, symq, dmb, rank, nph, ph_coefficients(1:nph),&
                    nq1, nq2, nq3, nqmax, 2, fcout)
  IF(Si%lrigid) CALL recompose_zstar(Si%nat, zrank, zbasis, ph_coefficients(nph+1:nph+zrank), Si%zeu)
  !CALL impose_asr2("simple", Si%nat, fcout, Si%zeu)
  !Si%lrigid = lrigid_save
  CALL write_fc2("matOUT.centered", Si, fcout)

  CLOSE(ulog)
  ioWRITE(stdout,'("   * WALL : ",f12.4," s")') get_wall()
  CALL print_timers_header()
  CALL t_read%print()
  CALL t_init%print()
  CALL t_minim%print()
  ioWRITE(stdout,'("   * inside minimization ")')
  CALL t_chi2%print()
  ioWRITE(stdout,'("   * inside chi2 ")')
  CALL t_comm%print()
  CALL t_force%print()
  CALL t_recom%print()
  IF(Si%lrigid) THEN
    CALL t_zstar%print()
    ioWRITE(stdout,'("   * inside zstar ")')
    CALL t_rigid%print()
  ENDIF

  CALL stop_mpi()

 ! ---- the program is over ---- !
 !
 CONTAINS
  !----------------------------------------------------------------------------
 END PROGRAM tdph
!------------------------------------------------------------------------------
