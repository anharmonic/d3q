!
! Written by Lorenzo Paulatto (2014-2016) IMPMC @ UPMC / CNRS UMR7590
!  Dual licenced under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!  and under the GPLv2 licence and following, see
!  <http://www.gnu.org/copyleft/gpl.txt>
!
MODULE r2q_program

  USE timers

#include "mpi_thermal.h"
  CONTAINS

  SUBROUTINE JOINT_DOS(input, out_grid, S, fc)
    USE code_input,       ONLY : code_input_type
    USE kinds,            ONLY : DP
    USE input_fc,         ONLY : forceconst2_grid, ph_system_info
    USE q_grids,          ONLY : q_grid, q_basis, setup_grid, prepare_q_basis
    USE constants,        ONLY : RY_TO_CMM1, pi
    USE functions,        ONLY : f_bose, f_gauss
    USE fc2_interpolate,  ONLY : freq_phq, bose_phq
    USE mpi_thermal,      ONLY : mpi_bsum, start_mpi, stop_mpi
    USE random_numbers,   ONLY : randy
    USE mpi_thermal,      ONLY :  mpi_bsum
    USE more_constants,   ONLY : write_conf
    IMPLICIT NONE
    TYPE(code_input_type)    :: input
    TYPE(q_grid), INTENT(in) :: out_grid
    TYPE(ph_system_info)     :: S
    TYPE(forceconst2_grid),INTENT(in) :: fc
    !
    TYPE(q_grid)  :: in_grid
    !
    COMPLEX(DP) :: U(S%nat3, S%nat3)
    !
    REAL(DP) :: nrg(input%ne), jd_C(input%ne), jd_X(input%ne), &
                xq_i(3), xq_j(3), xq_k(3)
    REAL(DP) :: sigma_ry, weight, e0_ry, de_ry, sigma_e_ry
    REAL(DP) :: freqi(S%nat3), freqj(S%nat3), freqk(S%nat3)
    REAL(DP) :: bosej(S%nat3), bosek(S%nat3), bose_C, bose_X
    REAL(DP) :: dom(input%ne), ctm(input%ne), jdos_X(input%ne), jdos_C(input%ne)
    REAL(DP) :: dom_C(S%nat3), dom_X(S%nat3)
    REAL(DP) :: ctm_C(S%nat3), ctm_X(S%nat3)
    REAL(DP) :: herring_C(S%nat3), aux(S%nat3)
    INTEGER :: iq, jq, k,j,i, SLOW, FAST, fmin, fmax
    CHARACTER (LEN=6), EXTERNAL :: int_to_char
    CHARACTER(len=256) :: filename
    !
    !
    FORALL(i=1:input%ne) nrg(i) = input%de * (i-1) + input%e0
    nrg = nrg/RY_TO_CMM1
    !
    sigma_ry = input%sigma(1)/RY_TO_CMM1
    sigma_e_ry = input%sigma_e/RY_TO_CMM1
    e0_ry = input%e0/RY_TO_CMM1
    de_ry = input%de/RY_TO_CMM1
    !xq0 = input%q_initial

    CALL setup_grid(input%grid_type, S%bg, input%nk(1),input%nk(2),input%nk(3), &
                in_grid, scatter=.true., quiet=.false.)
    jdos_X = 0._dp
    jdos_C = 0._dp
    !
    filename=TRIM(input%outdir)//"/"//&
            TRIM(input%prefix)//"_T"//TRIM(write_conf(1,input%nconf,input%T))//&
                  "_s"//TRIM(write_conf(1,input%nconf,input%sigma))//".out"
    OPEN(unit=10000, file=filename, status="UNKNOWN")
    WRITE(10000, '(a)') "# Total JDOS"
    WRITE(10000, '(a)') "# Energy, tot JDOS, X JDOS, C JDOS"

    filename=TRIM(input%outdir)//"/"//&
            TRIM(input%prefix)//"-q_T"//TRIM(write_conf(1,input%nconf,input%T))//&
                  "_s"//TRIM(write_conf(1,input%nconf,input%sigma))//".out"
    OPEN(unit=30000, file=filename, status="UNKNOWN")
    WRITE(30000, '(a)') "# Q-dependent JDOS"
    WRITE(30000, '(a)') "# i,wq, q, freq, tot JDOS, X JDOS, C JDOS, Herring JDOS (cubic only)"
    !
    ! out_grid should be a proper grid when computing the total energy-dependent JDOS, 
    ! or a path in reciprocal space when doing the local JDOS
    IQ_LOOP : &
    DO iq = 1, out_grid%nq
      !
      CALL print_percent_wall(10._dp, 300._dp, iq, out_grid%nq, (iq==1))
      !
      xq_i = out_grid%xq(:,iq)
      CALL freq_phq(xq_i, S, fc, freqi, U)
      !
      ctm_C = 0._dp
      ctm_X = 0._dp
      herring_C = 0._dp
      !
      JQ_LOOP : &
      DO jq = 1, in_grid%nq
        xq_j = in_grid%xq(:,jq)
        CALL freq_phq(xq_j, S, fc, freqj, U)
        
        xq_k = -(xq_i + xq_j)
        CALL freq_phq(xq_k, S, fc, freqk, U)
        
        IF(input%T(1)>0)THEN
          CALL bose_phq(input%T(1),S%nat3, freqj(:), bosej(:))
          CALL bose_phq(input%T(1),S%nat3, freqk(:), bosek(:))
        ENDIF
        !
        DO k = 1,S%nat3
        DO j = 1,S%nat3
          !
          IF(input%T(1)>0)THEN
            bose_C = 2* (bosej(j) - bosek(k))
            bose_X = bosej(j) + bosek(k) + 1
          ELSE
            bose_C = 1._dp
            bose_X = 1._dp
          ENDIF
          
          dom_C(:) = freqi(:)+freqj(j)-freqk(k) ! cohalescence
          aux(:) = in_grid%w(jq)*bose_C*f_gauss(dom_C, sigma_ry) !delta 
          ctm_C(:) = ctm_C(:)+ aux
          !
          dom_X(:) = freqi(:)-freqj(j)-freqk(k) ! scattering/decay
          ctm_X(:) = ctm_X(:)+ in_grid%w(jq)*bose_X*f_gauss(dom_X, sigma_ry) !delta
          !
          ! Herring processes are important at very low frequency: when the LA mode
          ! collides with a slow TA to form a fast TA. This code only works
          ! assumes that modes 1 and 2 are TA and mode 3 is LA (i.e. GaAs and Silicon)
! ! !           IF(freqk(k)>=freqj(j) .and. ((j==1.and.k==2).or.(j==2.and.k==1))) THEN
! ! !             SLOW=j
! ! !             FAST=k
! ! !           ELSE
! ! !             SLOW=-99
! ! !             FAST=-99
! ! !           ENDIF
! This condition is equivalent, because bands are ordered:
          SLOW = 1
          FAST = 2
          !
          IF( ALL((/j,k/)==(/SLOW,FAST/)) ) THEN
            herring_C(:) = herring_C(:) + aux
          ENDIF
          !
        ENDDO
        ENDDO
        !
      ENDDO &
      JQ_LOOP 
      !
      CALL mpi_bsum(S%nat3,ctm_X)
      CALL mpi_bsum(S%nat3,ctm_C)
      CALL mpi_bsum(S%nat3,herring_C)
      !
      WRITE(30000, '(i9,4f12.6,99e14.6)') iq, out_grid%w(iq), out_grid%xq(:,iq), &
              freqi*RY_TO_CMM1, (ctm_X+ctm_C), ctm_X, ctm_C, herring_C
      !
      weight = out_grid%w(iq)*(input%de/RY_TO_CMM1)
      ctm_X = ctm_X*weight*de_ry
      ctm_C = ctm_C*weight*de_ry
      
      DO i = 1, S%nat3
        !
        ! Only sum a +/- 5 sigma range around the phonon energy
        fmin = MAX(1,        NINT((freqi(i)-e0_ry-5*sigma_e_ry)/de_ry))
        fmax = MIN(input%ne, NINT((freqi(i)-e0_ry+5*sigma_e_ry)/de_ry))
        dom(fmin:fmax) = freqi(i)-nrg(fmin:fmax)
        !ctm(:) = ctm(:) + (ctm_C(i)+ctm_X(i))*weight*f_gauss(dom, sigma_e_ry) 
        ctm(fmin:fmax) = f_gauss(dom(fmin:fmax), sigma_e_ry) 
        jdos_X(fmin:fmax) = jdos_X(fmin:fmax)+ctm(fmin:fmax)*ctm_X(i)
        jdos_C(fmin:fmax) = jdos_C(fmin:fmax)+ctm(fmin:fmax)*ctm_C(i)
        !
!         OPEN(unit=10000, file=TRIM(input%prefix)//"_nu"//&
!                          trim(int_to_char(i))//".out", status="UNKNOWN")
!         WRITE(10000,'(a)') " # energy (cmm1)       total jdos"//&
!                            "            jdos (scattering)     jdos (cohalescence)"
!         DO j = 1,input%ne
!           WRITE(10000,'(4ES27.15E3)') RY_TO_CMM1*nrg(j),ctm(j)*(ctm_X(i)+ctm_C(i)),&
!                                    ctm(j)*ctm_X(i),ctm(j)*ctm_C(i)
!         ENDDO
!         CLOSE(10000)
      ENDDO
      !
    ENDDO IQ_LOOP
    CLOSE(30000)
    !
!     OPEN(unit=10000, file=TRIM(input%prefix)//".out", status="UNKNOWN")
!     WRITE(10000,'(a)') " # energy (cmm1)       total jdos            jdos (scattering)     jdos (cohalescence)"
    DO i = 1,input%ne
      WRITE(10000,'(4ES27.15E3)') RY_TO_CMM1*nrg(i),jdos_X(i)+jdos_C(i),jdos_X(i),jdos_C(i)
    ENDDO
    CLOSE(10000)
    !
  END SUBROUTINE JOINT_DOS
  !
  SUBROUTINE TWOPH_DOS(input, out_grid, S, fc)
    USE code_input,       ONLY : code_input_type
    USE kinds,            ONLY : DP
    USE input_fc,         ONLY : forceconst2_grid, ph_system_info
    USE q_grids,          ONLY : q_grid, q_basis, setup_grid, prepare_q_basis
    USE constants,        ONLY : RY_TO_CMM1, pi
    USE functions,        ONLY : f_bose, f_gauss
    USE fc2_interpolate,  ONLY : freq_phq, bose_phq
    USE mpi_thermal,      ONLY : mpi_bsum, start_mpi, stop_mpi
    USE random_numbers,   ONLY : randy
    USE mpi_thermal,      ONLY : mpi_bsum
    USE more_constants,   ONLY : write_conf
    IMPLICIT NONE
    TYPE(code_input_type)    :: input
    TYPE(q_grid), INTENT(in) :: out_grid
    TYPE(ph_system_info)     :: S
    TYPE(forceconst2_grid),INTENT(in) :: fc
    !
    TYPE(q_grid)  :: in_grid
    !
    COMPLEX(DP) :: U(S%nat3, S%nat3)
    !
    REAL(DP) :: nrg(input%ne), jd_C(input%ne), jd_X(input%ne), &
                xq_i(3), xq_j(3), xq_k(3)
    REAL(DP) :: sigma_ry, weight, e0_ry, de_ry, sigma_e_ry
    REAL(DP) :: freqj(S%nat3), freqk(S%nat3)
    REAL(DP) :: tphdos_X(input%ne), tphdos_C(input%ne)
    REAL(DP) :: dom_C(input%ne), dom_X(input%ne)
    REAL(DP) :: ctm_C(input%ne), ctm_X(input%ne)
    REAL(DP) :: herring_C(S%nat3), aux(S%nat3)
    INTEGER :: iq, jq, k,j,i, SLOW, FAST, fmin, fmax
    CHARACTER (LEN=6), EXTERNAL :: int_to_char
    CHARACTER(len=256) :: filename
    !
    !
    FORALL(i=1:input%ne) nrg(i) = input%de * (i-1) + input%e0
    nrg = nrg/RY_TO_CMM1
    !
    sigma_ry = input%sigma(1)/RY_TO_CMM1
    sigma_e_ry = input%sigma_e/RY_TO_CMM1
    e0_ry = input%e0/RY_TO_CMM1
    de_ry = input%de/RY_TO_CMM1
    !xq0 = input%q_initial

    CALL setup_grid(input%grid_type, S%bg, input%nk(1),input%nk(2),input%nk(3), &
                in_grid, scatter=.true., quiet=.false.)
    !
    filename=TRIM(input%outdir)//"/"//&
            TRIM(input%prefix)//"_T"//TRIM(write_conf(1,input%nconf,input%T))//&
                  "_s"//TRIM(write_conf(1,input%nconf,input%sigma))//".out"
    OPEN(unit=10000, file=filename, status="UNKNOWN")
    WRITE(10000, '(a)') "# Total JDOS"
    WRITE(10000, '(a)') "# Energy, tot JDOS, X JDOS, C JDOS"
    !
    ! out_grid should be a proper grid when computing the total energy-dependent JDOS, 
    ! or a path in reciprocal space when doing the local JDOS
    IQ_LOOP : &
    DO iq = 1, out_grid%nq
      !
      tphdos_X = 0._dp
      tphdos_C = 0._dp
      WRITE(10000, '(a,i3,3f12.6)') "# ", iq, out_grid%xq(:,iq)
      !
      CALL print_percent_wall(10._dp, 300._dp, iq, out_grid%nq, (iq==1))
      !
      xq_i = out_grid%xq(:,iq)
!       CALL freq_phq(xq_i, S, fc, freqi, U)
      !
      ctm_C = 0._dp
      ctm_X = 0._dp
      !
      JQ_LOOP : &
      DO jq = 1, in_grid%nq
        xq_j = in_grid%xq(:,jq)
        CALL freq_phq(xq_j, S, fc, freqj, U)
        
        xq_k = -(xq_i + xq_j)
        CALL freq_phq(xq_k, S, fc, freqk, U)
        
!         IF(input%T(1)>0)THEN
!           CALL bose_phq(input%T(1),S%nat3, freqj(:), bosej(:))
!           CALL bose_phq(input%T(1),S%nat3, freqk(:), bosek(:))
!         ENDIF
        !
        weight = in_grid%w(jq)*(input%de/RY_TO_CMM1)
        !
        DO k = 1,S%nat3
        DO j = 1,S%nat3
          !
          fmin = MAX(1,        NINT((freqj(j)-freqk(k)-e0_ry-5*sigma_e_ry)/de_ry))
          fmax = MIN(input%ne, NINT((freqj(j)-freqk(k)-e0_ry+5*sigma_e_ry)/de_ry))
          dom_C(fmin:fmax) = freqj(j)-freqk(k)-nrg(fmin:fmax) ! cohalescence
          ctm_C(fmin:fmax) = f_gauss(dom_C(fmin:fmax), sigma_e_ry) 
          tphdos_C(fmin:fmax) = tphdos_C(fmin:fmax)+ctm_C(fmin:fmax)*weight
          !
          fmin = MAX(1,        NINT((freqj(j)+freqk(k)-e0_ry-5*sigma_e_ry)/de_ry))
          fmax = MIN(input%ne, NINT((freqj(j)+freqk(k)-e0_ry+5*sigma_e_ry)/de_ry))
          dom_X(fmin:fmax) = freqj(j)+freqk(k)-nrg(fmin:fmax) ! decay
          ctm_X(fmin:fmax) = f_gauss(dom_X(fmin:fmax), sigma_e_ry) 
          tphdos_X(fmin:fmax) = tphdos_X(fmin:fmax)+ctm_X(fmin:fmax)*weight
          !
        ENDDO
        ENDDO
        !
      ENDDO &
      JQ_LOOP 
      !
      CALL mpi_bsum(input%ne,tphdos_X)
      CALL mpi_bsum(input%ne,tphdos_C)
      !
      weight = out_grid%w(iq)*(input%de/RY_TO_CMM1)
      ctm_X = ctm_X*weight*de_ry
      ctm_C = ctm_C*weight*de_ry
      
!      DO i = 1, S%nat3
!        !
!        ! Only sum a +/- 5 sigma range around the phonon energy
!        !
!      ENDDO
      !
      DO i = 1,input%ne
        WRITE(10000,'(1f12.6, 4ES27.15E3)') out_grid%w(iq), &
            RY_TO_CMM1*nrg(i),tphdos_X(i)+tphdos_C(i), tphdos_X(i),tphdos_C(i)
      ENDDO
      CLOSE(10000)
    ENDDO IQ_LOOP
    !
    !
  END SUBROUTINE TWOPH_DOS  !
  !
  SUBROUTINE EXTRINSIC_LW(input, qpath, S, fc2)
    USE kinds,              ONLY : DP
    USE fc2_interpolate,    ONLY : fftinterp_mat2, mat2_diag, freq_phq_path
    USE linewidth,          ONLY : linewidth_q, selfnrg_q, spectre_q
    USE constants,          ONLY : RY_TO_CMM1
    USE q_grids,            ONLY : q_grid, setup_grid
    USE more_constants,     ONLY : write_conf
    USE isotopes_linewidth, ONLY : isotopic_linewidth_q
    USE casimir_linewidth,  ONLY : casimir_linewidth_q
    USE input_fc,           ONLY : forceconst2_grid, ph_system_info
    USE code_input,         ONLY : code_input_type
    USE overlap,            ONLY : order_type
    USE timers
    IMPLICIT NONE
    !
    TYPE(code_input_type),INTENT(in)     :: input
    TYPE(forceconst2_grid),INTENT(in) :: fc2
    TYPE(ph_system_info),INTENT(in)   :: S
    TYPE(q_grid),INTENT(in)      :: qpath
    !
    TYPE(order_type) :: order
    COMPLEX(DP) :: D(S%nat3, S%nat3)
    REAL(DP) :: w2(S%nat3)
    INTEGER :: iq, it
    TYPE(q_grid) :: grid
    REAL(DP)   :: lw(S%nat3,input%nconf), lw_isot(S%nat3,input%nconf)
    REAL(DP)   :: lw_casimir(S%nat3)
    REAL(DP)   :: sigma_ry(input%nconf)
    CHARACTER(len=32) :: f1, f2
    CHARACTER(len=256) :: filename
    CHARACTER(len=6) :: pos 
    !
    CALL setup_grid(input%grid_type, S%bg, input%nk(1), input%nk(2), input%nk(3), &
                    grid, scatter=.true., xq0=input%xk0, quiet=.false.)
    !
    IF(ionode)THEN
      IF(input%skip_q>0) THEN; pos="append"; ELSE; pos = "asis"; ENDIF
      DO it = 1,input%nconf
        filename=TRIM(input%outdir)//"/"//&
                TRIM(input%prefix)//"_T"//TRIM(write_conf(it,input%nconf,input%T))//&
                  "_s"//TRIM(write_conf(it,input%nconf,input%sigma))//".out"
        OPEN(unit=1000+it, file=filename, position=pos)
        ioWRITE(1000+it, *) "# calculation of extrinsic linewidth: total ... isotopes ... casimir"
        ioWRITE(1000+it, '(a,i6,a,f6.1,a,100f6.1)') "# ", it, &
                         "     T=",input%T(it), "    sigma=", input%sigma(it)
        ioFLUSH(1000+it)
      ENDDO
      ! Prepare formats to write out data
      ioWRITE(f1,'(i6,a)') S%nat3, "f12.6,6x,"
      ioWRITE(f2,'(i6,a)') S%nat3, "ES27.15E3,6x,"
    ENDIF
    !
    sigma_ry = input%sigma/RY_TO_CMM1
    !
    ioWRITE(*,'(2x,a,i6,a)') "Going to compute", qpath%nq, " points (1)"
    !
    DO iq = input%skip_q +1,qpath%nq
      !
      CALL print_percent_wall(10._dp, 300._dp, iq, qpath%nq, (iq==1))
      !
      CALL freq_phq_path(qpath%nq, iq, qpath%xq, S, fc2, w2, D)
      !
      ! If necessary, compute the ordering of the bands to assure modes continuity;
      ! on first call, it just returns the trivial 1...3*nat order
      !IF(input%sort_freq=="overlap" .or. iq==1) CALL order%set(S%nat3, w2, D)
      IF(input%sort_freq=="overlap" .or. iq==1) &
          CALL order%set_path(S%nat3, w2, D, iq, qpath%nq, qpath%w, qpath%xq)
      !
      ! Compute isotopic linewidth
      IF(input%isotopic_disorder)THEN
          timer_CALL t_lwisot%start()
        lw_isot = isotopic_linewidth_q(qpath%xq(:,iq), input%nconf, input%T,&
                                        sigma_ry, S, grid, fc2, w2, D)
          timer_CALL t_lwisot%stop()
      ENDIF
      !
      ! Compute Casimir finite sample size linewidth
      IF(input%casimir_scattering)THEN
          timer_CALL t_lwcasi%start()
        lw_casimir = casimir_linewidth_q(qpath%xq(:,iq), input%sample_length, &
                                      input%sample_dir, S, fc2)
        timer_CALL t_lwcasi%stop()
      ENDIF
      !
      DO it = 1,input%nconf
        !
        ! Add isotopic and casimir linewidths
        lw(:,it) = lw_isot(:,it) + lw_casimir
        !
        ioWRITE(1000+it, '(i4,f12.6,2x,3f12.6,2x,'//f1//f2//f2//f2//'x)') &
              iq,qpath%w(iq),qpath%xq(:,iq), &
              w2(order%idx(:))*RY_TO_CMM1, lw(order%idx(:),it)*RY_TO_CMM1, &
              lw_isot(order%idx(:),it)*RY_TO_CMM1, lw_casimir(order%idx(:))*RY_TO_CMM1
        ioFLUSH(1000+it)
      ENDDO
      !
      !
    ENDDO
    !
    !
    IF(ionode)THEN
    DO it = 1,input%nconf
      CLOSE(unit=1000+it)
    ENDDO
    ENDIF
    !
    CALL print_all_timers()
    !
  END SUBROUTINE EXTRINSIC_LW
  !
  !
  !
  SUBROUTINE PH_DOS(input, S, fc)
    USE code_input,       ONLY : code_input_type
    USE kinds,            ONLY : DP
    USE input_fc,         ONLY : forceconst2_grid, ph_system_info
    USE q_grids,          ONLY : q_grid, q_basis, setup_grid, prepare_q_basis
    USE constants,        ONLY : RY_TO_CMM1, pi
    USE functions,        ONLY : f_bose, f_gauss
    USE fc2_interpolate,  ONLY : freq_phq
    USE mpi_thermal,      ONLY : mpi_bsum, start_mpi, stop_mpi
    USE random_numbers,   ONLY : randy
    IMPLICIT NONE
    TYPE(code_input_type) :: input
    TYPE(ph_system_info)   :: S
    TYPE(forceconst2_grid),INTENT(in) :: fc
    !
    TYPE(q_grid)  :: qgrid
    !
    REAL(DP) :: freqj(S%nat3)
    COMPLEX(DP) :: U(S%nat3, S%nat3)
    !
    !REAL(DP) :: xq_random(3)
    !
    REAL(DP) :: nrg(input%ne), xq_j(3)
    REAL(DP) :: sigma_ry, weight, e0_ry, de_ry
    REAL(DP) :: dos(input%ne), dom(input%ne)
    INTEGER :: jq, k,j,i, fmin, fmax
    !
    !
    FORALL(i=1:input%ne) nrg(i) = input%de * (i-1) + input%e0
    nrg = nrg/RY_TO_CMM1
    !
    !sigma_ry = input%sigma(1)/RY_TO_CMM1
    sigma_ry = input%sigma_e/RY_TO_CMM1
    e0_ry = input%e0/RY_TO_CMM1
    de_ry = input%de/RY_TO_CMM1
    
    dos = 0._dp

    !xq_random  = (/ randy(), randy(), randy() /)
    CALL setup_grid(input%grid_type, S%bg, input%nk(1),input%nk(2),input%nk(3), &
                qgrid, scatter=.true.)
    
    DO jq = 1, qgrid%nq
      xq_j = qgrid%xq(:,jq)
      CALL freq_phq(xq_j, S, fc, freqj, U)
      !
      weight = qgrid%w(jq)*(input%de/RY_TO_CMM1)
      !
      DO j = 1,S%nat3
        !
!         dom(:) =freqj(j)-nrg(:)
!         dos = dos + weight * f_gauss(dom, sigma_ry) 
        fmin = MAX(1,        NINT((freqj(j)-e0_ry-5*sigma_ry)/de_ry))
        fmax = MIN(input%ne, NINT((freqj(j)-e0_ry+5*sigma_ry)/de_ry))
        dom(fmin:fmax) =freqj(j)-nrg(fmin:fmax)
        dos(fmin:fmax) = dos(fmin:fmax) + weight * f_gauss(dom(fmin:fmax), sigma_ry) 
        !
      ENDDO
      !
    ENDDO
    CALL mpi_bsum(input%ne,dos)
    !
    OPEN(unit=10000, file=TRIM(input%prefix)//".out", status="UNKNOWN")
    DO i = 1,input%ne
      WRITE(10000,'(4ES27.15E3)') RY_TO_CMM1*nrg(i),dos(i)
    ENDDO
    CLOSE(10000)
    !
  END SUBROUTINE PH_DOS
  !
  !
  FUNCTION f_factor(input, S, xq) RESULT(f)
    USE code_input,       ONLY : code_input_type
    USE kinds,            ONLY : DP
    USE input_fc,         ONLY : forceconst2_grid, ph_system_info
    IMPLICIT NONE
    TYPE(code_input_type) :: input
    TYPE(ph_system_info)   :: S
    REAL(DP) :: f(S%ntyp)
    !
    REAL(DP) :: xq(3), d2,d3,d4, integral
    REAL(DP),ALLOCATABLE :: r(:), charge(:), bes(:)
    INTEGER :: i, j, nlines
    !
    nlines = 1095
    DO i =  1,S%ntyp
      !
      ALLOCATE(r(nlines), charge(nlines), bes(nlines))
      IF(i==1)THEN
        OPEN(unit=400, file="O_core.dat", status="old", action="read")
      ELSEIF(i==2)THEN
        OPEN(unit=400, file="H_core.dat", status="old", action="read")
      ENDIF
      DO j = 1,nlines
        READ(400,*) r(j), d2, d3, d4
        charge(j) = d3+d4
      ENDDO
      CALL sph_bes (nlines, r, xq, 0, bes)
      charge = bes*charge
      CALL simpson(nlines, charge, r, integral)
      !
      f(i) = integral
      !      
      DEALLOCATE(r,charge,bes)
      CLOSE(400)
    ENDDO
    
  END FUNCTION
  !
  SUBROUTINE DYNFF_S(input,listq, S, fc)
    USE code_input,       ONLY : code_input_type
    USE kinds,            ONLY : DP
    USE input_fc,         ONLY : forceconst2_grid, ph_system_info
    USE q_grids,          ONLY : q_grid, q_basis, setup_grid, prepare_q_basis
    USE constants,        ONLY : RY_TO_CMM1, pi, tpi
    USE functions,        ONLY : f_bose
    USE fc2_interpolate,  ONLY : freq_phq
    USE mpi_thermal,      ONLY : mpi_bsum, start_mpi, stop_mpi
    IMPLICIT NONE
    TYPE(code_input_type) :: input
    TYPE(q_grid), INTENT(in) :: listq
    TYPE(ph_system_info)   :: S
    TYPE(forceconst2_grid),INTENT(in) :: fc
    !
    TYPE(q_grid)  :: qgrid
    !
    COMPLEX(DP) :: U(S%nat3, S%nat3), phase, Fin
    REAL(DP) :: ff(S%ntyp), freq(S%nat3), rr(3,S%nat), arg, e_jat(3), &
                Fin2, loss, gain, n, xq(3)
    !
    !
    INTEGER :: i, j, iq, jat, nu
    !
!     CALL setup_grid(input%grid_type, S%bg, input%nk(1),input%nk(2),input%nk(3), &
!                 qgrid, scatter=.true.)

     ! take r vectors inside the first unit cell
     rr = S%tau
!      CALL cryst_to_cart(S%nat, rr, S%bg, -1)
!      rr = rr - NINT(rr)
!      CALL cryst_to_cart(S%nat, rr, S%at, +1)
    
    OPEN(unit=10000, file=TRIM(input%prefix)//".out", status="UNKNOWN")
    DO iq = 1, listq%nq
      !
      xq = listq%xq(:,iq)
      ff = f_factor(input, S, xq)
      !
      CALL freq_phq(xq, S, fc, freq, U)
      WRITE(10000,'(3f12.6)',advance='no') xq
      !
      DO nu = 1,S%nat3
        !
        n = f_bose(freq(nu), input%T(1)) 
        !
        Fin = 0._dp
        DO jat = 1,S%nat
                !
                !jat = (j-1)/3+1
                !jalpha = j-(jat-1)*3
                arg = tpi * SUM(xq*rr(:,jat))
                phase = CMPLX(Cos(arg),-Sin(arg), kind=DP)
                !
                e_jat = U((jat-1)*3+1 : (jat-1)*3+3,  nu)
                Fin = Fin + ff(S%ityp(jat)) * phase * SUM(xq*e_jat) / DSQRT(S%amass(S%ityp(jat)))
                !
        ENDDO
        Fin2 = DBLE(Fin * CONJG(Fin))
        !
        loss = (n+1)/freq(nu) * Fin2
        gain = n/freq(nu) * Fin2
        WRITE(10000,'(5x,3ES27.15E3)',advance='no') freq(nu), loss, gain
        !
      ENDDO
      !
      WRITE(10000,*)
      !
    ENDDO
    !
    CLOSE(10000)
    !
  END SUBROUTINE DYNFF_S
  !
  !
  !
  SUBROUTINE RMS(input, S, fc)
    USE code_input,       ONLY : code_input_type
    USE kinds,            ONLY : DP
    USE input_fc,         ONLY : forceconst2_grid, ph_system_info
    USE q_grids,          ONLY : q_grid, setup_grid
    !USE constants,        ONLY : RY_TO_CMM1, pi
    USE functions,        ONLY : f_wtoa !f_bose, f_gauss
    USE fc2_interpolate,  ONLY : freq_phq_safe, set_nu0
    !USE random_numbers,   ONLY : randy
    IMPLICIT NONE
    TYPE(code_input_type) :: input
    TYPE(ph_system_info)   :: S
    TYPE(forceconst2_grid),INTENT(in) :: fc
    !
    TYPE(q_grid)  :: qgrid
    REAL(DP) :: freqj(S%nat3), aq(S%nat3), arms(S%nat3), xq_j(3)
    COMPLEX(DP) :: U(S%nat3, S%nat3)
    INTEGER :: jq, ia, nu, mu, mu0
    
    CALL setup_grid(input%grid_type, S%bg, input%nk(1),input%nk(2),input%nk(3), &
                qgrid, scatter=.false.)
    !
    arms = 0._dp
    DO jq = 1, qgrid%nq
      xq_j = qgrid%xq(:,jq)
      CALL freq_phq_safe(xq_j, S, fc, freqj, U)
      
      !bosej(:) = f_bose(freqj, input%T(1))
      
      aq(:) = f_wtoa(freqj, input%T(1))
      
      mu0  = set_nu0(xq_j, S%at)
      
      DO ia = 1, S%nat
        nu = (ia-1)*3 + 1
        DO mu = mu0, S%nat3
          
          arms(ia) = arms(ia) + aq(mu)**2 &
                       *DBLE(SUM(U(nu:nu+2,mu)*CONJG(U(nu:nu+2,mu)))) &
                              /S%amass(S%ityp(ia)) * qgrid%w(jq)
        ENDDO
      ENDDO
    ENDDO
    !
    WRITE(*,'(a)') "  atm    sqrt(rms) [bohr]"
    DO ia = 1, S%nat
      WRITE(*,'(i3,x,a3,2f12.6)') ia, S%atm(S%ityp(ia)),&
                                  DSQRT(arms(ia))
    ENDDO
    
  END SUBROUTINE RMS
  !
END MODULE r2q_program
!
!
!
PROGRAM r2q 

  USE kinds,            ONLY : DP
  USE r2q_program
  USE input_fc,         ONLY : read_fc2, aux_system, div_mass_fc2, &
                              forceconst2_grid, ph_system_info, &
                              multiply_mass_dyn, write_dyn
  USE asr2_module,      ONLY : impose_asr2
  USE constants,        ONLY : RY_TO_CMM1
  USE fc2_interpolate,  ONLY : freq_phq, freq_phq_path, fftinterp_mat2
  USE q_grids,          ONLY : q_grid
  USE code_input,       ONLY : code_input_type, READ_INPUT
  USE ph_velocity,      ONLY : velocity
  USE more_constants,   ONLY : print_citations_linewidth
  USE overlap,          ONLY : order_type
  USE mpi_thermal,      ONLY : start_mpi, stop_mpi, ionode
  USE nanoclock,        ONLY : init_nanoclock
  USE neutrons,         ONLY : neutron_cross_section
  IMPLICIT NONE
  !
  TYPE(forceconst2_grid) :: fc2
  TYPE(ph_system_info)   :: S
  TYPE(q_grid)           :: qpath
  TYPE(code_input_type)  :: input
  !
  CHARACTER(len=512) :: filename
  !
  !REAL(DP) :: xq(3)
  REAL(DP),ALLOCATABLE :: freq(:), vel(:,:), proj(:,:)
  REAL(DP) :: aux
  COMPLEX(DP),ALLOCATABLE :: U(:,:), D(:,:)
  INTEGER :: i, output_unit=10000, nu, ia
  TYPE(order_type) :: order
  CHARACTER (LEN=6),  EXTERNAL :: int_to_char
  !
  CALL start_mpi()
  CALL init_nanoclock()
  IF(ionode) CALL print_citations_linewidth()
  !  
  CALL READ_INPUT("R2Q", input, qpath, S, fc2)
  !
!   IF(input%nconf>1) THEN
!     CALL errore("R2Q", "r2q.x only supports one configuration at a time.",1)
!   ENDIF

  IF( input%calculation=="dos") THEN
    CALL PH_DOS(input,S,fc2)
  ELSE IF ( input%calculation=="dynff") THEN
    CALL DYNFF_S(input,qpath, S, fc2)
  ELSE IF( input%calculation=="jdos") THEN
    CALL JOINT_DOS(input,qpath,S,fc2)
  ELSE IF( input%calculation=="2dos") THEN
    CALL TWOPH_DOS(input,qpath,S,fc2)
  ELSE IF ( input%calculation=="rms") THEN
    CALL RMS(input, S, fc2)
  ELSE IF ( input%calculation=="extr") THEN
    CALL EXTRINSIC_LW(input, qpath, S, fc2)
  ELSE
    ALLOCATE(freq(S%nat3))
    ALLOCATE(U(S%nat3,S%nat3))
    IF(input%print_dynmat) ALLOCATE(D(S%nat3,S%nat3))

    filename=TRIM(input%outdir)//"/"//TRIM(input%prefix)//".out"
    OPEN(unit=output_unit, file=filename)

    IF(input%print_velocity) THEN
      filename=TRIM(input%outdir)//"/"//TRIM(input%prefix)//"_vel.out"
      OPEN(unit=output_unit+1, file=filename)
      ALLOCATE(vel(3,S%nat3))
    ENDIF

    IF(input%print_neutron_cs) THEN
      filename=TRIM(input%outdir)//"/"//TRIM(input%prefix)//"_ins.out"
      OPEN(unit=output_unit+2, file=filename)
    ENDIF
    !
    ALLOCATE(proj(S%ntyp,S%nat3))
    DO i = 1,qpath%nq
      !CALL freq_phq(qpath%xq(:,i), S, fc2, freq, U)
      CALL freq_phq_path(qpath%nq, i, qpath%xq, S, fc2, freq, U)
      IF(input%sort_freq=="overlap" .or. i==1) & !CALL order%set(S%nat3, freq, U)
          CALL order%set_path(S%nat3, freq, U, i, qpath%nq, qpath%w, qpath%xq)

      ! project on atomic types
      proj=0._dp
      DO nu = 1,S%nat3
      DO ia = 1, S%nat
        aux = SUM(U((ia-1)*3+1:(ia-1)*3+3, nu)*CONJG(U((ia-1)*3+1:(ia-1)*3+3, nu)))
        proj(S%ityp(ia),nu) = proj(S%ityp(ia),nu) + aux
      ENDDO
      ENDDO

      ioWRITE(output_unit, '(i6,f12.6,3x,3f12.6,999e16.6)') &
        i, qpath%w(i), qpath%xq(:,i), freq(order%idx(:))*RY_TO_CMM1, proj(:,order%idx(:))
      ioFLUSH(output_unit)
      
      IF(input%print_dynmat) THEN
        CALL fftinterp_mat2(qpath%xq(:,i), S, fc2, D)
        D = multiply_mass_dyn(S, D)
        filename = TRIM(input%outdir)//"/"//TRIM(input%prefix)//"_dyn"//TRIM(int_to_char(i))
        CALL write_dyn(filename, qpath%xq(:,i), U, S)
!         DO nu = 1,S%nat3
!           WRITE(stdout, '(99(2f10.4,2x))') U(:,nu)
!         ENDDO
      ENDIF

      IF(input%print_velocity) THEN
        vel = velocity(S, fc2, qpath%xq(:,i))
        ioWRITE(output_unit+1, '(i6,f12.6,3x,3f12.6,999(3e16.8,3x))') &
          i, qpath%w(i), qpath%xq(:,i), vel(:,order%idx(:))
        ioFLUSH(output_unit+1)
      ENDIF
      
      IF(input%print_neutron_cs)THEN
        ioWRITE(output_unit+2, '(i6,f12.6,3x,3f12.6,999e16.8)') &
          i, qpath%w(i), qpath%xq(:,i), neutron_cross_section( qpath%xq(:,i), DBLE((/0,0,0/)), S, fc2)
      ENDIF

    ENDDO
    !
    CLOSE(output_unit)
    DEALLOCATE(freq, U)
    IF(input%print_dynmat) DEALLOCATE(D)
    IF(input%print_velocity) THEN
      CLOSE(output_unit+1)
      DEALLOCATE(vel)
    ENDIF
    !
  ENDIF
  CALL stop_mpi()
  
END PROGRAM r2q
