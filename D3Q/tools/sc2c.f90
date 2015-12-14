!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! 
!
!
!----------------------------------------------------------------------------
PROGRAM sc2c
  !----------------------------------------------------------------------------
  !  The name and order of files is not important as long as q=0 is the first
  !
  USE sc2c_params
  USE mp,         ONLY : mp_start, mp_end, mp_barrier
  USE mp_global,  ONLY : nproc, mpime
  USE fft_scalar, ONLY : cfft3d
  !
  IMPLICIT NONE
  !
  CHARACTER(len=4),PARAMETER :: prog = 'sc2c'
  CHARACTER(len=256) :: fildyn, fildynSC !filin, filj, filf, flfrc
  !
  TYPE(ph_system_info),TARGET  :: CQ, SG
  TYPE(ph_system_info),POINTER :: PT

  INTEGER,ALLOCATABLE :: atoms_map(:)
  INTEGER :: ierr, i_c, i_sc, n_translations_found = 0
  INTEGER,PARAMETER :: u = 1
  INTEGER :: index_in_list
  REAL(DP),ALLOCATABLE :: translations(:,:)
  REAL(DP) :: delta(3), q(3), qdr, work
  INTEGER  :: volume_ratio
  INTEGER :: i,j,k, at_i, at_j, c_at_i, c_at_j
  COMPLEX(DP),ALLOCATABLE :: refolded_D(:,:,:,:)
  COMPLEX(DP),PARAMETER :: ii = (0._dp, 1._dp)
  !
  NAMELIST / input / fildyn, fildynSC
  CALL input_from_file ( ) 
  READ ( 5, input )

  DO i = 1,2
    IF( i==1 ) THEN
      PT => CQ
      OPEN (unit=u, file=TRIM(fildyn),   status='old', form='formatted', iostat=ierr)
      CALL errore(prog, 'Cannot open '//TRIM(fildyn), ierr)
      WRITE(6,'(2x,a,i3)') 'Reading file "'//TRIM(fildyn)//'"', i
    ELSE
      PT => SG
      OPEN (unit=u, file=TRIM(fildynSC), status='old', form='formatted', iostat=ierr)
      CALL errore(prog, 'Cannot open '//TRIM(fildynSC), ierr)
      WRITE(6,'(2x,a,i3)') 'Reading file "'//TRIM(fildynSC)//'"', i
    ENDIF

    CALL read_file_sc2c (u,   PT%nqs,  PT%q,   PT%epsil,  PT%lrigid, &
                    PT%ntyp,  PT%ityp, PT%nat, PT%tau,    PT%zeu,    &
                    PT%ibrav, PT%symm_type,    PT%celldm, PT%at,     &
                    PT%atm,   PT%amass,PT%phiq )
    CLOSE(u )
    !
    CALL latgen(PT%ibrav,PT%celldm,PT%at(1,1),PT%at(1,2),PT%at(1,3),PT%omega)
    PT%at = PT%at / PT%celldm(1)  !  bring at in units of alat 
    !
    CALL volume(PT%celldm(1),PT%at(1,1),PT%at(1,2),PT%at(1,3),PT%omega)
    CALL recips(PT%at(1,1),PT%at(1,2),PT%at(1,3),PT%bg(1,1),PT%bg(1,2),PT%bg(1,3))
    !
    PT%at  = PT%at  * PT%celldm(1)  ! bring at back to bohr units
    PT%tau = PT%tau * PT%celldm(1)  ! also atomic positions in bohr units
    PT%bg  = PT%bg  / PT%celldm(1)  ! reciprocal lattice in 2pi/bohr
    PT%q   = tpi*PT%q/PT%celldm(1)  ! q-point in units of 2pi/bohr
  ENDDO
  !
  IF( ANY(ABS(SG%q(:,1)) > eps8 ) ) &
    CALL errore(prog, 'Supercell must be at Gamma', 1)
  !
  ! Find the translation vectors
  volume_ratio = INT(SG%omega/CQ%omega +0.5_dp)
  WRITE(6, '(2x,a,i3)') "Supercell/cell volume ratio:", volume_ratio
  ALLOCATE( translations(3, volume_ratio))
  ALLOCATE( atoms_map(SG%nat) )
  !
  WRITE(6, '(2x,a)') "Internal SC fractional traslations found (bohr):"
  Do i_c = 1, CQ%nat
    Do i_sc = 1, SG%nat
      !
      delta = SG%tau(:,i_sc) - CQ%tau(:,i_c)
      IF(check_int_linearcombination(delta, CQ%bg)) THEN
        !
        atoms_map(i_sc) = i_c
        !
        index_in_list = check_which_in_list( delta, translations, n_translations_found, SG%at)
        !
        IF(index_in_list<0) THEN
          write(6,'(4x,3f12.6,2x,2i3)') delta, i_c, i_sc
          n_translations_found = n_translations_found+1
          !
          IF(n_translations_found > volume_ratio) &
            CALL errore(prog, 'Found more translations than expected', n_translations_found)
          !
          translations(:,n_translations_found) = delta
          !
        ENDIF 
        ! 
      ENDIF
    ENDDO
  ENDDO
  !
  IF(n_translations_found < volume_ratio) &
    CALL errore(prog, 'Found less translations than expected', n_translations_found)
  !
  WRITE(*, '("Possible choices for good q-points:")')
  WRITE(*, '("(Given as reference only! (2pi/alat <- of unit cell)")')
  WRITE(*, '(4x,3f14.9)') SG%bg(:,1) * CQ%celldm(1)
  WRITE(*, '(4x,3f14.9)') SG%bg(:,2) * CQ%celldm(1)
  WRITE(*, '(4x,3f14.9)') SG%bg(:,3) * CQ%celldm(1)
  WRITE(*, '(4x,3f14.9)') (SG%bg(:,1)+SG%bg(:,2)) * CQ%celldm(1) 
  WRITE(*, '(4x,3f14.9)') (SG%bg(:,1)+SG%bg(:,3)) * CQ%celldm(1) 
  WRITE(*, '(4x,3f14.9)') (SG%bg(:,2)+SG%bg(:,3)) * CQ%celldm(1) 
  WRITE(*, '(4x,3f14.9)') (SG%bg(:,1)+SG%bg(:,2)+SG%bg(:,3)) * CQ%celldm(1)
  !
!  write(*,*)
!  WRITE(*, '(4x,3f12.7)') CQ%bg(:,1) * CQ%celldm(1)
!  WRITE(*, '(4x,3f12.7)') CQ%bg(:,2) * CQ%celldm(1)
!  WRITE(*, '(4x,3f12.7)') CQ%bg(:,3) * CQ%celldm(1)
!  WRITE(*, '(4x,3f12.7)') (CQ%bg(:,1)+CQ%bg(:,2)) * CQ%celldm(1)
!  WRITE(*, '(4x,3f12.7)') (CQ%bg(:,1)+CQ%bg(:,3)) * CQ%celldm(1)
!  WRITE(*, '(4x,3f12.7)') (CQ%bg(:,2)+CQ%bg(:,3)) * CQ%celldm(1)
!  WRITE(*, '(4x,3f12.7)') (CQ%bg(:,1)+CQ%bg(:,2)+CQ%bg(:,3)) * CQ%celldm(1)
 
  !
  WRITE(*, '("Refolding dynamical matrix at q-point (2pi/bohr):")')
  q = CQ%q(:,1) !SG%bg(:,1)
  WRITE(*, '(4x,3f12.5)') q
  !
  ! Check if the q point respects the translational symmetry of the supercell
  IF( ABS(EXP(ii* SUM(q*SG%at(:,1)) )-1._dp) &
     +ABS(EXP(ii* SUM(q*SG%at(:,2)) )-1._dp) &
     +ABS(EXP(ii* SUM(q*SG%at(:,3)) )-1._dp) > eps8 ) &
  THEN
    CALL errore(prog, 'Supercell does not match the q-point', 2)
  ENDIF
  !
  ALLOCATE(refolded_D(3,3,CQ%nat,CQ%nat))
  refolded_D = 0._dp
  !
  DO at_i = 1,SG%nat
  DO at_j = 1,SG%nat
    c_at_i = atoms_map(at_i)
    c_at_j = atoms_map(at_j)
    DO i = 1,3
    DO j = 1,3
      !
!       delta = (SG%tau(:,at_j)-SG%tau(:,at_i)) - (CQ%tau(:,c_at_j)-CQ%tau(:,c_at_i))
!       qdr= (q(1)*delta(1) + q(2)*delta(2) + q(3)*delta(3))
      delta = (SG%tau(:,at_j)-CQ%tau(:,c_at_j)) - (SG%tau(:,at_i)-CQ%tau(:,c_at_i))
      qdr= (q(1)*delta(1) + q(2)*delta(2) + q(3)*delta(3))
!       write(*,'(4i3,4f12.6)') at_i, at_j, c_at_i, c_at_j, delta, qdr/tpi
      !
      refolded_D(i,j,c_at_i,c_at_j) = refolded_D(i,j,c_at_i,c_at_j) &
        + SG%phiq(i,j,at_i,at_j,1) * exp(ii * qdr )/volume_ratio
      !
    ENDDO
    ENDDO
  ENDDO
  ENDDO

  WRITE(*, '("Dynamical matrix:")')
  WRITE(*, '(a)') &
'alp1 alp2 at1 at2     refolded D                original D                diff        diff%'
  DO at_i = 1,CQ%nat
  DO at_j = 1,CQ%nat
    DO i = 1,3
    DO j = 1,3
      !
      work = ABS(refolded_D(i,j,at_i,at_j)-CQ%phiq(i,j,at_i,at_j,1))
      work = work /(ABS(refolded_D(i,j,at_i,at_j))+ABS(CQ%phiq(i,j,at_i,at_j,1)))
      work = 100*work
      WRITE(*,'(2(2i3,2x),2x,2(2f12.6,2x),f12.6,f10.4)') i,j,at_i,at_j,refolded_D(i,j,at_i,at_j), CQ%phiq(i,j,at_i,at_j,1), &
                                                                    ABS(refolded_D(i,j,at_i,at_j)-CQ%phiq(i,j,at_i,at_j,1)), &
                                                                    work
      !
    ENDDO
    ENDDO
  ENDDO
  ENDDO


END PROGRAM sc2c


