!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
PROGRAM d3_sc2c
  !----------------------------------------------------------------------------
  !  wat20100225:
  !  Sketch of a program that 
  !    - reads the third order dynamical matrix of a supercell for q1=0, q2=-q, q3=q
  !         (generated with the old D3(0,-q,q) code)
  !    - reads the third order dymanical matrix of the unitcell for some q1, q2, q3
  !         (generated with the new D3(q1,q2,q3) code)
  !    - checks that the unitcell indeed is a unit cell of the supercell
  !    - checks that the dynamical matrix of the supercell can be backfolded to
  !         the unit cell for the (unit cell) q-vectors q1, q2, q3
  !    - if so, fold back the dynamical matrix of the supercell to the unit cell
  !    - write a comparison of the original unit cell dynamical matrix and the 
  !         backfolded one to the output
  !
  !   wat20100225:
  !   Using this code, one should take into consideration, that so far only the
  !   Ewald contribution has been implemented in the new code.
  !
  USE d3_tools
  !
  ! wat20100225: Clean up the things we don't use
  ! wat20100225: USE sc2c_params, ONLY : read_file_sc2c, ph_system_info
  !USE mp,         ONLY : mp_start, mp_end, mp_barrier
!  USE fft_scalar, ONLY : cfft3d
  !USE iotk_module
  USE d3matrix_io2
  !
  IMPLICIT NONE
  !
  CHARACTER(len=4),PARAMETER :: prog = 'd3_sc2c'
  CHARACTER(len=256) :: fileD3SC, fileD3UC, file
  !
  TYPE(d3_system_info),TARGET  :: D3SC, D3UC
  TYPE(d3_system_info),POINTER :: PT


  INTEGER,ALLOCATABLE :: atoms_map(:)
  INTEGER :: ierr, i_c, i_sc, n_translations_found = 0
  INTEGER,PARAMETER :: u = 1
  INTEGER :: index_in_list
  REAL(DP),ALLOCATABLE :: translations(:,:)
  REAL(DP) :: delta(3), delta1(3), delta2(3), delta3(3), q1r1, q2r2, q3r3, &
              work, MAXWORK,MAXWORKd, diff, MAXDIFF, MAXDIFFp
  INTEGER  :: volume_ratio
  INTEGER :: i,j,k, at_i, at_j, at_k, c_at_i, c_at_j, c_at_k
  COMPLEX(DP),ALLOCATABLE :: refolded_D3(:,:,:,:,:,:)
  COMPLEX(DP),PARAMETER :: ii = (0._dp, 1._dp)
  COMPLEX(DP) :: phase
  !
  INTEGER,INTRINSIC :: iargc
  INTEGER :: nargs
  !
  nargs = iargc()
  !
  IF(nargs>=1) THEN
    CALL getarg(1, fileD3UC)
  ELSE
    PRINT*, "Note: you can also pass the file names on command line"
    PRINT*, "enter unit cell file"
    READ(*,"(a256)") fileD3UC
  ENDIF
  IF(nargs>=2) THEN
    CALL getarg(2, fileD3SC)
  ELSE
    PRINT*, "enter super-cell file"
    READ(*,"(a256)") fileD3SC
  ENDIF

  ! Read the third order dynamical matrices of the supercell and the unit cell
  DO i = 1,2
    IF( i==1 ) THEN
      PT => D3SC
      file = fileD3SC
!       OPEN (unit=u, file=TRIM(fileD3SC),   status='old', form='formatted', iostat=ierr)
!       CALL errore(prog, 'Cannot open '//TRIM(fileD3SC), ierr)
!       WRITE(6,'(2x,a,i3)') 'Reading third order super cell dynamical matrix from file "'//TRIM(fileD3SC)//'"', i
    ELSE
      PT => D3UC
      file = fileD3UC
!       OPEN (unit=u, file=TRIM(fileD3UC), status='old', form='formatted', iostat=ierr)
!       CALL errore(prog, 'Cannot open '//TRIM(fileD3UC), ierr)
!       WRITE(6,'(2x,a,i3)') 'Reading third order unit cell  dynamical matrix from file "'//TRIM(fileD3UC)//'"', i
    ENDIF

    CALL read_d3dyn_xml2(file, PT%q1,PT%q2,PT%q3, PT%phiq, PT%ntyp, PT%nat, PT%ibrav, PT%celldm, &
                        PT%at, PT%ityp, PT%tau, PT%atm, PT%amass)

    CLOSE(u )
    !
    CALL latgen(PT%ibrav,PT%celldm,PT%at(1,1),PT%at(1,2),PT%at(1,3),PT%omega)
    PT%at = PT%at / PT%celldm(1)  !  bring at in units of alat
    CALL volume(PT%celldm(1),PT%at(1,1),PT%at(1,2),PT%at(1,3),PT%omega)
    CALL recips(PT%at(1,1),PT%at(1,2),PT%at(1,3),PT%bg(1,1),PT%bg(1,2),PT%bg(1,3))
    !
    ! wat20100225: ATTENTION, for the second order, we had:
    PT%at  = PT%at  * PT%celldm(1)  ! bring at back to bohr units
    PT%tau = PT%tau * PT%celldm(1)  ! also atomic positions in bohr units
    PT%bg  = PT%bg  / PT%celldm(1)  ! reciprocal lattice in 2pi/bohr
    write(*, '(a)') TRIM(file)
    write(*,'("  q vectors:")')
    write(*, '(4x,a,3f12.6)') "p   ", PT%q1(:,1)
    write(*, '(4x,a,3f12.6)') "q   ", PT%q2(:,1)
    write(*, '(4x,a,3f12.6)') "-q-p", PT%q3(:,1)
    write(*,'("  unit cell:")')
    write(*, '(4x,a,3f12.6)') "a1  ", PT%at(:,1)
    write(*, '(4x,a,3f12.6)') "a2  ", PT%at(:,2)
    write(*, '(4x,a,3f12.6)') "a3  ", PT%at(:,3)
    PT%q1   = tpi*PT%q1/PT%celldm(1)  ! q-point in units of 2pi/bohr
    PT%q2   = tpi*PT%q2/PT%celldm(1)  ! q-point in units of 2pi/bohr
    PT%q3   = tpi*PT%q3/PT%celldm(1)  ! q-point in units of 2pi/bohr

  ENDDO
  write(*,'("conditions:")')
  write(*, '(a,3f12.6)') "p (uc)            ", D3UC%q1(:,1)
  write(*, '(a,3f12.6)') "q (uc) - q (sc)   ", D3UC%q2(:,1)-D3SC%q2(:,1)
  write(*, '(a,3f12.6)') "-q-p (uc) + q (sc)", D3UC%q3(:,1)+D3SC%q2(:,1)

  write(*, '("check conditions (q,p from uc; at from sc):")')
  write(*, '(a,3f14.9)') "exp( ii p * at1 ) = ", EXP(ii* SUM(D3UC%q1(:,1) * D3SC%at(:,1)) )
  write(*, '(a,3f14.9)') "exp( ii p * at2 ) = ", EXP(ii* SUM(D3UC%q1(:,1) * D3SC%at(:,2)) )
  write(*, '(a,3f14.9)') "exp( ii p * at3 ) = ", EXP(ii* SUM(D3UC%q1(:,1) * D3SC%at(:,3)) )

  write(*, '(a,3f14.9)') "exp( ii (q_uc - q_sc) * at1 ) = ", EXP(ii* SUM((D3UC%q2(:,1)-D3SC%q2(:,1)) * D3SC%at(:,1)) )
  write(*, '(a,3f14.9)') "exp( ii (q_uc - q_sc) * at2 ) = ", EXP(ii* SUM((D3UC%q2(:,1)-D3SC%q2(:,1)) * D3SC%at(:,2)) )
  write(*, '(a,3f14.9)') "exp( ii (q_uc - q_sc) * at3 ) = ", EXP(ii* SUM((D3UC%q2(:,1)-D3SC%q2(:,1)) * D3SC%at(:,3)) )

  write(*, '(a,3f14.9)') "exp( ii (-q_uc - p + q_sc) * at1 ) = ", EXP(ii* SUM((D3UC%q3(:,1)+D3SC%q2(:,1)) * D3SC%at(:,1)) )
  write(*, '(a,3f14.9)') "exp( ii (-q_uc - p + q_sc) * at2 ) = ", EXP(ii* SUM((D3UC%q3(:,1)+D3SC%q2(:,1)) * D3SC%at(:,2)) )
  write(*, '(a,3f14.9)') "exp( ii (-q_uc - p + q_sc) * at3 ) = ", EXP(ii* SUM((D3UC%q3(:,1)+D3SC%q2(:,1)) * D3SC%at(:,3)) )

  volume_ratio = NINT(D3SC%omega/D3UC%omega)
  WRITE(6, '(2x,a,i3)') "Ratio of volumes, supercell/unitcell:", volume_ratio
  !
  ! wat20100225: Perform checks on the q-points
  ! REMARK: In the old code, only one q-vector is written to the file,
  !         expecting the derivation to be of the form D3(0,-q,q) anyway.
  !         In the new code, we should write down all three q-vectors q1, q2, q3.
  !         Are we doing this ? Are they stored in D3UC%q(:,1-3) ?
  !
  ! wat20100225: How can we make sure the third order dynamical matrix was calculated as
  !              D3(0,-q,q) in the old code ? Is D3SC%nqs==1 sufficient ?
  !
  IF( ABS(SUM(D3SC%q1**2)) > eps8 )  &
       CALL errore(prog, 'Supercell must be calculated as D3(0,q,-q)', 1)

  ! wat20100225: Check if the q points q1, q2, q3 of the unit cell respects the
  !              translational symmetry of the supercell
  !
  IF( ABS(EXP(ii* SUM(D3UC%q1(:,1) * D3SC%at(:,1)) )-1._dp) &
     +ABS(EXP(ii* SUM(D3UC%q1(:,1) * D3SC%at(:,2)) )-1._dp) &
     +ABS(EXP(ii* SUM(D3UC%q1(:,1) * D3SC%at(:,3)) )-1._dp) &
     +ABS(EXP(ii* SUM((D3UC%q2(:,1)-D3SC%q2(:,1)) * D3SC%at(:,1)) )-1._dp) &
     +ABS(EXP(ii* SUM((D3UC%q2(:,1)-D3SC%q2(:,1)) * D3SC%at(:,2)) )-1._dp) &
     +ABS(EXP(ii* SUM((D3UC%q2(:,1)-D3SC%q2(:,1)) * D3SC%at(:,3)) )-1._dp) & 
     +ABS(EXP(ii* SUM((D3UC%q3(:,1)+D3SC%q2(:,1)) * D3SC%at(:,1)) )-1._dp) &
     +ABS(EXP(ii* SUM((D3UC%q3(:,1)+D3SC%q2(:,1)) * D3SC%at(:,2)) )-1._dp) &
     +ABS(EXP(ii* SUM((D3UC%q3(:,1)+D3SC%q2(:,1)) * D3SC%at(:,3)) )-1._dp) > eps8 ) &
  THEN
    CALL errore(prog, 'Supercell does not match the q-points', 2) 
  ENDIF 
  !
  ! Find the translation vectors
  ! wat20100225: This is also a good way to check whether the UC is indeed a unit cell of SC
  ALLOCATE( translations(3, volume_ratio))
  ALLOCATE( atoms_map(D3SC%nat) )
  !
  WRITE(6, '(2x,a)') "Internal SC fractional traslations found (bohr):"
  Do i_c = 1, D3UC%nat
    Do i_sc = 1, D3SC%nat
      !
      delta = D3SC%tau(:,i_sc) - D3UC%tau(:,i_c)
        write(6,'(a,2i3,3f12.6)') "check", i_c, i_sc, delta
      IF(check_int_linearcombination(delta, D3UC%bg)) THEN
      !IF(check_int_linearcombination(delta, D3UC%at)) THEN
        !
        atoms_map(i_sc) = i_c
        !
        index_in_list = check_which_in_list( delta, translations, n_translations_found, D3SC%at)
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
  ALLOCATE(refolded_D3(3,3,3,D3UC%nat,D3UC%nat,D3UC%nat))
  refolded_D3 = 0._dp
  !


  DO at_i = 1,D3SC%nat
     write(*, '(4x, 2i3, 2x, 3f14.10 )') at_i, atoms_map(at_i), D3SC%tau(1,at_i), D3SC%tau(2,at_i), D3SC%tau(3,at_i)
  ENDDO


  DO at_i = 1,D3SC%nat
  DO at_j = 1,D3SC%nat
  DO at_k = 1,D3SC%nat
    c_at_i = atoms_map(at_i)
    c_at_j = atoms_map(at_j)
    c_at_k = atoms_map(at_k)
    DO i = 1,3
    DO j = 1,3
    DO k = 1,3
      ! wat20100225: Okay, these deltas are probably wrong...
      delta1 = D3SC%tau(:,at_i)-D3UC%tau(:,c_at_i)
      delta2 = D3SC%tau(:,at_j)-D3UC%tau(:,c_at_j)
      delta3 = D3SC%tau(:,at_k)-D3UC%tau(:,c_at_k)
      !
      q1r1 = SUM(D3UC%q1(:,1) *delta1) !  + D3UC%p(2,1)*delta1(2)  + D3UC%p(3,1)*delta1(3))
      q2r2 = SUM(D3UC%q2(:,1) *delta2) ! + D3UC%q(2,1)*delta2(2)  + D3UC%q(3,1)*delta2(3))
      q3r3 = SUM(D3UC%q3(:,1)*delta3)! + D3UC%pq(2,1)*delta3(2) + D3UC%pq(3,1)*delta3(3))
      !
!       phase = exp(ii * q1r1) * exp(ii * q2r2) * exp(ii * q3r3)/volume_ratio
      phase = exp(ii * (q1r1+q2r2+q3r3))/volume_ratio
      !
      refolded_D3(i,j,k,c_at_i,c_at_j,c_at_k) = refolded_D3(i,j,k,c_at_i,c_at_j,c_at_k) &
        + D3SC%phiq(i,j,k,at_i,at_j,at_k) * phase
      !
    ENDDO
    ENDDO
    ENDDO
  ENDDO
  ENDDO
  ENDDO

  ! wat20100225: write comparison
  MAXWORK = 0._dp
  MAXDIFF = 0._dp
  
  WRITE(*, '("Third order dynamical matrix:")')
  WRITE(*, '(a)') &
'alp1 alp2 alp3 at1 at2 at3     refolded D3                original D3                diff        diff%'
  DO at_i = 1,D3UC%nat
  DO at_j = 1,D3UC%nat
  DO at_k = 1,D3UC%nat
    DO i = 1,3
    DO j = 1,3
    DO k = 1,3
      !
      if (abs(refolded_D3(i,j,k,at_i,at_j,at_k)-D3UC%phiq(i,j,k,at_i,at_j,at_k)) < 1.d-6) then
        work = 0._dp
      else
        work = 100*(abs(refolded_D3(i,j,k,at_i,at_j,at_k))-abs(D3UC%phiq(i,j,k,at_i,at_j,at_k))) &
                 /(abs(refolded_D3(i,j,k,at_i,at_j,at_k))+abs(D3UC%phiq(i,j,k,at_i,at_j,at_k)))
      endif
      diff = ABS(refolded_D3(i,j,k,at_i,at_j,at_k)-D3UC%phiq(i,j,k,at_i,at_j,at_k))
      IF (work>MAXWORK) THEN
        MAXWORK  = work
        MAXWORKd = diff
      ENDIF
      IF (diff>MAXDIFF) THEN
        MAXDIFF  = diff
        MAXDIFFp = work
      ENDIF
      !
      WRITE(*,'(2(3i3,2x),2x,2(2f13.8,2x),f13.8,f10.5)') i,j,k, &
           at_i,at_j,at_k, &
           refolded_D3(i,j,k,at_i,at_j,at_k), &
           D3UC%phiq(i,j,k,at_i,at_j,at_k), &
           diff, work
      !
    ENDDO
    ENDDO
    ENDDO
  ENDDO
  ENDDO
  ENDDO

  WRITE(*,'(a,f10.1," (",f6.2,")",4x,f6.2," (",f8.1,")")') "maxdiff 10e-6 (%) / maxdiff % (1.d-3):",&
         MAXDIFF*1.d+6, MAXDIFFp, MAXWORK, MAXWORKd*1.d+3
                           


END PROGRAM d3_sc2c


