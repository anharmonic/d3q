!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Compute a D3 matrix from finite differences of D2 matrices and compares this
! with a D3 matrix from input. 
!
! MUST BE REFACTORED BECAUSE IT DOES NOT READ THE CURRENT FILE FORMAT!!
!
!----------------------------------------------------------------------------
PROGRAM sc3c
  !----------------------------------------------------------------------------
  !  The name and order of files is not important as long as q=0 is the first
  !
  USE d3_tools
  USE sc2c_params, ONLY : read_file_sc2c, ph_system_info
  USE mp,         ONLY : mp_start, mp_end, mp_barrier
  USE mp_global,  ONLY : nproc, mpime
  USE fft_scalar, ONLY : cfft3d
  !USE iotk_module
  !
  IMPLICIT NONE
  !
  CHARACTER(len=4),PARAMETER :: prog = 'sc3c'
  CHARACTER(len=256) :: fildisp1, fildisp2, fildynSC !filin, filj, filf, flfrc
  !
  TYPE(ph_system_info),TARGET  :: DM1, DM2
  TYPE(ph_system_info),POINTER :: PT
  TYPE(d3_system_info)         :: D3

  INTEGER,ALLOCATABLE :: atoms_map(:)
  INTEGER :: ierr, i_c, i_sc, n_translations_found = 0
  INTEGER,PARAMETER :: u = 1
  INTEGER :: index_in_list
  REAL(DP),ALLOCATABLE :: translations(:,:), dtau(:,:)
  REAL(DP) :: delta(3), q(3), qdr, work, h
  INTEGER  :: volume_ratio
  INTEGER :: i,j,k,l,m,n, at_i, at_j, c_at_i, c_at_j, nat
  COMPLEX(DP),ALLOCATABLE :: refolded_D(:,:,:,:,:,:), finite_diff(:,:,:,:), DD(:,:)
  COMPLEX(DP),PARAMETER :: ii = (0._dp, 1._dp)
  COMPLEX(DP) :: zwork
  !
  NAMELIST / input / fildisp1, fildisp2, fildynSC
  CALL input_from_file ( ) 
  READ ( 5, input )

  

  ! Read the two second order dynamical matrices for the finite difference
  DO i = 1,2
    IF( i==1 ) THEN
      PT => DM1
      OPEN (unit=u, file=TRIM(fildisp1),   status='old', form='formatted', iostat=ierr)
      CALL errore(prog, 'Cannot open '//TRIM(fildisp1), ierr)
      WRITE(6,'(2x,a,i3)') 'Reading second order dynamical matrix from file "'//TRIM(fildisp1)//'"', i
    ELSE
      PT => DM2
      OPEN (unit=u, file=TRIM(fildisp2), status='old', form='formatted', iostat=ierr)
      CALL errore(prog, 'Cannot open '//TRIM(fildisp2), ierr)
      WRITE(6,'(2x,a,i3)') 'Reading second order dynamical matrix from file "'//TRIM(fildisp2)//'"', i
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
  !
  ! Read the third order dynamical matrix
  OPEN (unit=u, file=TRIM(fildynSC),   status='old', form='formatted', iostat=ierr)
  CALL errore(prog, 'Cannot open '//TRIM(fildynSC), ierr)
  WRITE(6,'(2x,a)') 'Reading third order dynamical matrix from file "'//TRIM(fildynSC)//'"'

  CALL read_d3(u, D3%nat, D3%ntyp, D3%ityp, D3%atm, D3%amass, D3%tau, D3%ibrav, D3%at, D3%celldm, &
                     D3%q1,D3%q2,D3%q3, D3%phiq)
  write(*, '(a,3f12.6)') "p   ", D3%q1(:,1)
  write(*, '(a,3f12.6)') "q   ", D3%q2(:,1)
  write(*, '(a,3f12.6)') "-q-p", D3%q3(:,1)

  CALL latgen(D3%ibrav,D3%celldm,D3%at(1,1),D3%at(1,2),D3%at(1,3),D3%omega)
    D3%at = D3%at / D3%celldm(1)  !  bring at in units of alat
    !
  CALL volume(D3%celldm(1),D3%at(1,1),D3%at(1,2),D3%at(1,3),D3%omega)
  CALL recips(D3%at(1,1),D3%at(1,2),D3%at(1,3),D3%bg(1,1),D3%bg(1,2),D3%bg(1,3))
  !
  ALLOCATE(dtau(3,DM1%nat))
  dtau = DM1%tau - DM2%tau
  h    = SQRT(SUM(dtau**2)) !/DM1%celldm(1)
  Write(*,'(2x,"Displacement (bohr, alat)",1f12.6,",",1e10.1)') h, h/DM1%celldm(1)

  nat = DM1%nat
  ALLOCATE(DM1%D(3*nat,3*nat))
  open(unit=1077, file='d2ionq.xml+', action='read', status='old')
!   CALL iotk_scan_dat(1077, 'd2ionq', DM1%D)
  do i = 1,3*nat
  do j = 1,3*nat
    read(1077, '(2i4,2f32.16)') l,m, zwork
    DM1%D(l,m) = zwork
  enddo
  enddo

  close(1077)
  !
  ALLOCATE(DM2%D(3*nat,3*nat))
  open(unit=1077, file='d2ionq.xml-', action='read', status='old')
!   CALL iotk_scan_dat(1077, 'd2ionq', DM2%D)
  do i = 1,3*nat
  do j = 1,3*nat
    read(1077, '(2i4,2f32.16)') l,m, zwork
    DM2%D(l,m) = zwork
  enddo
  enddo
  close(1077)
  !
  ALLOCATE(D3%Dn(3*nat,3*nat,3*nat))
  open(unit=1077, file='d3ionq-n.xml', action='read', status='old')
!   CALL iotk_scan_dat(1077, 'd3ionq', D3%D)
  do i = 1,3*nat
  do j = 1,3*nat
  do k = 1,3*nat
    read(1077, '(3i4,2f32.16)') l,m,n, zwork
    D3%Dn(l,m,n) = zwork
  enddo
  enddo
  enddo
  close(1077)
  ALLOCATE(D3%Do(3*nat,3*nat,3*nat))
  open(unit=1077, file='d3ionq-n.xml', action='read', status='old')
!   CALL iotk_scan_dat(1077, 'd3ionq', D3%D)
  do i = 1,3*nat
  do j = 1,3*nat
  do k = 1,3*nat
    read(1077, '(3i4,2f32.16)') l,m,n, zwork
    D3%Do(l,m,n) = zwork
  enddo
  enddo
  enddo
  close(1077)
  !
  ALLOCATE(DD(3*nat,3*nat))
  DD=(DM2%D(:,:) - DM1%D(:,:))/h

  ALLOCATE(finite_diff(3,3,DM1%nat,DM1%nat))
  finite_diff=(DM2%phiq(:,:,:,:,1) - DM1%phiq(:,:,:,:,1))/h


  ! WHAT WOULD BE THE CORRESPONDING TEST ?
  !IF( ANY(ABS(SG%q(:,1)) > eps8 ) ) &
  !  CALL errore(prog, 'Supercell must be at Gamma', 1)

  volume_ratio = INT(DM1%omega/DM2%omega +0.5_dp)
  WRITE(6, '(2x,a,i3)') "Ratio of volumes, unit cells of second order matrices:", volume_ratio
  volume_ratio = INT(D3%omega/DM1%omega +0.5_dp)
  WRITE(6, '(2x,a,i3)') "Ratio of volumes, supercell/unitcell:", volume_ratio


  !
#define _A 3
#define _N 2
#define _INDEX_ _A,i,j,_N,at_i,at_j
! #define _INDEX_ i,j,_A,at_i,at_j,_N
  WRITE(*, '("Dynamical matrix:")')
  WRITE(*, '(a)') &
'alp1 alp2 at1 at2     finite diff               original D                diff        diff%'
  DO at_i = 1,DM1%nat
  DO at_j = 1,DM2%nat
    DO i = 1,3
    DO j = 1,3
      !
!       work = ABS(finite_diff(i,j,at_i,at_j)-CQ%phiq(i,j,at_i,at_j,1))
!       work = work /(ABS(refolded_D(i,j,at_i,at_j))+ABS(CQ%phiq(i,j,at_i,at_j,1)))
      work = 100*(1._dp-finite_diff(i,j,at_i,at_j)/D3%phiq(_INDEX_))
      WRITE(*,'(2(2i3,2x),2x,2(2f12.6,2x),f12.6,f12.6)') i,j,at_i,at_j, &
                                                         finite_diff(i,j,at_i,at_j), D3%phiq(_INDEX_), &
                                                         ABS(finite_diff(i,j,at_i,at_j)-D3%phiq(_INDEX_)), &
                                                         work
      !
    ENDDO
    ENDDO
  ENDDO
  ENDDO

#define _M 6
#define _INDEX2_ _M,i,j
  WRITE(*, '("Dynamical matrix (NEW ALGORITHM!!!!):")')
  WRITE(*, '(a)') &
'alp1 alp2     finit diff                original D                diff        diff%'
  DO i = 1,3*nat
  DO j = 1,3*nat
      work = 100*(1._dp-DD(i,j)/D3%Dn(_INDEX2_))
      WRITE(*,'((2i3,2x),2x,2(2f12.6,2x),f12.6,f12.6)') i,j,DD(i,j), D3%Dn(_INDEX2_), &
                                                        ABS(DD(i,j)-D3%Dn(_INDEX2_)), &
                                                        work
      !
  ENDDO
  ENDDO
  WRITE(*, '("Dynamical matrix (ORIGINAL ALGORITHM):")')
  WRITE(*, '(a)') &
'alp1 alp2     finit diff                original D                diff        diff%'
  DO i = 1,3*nat
  DO j = 1,3*nat
      work = 100*(1._dp-DD(i,j)/D3%Do(_INDEX2_))
      WRITE(*,'((2i3,2x),2x,2(2f12.6,2x),f12.6,f12.6)') i,j,DD(i,j), D3%Do(_INDEX2_), &
                                                        ABS(DD(i,j)-D3%Do(_INDEX2_)), &
                                                        work
      !
  ENDDO
  ENDDO


END PROGRAM sc3c


