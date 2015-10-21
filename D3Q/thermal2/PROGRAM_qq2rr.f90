!
! programmino di esempio...
!
PROGRAM qq2rr
  USE kinds,       ONLY : DP
  USE d3matrix_io, ONLY : read_d3dyn_xml
  USE d3_shuffle,      ONLY : nperms, d3perms_order
  USE parameters,  ONLY : ntypx
  USE input_fc,    ONLY : ph_system_info

  IMPLICIT NONE
  REAL(DP) :: xq(3,3)
  COMPLEX(DP),ALLOCATABLE :: p3(:,:,:, :,:,:), q3(:,:,:, :,:,:)
  INTEGER,ALLOCATABLE :: ityp(:)
  REAL(DP),ALLOCATABLE :: tau(:,:)
  CHARACTER(len=256) :: fname1, title
  INTEGER :: nat, i,j,k,a,b,c, ios, ntyp, nt, na, icar, ic,jc, iperm
  LOGICAL :: found, first
  CHARACTER(LEN=3) :: atm(ntypx)
  INTEGER :: nq(3)
  !
  TYPE(ph_system_info) :: S
  !
  first=.true.
  nq = (/ 4,4,4 /)
  fname1 = "anh"
  !
  DO k = 0, nq(1)-1
  DO j = 0, nq(2)-1
  DO i = 0, nq(3)-1
    !
    IF(first)THEN
      xq = 0._dp
      first=.false.
      CALL read_d3dyn_xml(fname1, xq(:,1), xq(:,2), xq(:,3), d3=p3, nat=S%nat, atm=S%atm, ntyp=S%ntyp, &
                          ityp=S%ityp, ibrav=S%ibrav, celldm=S%celldm, at=S%at, amass=S%amass,&
                          tau=S%tau, seek=.true.)

      CALL latgen( S%ibrav, S%celldm, S%at(:,1), S%at(:,2), S%at(:,3), S%omega )
      S%at=S%at/S%celldm(1)
      CALL recips(S%at(:,1), S%at(:,2), S%at(:,3), S%bg(:,1), S%bg(:,2), S%bg(:,3))
      !
    ELSE
      xq(:,2) = (/ DBLE(i)/nq(1), DBLE(j)/nq(2), DBLE(k)/nq(3) /)
      xq(:,3) = -xq(:,1)-xq(:,2)
      CALL cryst_to_cart (3,xq,S%bg,+1)
      WRITE(*,'(3(3f10.4,3x))') xq
      found = .false.
      iperm = 0
      DO WHILE(.not. found)
        iperm=iperm+1
        IF(iperm>nperms) STOP 777
        a = d3perms_order(1,iperm)
        b = d3perms_order(2,iperm)
        c = d3perms_order(3,iperm)
        WRITE(*,'(5x,3(3f10.4,3x))') xq(:,a), xq(:,b), xq(:,c)
        CALL read_d3dyn_xml(fname1, xq(:,a), xq(:,b), xq(:,c), at=S%at, d3=p3, seek=.true., found=found)
      ENDDO
    ENDIF
    !
  END DO 
  END DO 
  END DO 
 !
END PROGRAM qq2rr


