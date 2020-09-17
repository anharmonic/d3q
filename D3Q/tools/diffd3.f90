!
! Small example program, reads two D3 matrix files and print out the difference
!
PROGRAM read3
  USE kinds,        ONLY : DP
  USE d3matrix_io2, ONLY : read_d3dyn_xml2
  use cell_base,    ONLY : at, ibrav, celldm, omega


  IMPLICIT NONE
  REAL(DP) :: xq1(3), xq2(3), xq3(3), maxdiff, maxperc, perc, diff
  COMPLEX(DP),ALLOCATABLE :: p3(:,:,:, :,:,:), q3(:,:,:, :,:,:)
  CHARACTER(len=256) :: fname1, fname2, basename
  CHARACTER(len=5) :: format_version
  INTEGER :: natoms, i,j,k,a,b,c
  LOGICAL :: found, exst
  !
  ! Metodo 1: fornisco fname='anh' e i tre vettori q per trovare il file:
  !READ(*,*) fname1
  !READ(*,*) fname2

  INTEGER,INTRINSIC :: iargc
  INTEGER :: nargs
  !
  nargs = iargc()
  !
  IF(nargs>=1) THEN
    CALL getarg(1, fname1)
  ELSE
    PRINT*, "enter file1"
    READ(*,"(a256)") fname1
  ENDIF
  IF(nargs>=2) THEN
    CALL getarg(2, fname2)
  ELSE
    PRINT*, "enter file2 or directory"
    READ(*,"(a256)") fname2
  ENDIF
  
  IF(nargs>=3) THEN
    CALL getarg(3, basename)
    fname1 = TRIM(fname1)//"/"//TRIM(basename)
    fname2 = TRIM(fname2)//"/"//TRIM(basename)
  ENDIF
  !
  !CALL read_d3dyn_xml2(fname1, d3=p3, nat=natoms)
  !CALL read_d3dyn_xml2(fname2, d3=q3)
  INQUIRE(file=fname1, exist=exst)
  IF(.not. exst) STOP 1
  CALL read_d3dyn_xml2(fname1, xq1, xq2, xq3, d3=p3, nat=natoms, ibrav=ibrav, &
                      celldm=celldm, at=at, file_format_version=format_version)
  IF(format_version=="1.0.0") p3 = CONJG(p3)

!  write(0,'(3(3f8.4,2x))') xq1,xq2,xq3

  CALL latgen( ibrav, celldm, at(:,1), at(:,2), at(:,3), omega )
  at=at/celldm(1)


  INQUIRE(file=fname2, exist=exst)
  IF(.not. exst) STOP 2
  CALL read_d3dyn_xml2(fname2, xq1, xq2, xq3, d3=q3, file_format_version=format_version)
  IF(format_version=="1.0.0") q3 = CONJG(q3)
  !
  maxdiff = 0._dp
  maxperc = 0._dp
  !
  write(*,'(2a)') "p3: ", trim(fname1)
  write(*,'(2a)') "q3: ", trim(fname2)
  write(*,'(2a)') "i,j,k,a,b,c, i+3*(a-1),j+3*(b-1),k+3*(c-1), diff,",&
                  " ABS(p3(i,j,k,a,b,c)), ABS(q3(i,j,k,a,b,c)),  p3(i,j,k,a,b,c), q3(i,j,k,a,b,c), perc, maxperc"

  DO a = 1,natoms
  DO b = 1,natoms
  DO c = 1,natoms
    DO i = 1,3
    DO j = 1,3
    DO k = 1,3
        diff = ABS(p3(i,j,k,a,b,c)-q3(i,j,k,a,b,c))
        maxdiff = MAX(maxdiff,ABS(diff))
        IF(diff>1.d-5)THEN
            perc = 100*abs(p3(i,j,k,a,b,c)-q3(i,j,k,a,b,c))/(abs(p3(i,j,k,a,b,c))+abs(q3(i,j,k,a,b,c)))
            maxperc = MAX(maxperc,ABS(perc))
        ELSE
            perc=0._dp
        ENDIF
        !
        write(*, '(3(3i2,1x),3f14.5,2(2f22.7,3x),2(f10.4))') &
              i,j,k,a,b,c, i+3*(a-1), j+3*(b-1), k+3*(c-1),&
              diff, &
              ABS(p3(i,j,k,a,b,c)), ABS(q3(i,j,k,a,b,c)), &
              p3(i,j,k,a,b,c), q3(i,j,k,a,b,c), perc, maxperc

    ENDDO
    ENDDO
    ENDDO
  ENDDO
  ENDDO
  ENDDO
!  WRITE(*,'(a,f12.6,f10.4)') "maxdiff:", maxdiff*1.d+6, maxperc
  WRITE(*,'(a,f15.9,f10.4)') "maxdiff:", maxdiff, maxperc

!   fname='anh'
!   !
!   ! NOTA!!!!!!!!!!!!!
!   !  per il momento questo metodo funziona SOLO SE "at" del modulo cell_base contiene
!   !  già gli assi cartesiani della cella unitaria... è un casino ma per il momento non posso
!   !  fare altrimenti (MEMO: sistemare)
!   xq1 = (/ 0._dp,  0._dp,  0._dp /)
!   xq2 = (/ 0._dp,  0._dp, -1._dp /)
!   xq3 = (/ 0._dp,  0._dp,  1._dp /)
!   !
!   at(:,1) = (/ -0.5,  0.0,  0.5 /)
!   at(:,2) = (/  0.0,  0.5,  0.5 /)
!   at(:,3) = (/ -0.5,  0.5,  0.0 /)
!   !
!   ! in questo caso l'opzione seek=.true. DEVE essere specificata
!   CALL read_d3dyn_xml2(fname, xq1, xq2, xq3, nat=natoms, seek=.true.)
!   ! leggo solo nat (nota che devo specificare nat=argomento)
!   WRITE(*,'(5x,a,i3)') "File found, nat:", natoms
!   ! inoltre in output ricevo xq1, xq2 e xq3 in coordinate cartesiane:
!   !
!   ! Metodo 2: fornisco il nome completo del file
!   fname="anh_Q1.-1o2_-1o2_0_Q2.1o2_1o2_0_Q3.0_0_0"
!   ! Nota: non serve allocare p3, viene allocato automaticamente alla dimensione giusta all'interno
!   CALL read_d3dyn_xml2(fname, xq1, xq2, xq3, p3) 
!   ! ... in questo caso xq1, xq2 e xq3 vengono letti
!   WRITE(*,'(5x,a,3(3f10.5,2x))') "File found, :", xq1, xq2, xq3
!   !
!   ! Metodo 3: se il file non esiste e ho specificato il parametro "found" il codice non crasha:
!   CALL read_d3dyn_xml2("brfgzap", xq1, xq2, xq3, found=found) 
!   ! .. questo funziona sia con seek=.true. che con seek assente o falso
!   WRITE(*, '(5x,a,l8)') "File found?", found
!   !
!   ! In generale:
!   ! CALL read_d3dyn_xml2(basename, xq1,xq2,xq3, d3, ntyp, nat, ibrav, celldm, at, ityp, tau, atm, amass,found,seek)
!   !
!   ! I parametri sono:
! !     INTEGER,OPTIONAL,INTENT(out)          :: nat, ntyp       ! number of atoms and atomic species
! !     INTEGER,OPTIONAL,INTENT(out)          :: ibrav           ! lattice type
! !     !
! !     REAL(DP),OPTIONAL,INTENT(out)             :: celldm(6)   ! cell parameters depending on ibrav
! !     REAL(DP),OPTIONAL,INTENT(out)             :: at(3,3)     ! unit cell vectors NOTE: only filled if ibrav=0
! !     REAL(DP),ALLOCATABLE,OPTIONAL,INTENT(out) :: tau(:,:)    ! (3,nat) atomic positions, cartesians and alat
! !     REAL(DP),OPTIONAL,INTENT(out)             :: amass(ntypx)! (ntyp) mass of ions
! !     INTEGER,ALLOCATABLE,OPTIONAL,INTENT(out)  :: ityp(:)     ! (nat)  index of atomic types 
! !     CHARACTER(len=3),OPTIONAL,INTENT(out)     :: atm(ntypx)  ! (ntyp) atomic labels (es. Si)
!   !
!   ! ntypx viene dal modulo parameters (USE parameters, ONLY : ntypx), mentre ityp e tau sono allocabili (così è nel resto del codice)
!   !
!   ! NOTA: at (la cella del reticolo) viene letto SOLO SE ibrav=0, altrimenti devi generartelo a mano:
!   !  IF(ibrav/=0) CALL latgen( ibrav, celldm, at(:,1), at(:,2), at(:,3), omega )
!   !
  !
END PROGRAM read3


