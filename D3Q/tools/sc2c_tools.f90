MODULE sc2c_params
  !
  ! All variables read from file that need dynamical allocation
  !
  USE kinds,     ONLY: DP
  USE constants, ONLY : tpi
  IMPLICIT NONE

  INTEGER,PARAMETER  :: NTYPX     = 10
  REAL(DP),PARAMETER :: eps8      = 1.d-6
  INTEGER, PARAMETER :: seek      = 4
  REAL(DP),PARAMETER :: useek     = 0.25_dp ! 1._dp/DBLE(seek)
  !
  TYPE ph_system_info 
    ! atoms
    INTEGER              :: ntyp
    REAL(DP)             :: amass(ntypx)
    CHARACTER(len=3  )   :: atm(ntypx)
    ! atoms basis
    INTEGER              :: nat
    REAL(DP),ALLOCATABLE :: tau(:,:), zeu(:,:,:)
    INTEGER, ALLOCATABLE :: ityp(:)
    ! unit cell, and reciprocal
    INTEGER              :: ibrav
    CHARACTER(len=9)     :: symm_type
    REAL(DP)             :: celldm(6), at(3,3), bg(3,3)
    REAL(DP)             :: omega
    ! q points
    INTEGER              :: nqs
    REAL(DP)             :: q(3,48)
    ! phonon switches (mostly unused here)
    REAL(DP)             :: epsil(3,3)
    LOGICAL              :: lrigid
    ! dynamical matrix
    COMPLEX(DP), ALLOCATABLE :: phiq(:,:,:,:,:)
    COMPLEX(DP), ALLOCATABLE :: D(:,:)
  END TYPE
  !
CONTAINS
!
!----------------------------------------------------------------------------
SUBROUTINE read_file_sc2c( u, nqs, xq, epsil, lrigid, &
                      ntyp, ityp, nat, tau, zeu, ibrav, symm_type, celldm, at, atm, amass, phiq )
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
!   USE sc2c_params, ONLY: ntypx !, tau, zeu
  !
  IMPLICIT NONE
  !
  ! I/O variables
  INTEGER,INTENT(IN) :: u
  LOGICAL,INTENT(OUT) :: lrigid
  INTEGER,ALLOCATABLE, INTENT(OUT) :: ityp(:)
  REAL(DP),ALLOCATABLE,INTENT(OUT) :: tau(:,:), zeu(:,:,:)
  INTEGER,INTENT(OUT) :: nqs, ntyp, nat, ibrav
  COMPLEX(DP), ALLOCATABLE,INTENT(OUT) :: phiq(:,:,:,:,:) 


  REAL(DP) :: epsil(3,3)
  REAL(DP) :: xq(3,48), celldm(6), at(3,3), amass(ntypx)
  CHARACTER(LEN=3) :: atm(ntypx)
  CHARACTER(LEN=9) :: symm_type
  ! local variables
  INTEGER :: ntyp1,nat1,ibrav1,ityp1
  INTEGER :: i, j, na, nb, nt
  REAL(DP) :: tau1(3), amass1, at1(3,3), celldm1(6), q2
  REAL(DP) :: phir(3),phii(3)
  CHARACTER(LEN=75) :: line
  CHARACTER(LEN=3)  :: atm1
  CHARACTER(LEN=9) symm_type1
  !
  READ(u,*) 
  READ(u,*) 
  !
  ! read cell information from file
  !
  READ(u,*) ntyp,nat,ibrav,(celldm(i),i=1,6)
  if (ibrav==0) then
    read (1,'(a)') symm_type
    read (1,*) ((at(i,j),i=1,3),j=1,3)
  end if

  IF (ntyp.GT.nat) CALL errore('read_f','ntyp.gt.nat!!',ntyp)
  DO nt = 1,ntyp
    READ(u,*) i,atm(nt),amass(nt)
    IF (i.NE.nt) CALL errore('read_f','wrong data read',nt)
  END DO
  ALLOCATE ( ityp(nat), tau(3,nat) )
  DO na=1,nat
    READ(u,*) i,ityp(na),(tau(j,na),j=1,3)
    IF (i.NE.na) CALL errore('read_f','wrong data read',na)
  END DO
  !
  ALLOCATE ( phiq (3,3,nat,nat,48), zeu (3,3,nat) )
  !
  lrigid=.FALSE.
  !
  nqs = 0
100 CONTINUE
  READ(u,*)
  READ(u,'(a)') line
  IF (line(6:14).NE.'Dynamical') THEN
     IF (nqs.EQ.0) CALL errore('read',' stop with nqs=0 !!',1)
     q2 = xq(1,nqs)**2 + xq(2,nqs)**2 + xq(3,nqs)**2
     IF (q2.NE.0.d0) RETURN
     DO WHILE (line(6:15).NE.'Dielectric') 
        READ(u,'(a)',err=200, END=200) line
     END DO
     lrigid=.TRUE.
     READ(u,*) ((epsil(i,j),j=1,3),i=1,3)
     READ(u,*)
     READ(u,*)
     READ(u,*)
     WRITE (*,*) 'macroscopic fields =',lrigid
     WRITE (*,'(3f10.5)') ((epsil(i,j),j=1,3),i=1,3)
     DO na=1,nat
        READ(u,*)
        READ(u,*) ((zeu(i,j,na),j=1,3),i=1,3)
        WRITE (*,*) ' na= ', na
        WRITE (*,'(3f10.5)') ((zeu(i,j,na),j=1,3),i=1,3)
     END DO
     RETURN
200  WRITE (*,*) ' Dielectric Tensor not found'
     lrigid=.FALSE.     
     RETURN
  END IF
  !
  nqs = nqs + 1
  READ(u,*) 
  READ(u,'(a)') line
  READ(line(11:75),*) (xq(i,nqs),i=1,3)
  READ(u,*) 
  !
  DO na=1,nat
     DO nb=1,nat
        READ(u,*) i,j
        IF (i.NE.na) CALL errore('read_f','wrong na read',na)
        IF (j.NE.nb) CALL errore('read_f','wrong nb read',nb)
        DO i=1,3
           READ(u,*) (phir(j),phii(j),j=1,3)
           DO j = 1,3
              phiq (i,j,na,nb,nqs) = CMPLX(phir(j),phii(j),kind=DP)
           END DO
        END DO
     END DO
  END DO
  !
  go to 100
  !
END SUBROUTINE read_file_sc2c

FUNCTION check_int_linearcombination(test_vector, inverted_basis)
  IMPLICIT NONE
  REAL(DP),INTENT(IN) :: test_vector(3)
  REAL(DP),INTENT(IN) :: inverted_basis(3,3)
  REAL(DP) :: coefficients(3)
  LOGICAL :: check_int_linearcombination
  !
  coefficients = MATMUL(transpose(inverted_basis), test_vector)
  !
  check_int_linearcombination = ALL(ABS(coefficients-INT(coefficients+.5_dp))<eps8 )
  !
END FUNCTION check_int_linearcombination

FUNCTION check_which_in_list(test_vector, array_of_vectors, number_of_vectors, translations_basis)
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: test_vector(3)
  INTEGER, INTENT(IN)  :: number_of_vectors
  REAL(DP), INTENT(IN) :: array_of_vectors(3,number_of_vectors)
  REAL(DP), INTENT(IN) :: translations_basis(3,3)
  INTEGER :: i
  INTEGER :: check_which_in_list

  ! i.e it is not in the list
  check_which_in_list = -1
  !
  IF(number_of_vectors==0) RETURN
  !
  DO i = 1,number_of_vectors
    IF( equal_with_translation(test_vector, array_of_vectors(:,i), translations_basis) ) THEN
      check_which_in_list = i
      RETURN
    ENDIF
  ENDDO
END FUNCTION check_which_in_list

FUNCTION equal_with_translation(vector1, vector2, translations_basis)
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: vector1(3), vector2(3)
  REAL(DP), INTENT(IN) :: translations_basis(3,3)
  REAL(DP) :: translation(3)
  INTEGER :: i,j,k
  LOGICAL :: equal_with_translation

  equal_with_translation = .FALSE.
  DO i = -seek, seek
  DO j = -seek, seek
  DO k = -seek, seek
    translation = translations_basis(:,1) * i &
                 +translations_basis(:,2) * j &
                 +translations_basis(:,3) * k
    !
    IF(ALL(ABS(vector1+translation-vector2)<eps8)) THEN
      equal_with_translation = .TRUE.
      RETURN
    ENDIF
  ENDDO
  ENDDO
  ENDDO
END FUNCTION equal_with_translation

!
END MODULE sc2c_params
