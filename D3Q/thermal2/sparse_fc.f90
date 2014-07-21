!
! Copyright Lorenzo Paulatto, 2014 - released under the CeCILL licence v 2.1 
!   <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
MODULE sparse_fc
  USE kinds,            ONLY : DP
  USE parameters,       ONLY : ntypx
  USE io_global,        ONLY : stdout
  !
  TYPE forceconst3_sparse_helper
    REAL(DP),ALLOCATABLE :: fc(:)    ! the matrix elements
    INTEGER,ALLOCATABLE  :: idx(:,:) ! the index 3*nat,3*nat,3*nat
  END TYPE forceconst3_sparse_helper
  !
  TYPE forceconst3_sparse
    ! q points
    INTEGER :: n_R = 0
    INTEGER :: i_00 = -1 ! set to the index where both R2 and R3 are zero
    INTEGER :: nq(3) ! initial q-point grid size (unused)
    INTEGER,ALLOCATABLE  :: yR2(:,:), yR3(:,:) ! crystalline coords  3*n_R
    REAL(DP),ALLOCATABLE :: xR2(:,:), xR3(:,:) ! cartesian coords    3*n_R
    INTEGER,ALLOCATABLE  :: n_terms(:) ! number of non-zero elements for each R2,R3
    TYPE(forceconst3_sparse_helper),ALLOCATABLE :: dat(:)
    !
  END TYPE forceconst3_sparse
  !
!   TYPE forceconst3_grid
!     ! q points
!     INTEGER :: n_R = 0
!     INTEGER :: i_00 = -1 ! set to the index where both R2 and R3 are zero
!     INTEGER :: nq(3) ! initial q-point grid size (unused)
!     INTEGER,ALLOCATABLE  :: yR2(:,:), yR3(:,:) ! crystalline coords  3*n_R
!     REAL(DP),ALLOCATABLE :: xR2(:,:), xR3(:,:) ! cartesian coords    3*n_R
!     REAL(DP),ALLOCATABLE :: FC(:,:,:,:) ! 3*nat,3*nat,3*nat, n_R
!     !TYPE(index_r) :: idx2, idx3
!   END TYPE forceconst3_grid  
  CONTAINS
  !
  SUBROUTINE fc3_all_to_sparse(fc, sfc)
    !
    USE input_fc ONLY : forceconst3_grid
    IMPLICIT NONE
    !
    TYPE(forceconst3_grid),INTENT(in)    :: fc
    TYPE(forceconst3_sparse),INTENT(out) :: sfc
    !
    REAL(DP),PARAMETER :: eps = 1.d-12
    INTEGER :: iR
    INTEGER :: na1,  na2,  na3,  j1,  j2,  j3, jn1, jn2, jn3
    !
    sfc%n_R =fc%n_R
    sfc%i_00 =fc%i_00
    sfc%nq(3) =fc%nq(3)
    !
    ALLOCATE(sfc%yR2(3,sfc%n_R), sfc%yR3(3,sfc%n_R))
    ALLOCATE(sfc%xR2(3,sfc%n_R), sfc%xR3(3,sfc%n_R))
    sfc%yR2 = fc%yR2;    sfc%yR3 = fc%yR3
    sfc%xR2 = fc%xR2;    sfc%xR3 = fc%xR3
    ALLOCATE(sfc%dat(sfc%n_R))
    ALLOCATE(sfc%n_terms(sfc%n_R))
    
    DO iR = 1, sfc%n_R
      sfc%n_terms(iR) = COUNT(  ABS(fc%FC(:,:,:,iR))>eps )
      ALLOCATE(sfc%dat(iR)%idx(3,sfc%n_terms(iR)))
      ALLOCATE(sfc%dat(iR)%fc(sfc%n_terms(iR)))
    ENDDO
    
    DO iR = 1, sfc%n_R
      !
      n_found = 0
      !
      DO na3=1,S%nat 
      DO na2=1,S%nat 
      DO na1=1,S%nat
        DO j3=1,3     
        jn3 = j3 + (na3-1)*3
        DO j2=1,3     
        jn2 = j2 + (na2-1)*3
        DO j1=1,3
        jn1 = j1 + (na1-1)*3
          !
          IF(ABS(fc%FC(jn1,jn2,jn3,iR))>eps)THEN
            n_found = n_found + 1
            IF(n_found > sfc%n_terms(iR)) CALL errore("fc3_all_to_sparse", "too many terms", iR)
            sfc%dat(iR)%idx(1,n_found) = jn1
            sfc%dat(iR)%idx(2,n_found) = jn2
            sfc%dat(iR)%idx(3,n_found) = jn3
            sfc%dat(iR)%FC(n_found) = fc%FC(jn1,jn2,jn3,iR)
          ENDIF
          !
        ENDDO
        ENDDO
        ENDDO
      ENDDO
      ENDDO
      ENDDO
      !
      IF(n_found /= sfc%n_terms(iR)) CALL errore("fc3_all_to_sparse", "wrong number of terms", iR)
      !
    ENDDO
    !
  END SUBROUTINE fc3_all_to_sparse

END MODULE
