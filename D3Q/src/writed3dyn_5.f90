!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE writed3dyn_5 (d3dyn_x, filename, isw)
  !-----------------------------------------------------------------------
  !
  !     writes in a file the third derivative of dynamical matrix
  !     isw = +1  :  d3dyn_x is in cartesian axis
  !     isw = -1  :  rotates d3dyn_x from the basis of pattern to
  !                      cartesian axis
  USE kinds,       ONLY : DP
  USE ions_base,   ONLY : nat
  USE io_global,   ONLY : ionode
  USE d3_basis,    ONLY : patq
  USE kplus3q,     ONLY : kplusq
  !
  USE ions_base,   ONLY : nat
  USE cell_base,   ONLY : at, bg
  USE symm_base,   ONLY : s, irt, invs
  USE phcom,       ONLY : rtau,irgq,irotmq,nsymq,minus_q
  USE d3_basis,    ONLY : d3_pat2cart, patq
  USE d3_symmetry, ONLY : d3_symmetrize
  USE d3_shuffle,  ONLY : d3_shuffle_global
  USE d3matrix_module,    ONLY : d3matrix
  USE d3com,       ONLY : fild3dyn
  !
  IMPLICIT NONE
  !
  COMPLEX(DP),INTENT(in)      :: d3dyn_x (3 * nat, 3 * nat, 3 * nat)   ! the third derivative of the dynamical matrix
  CHARACTER(len=*),INTENT(in) :: filename ! input: the name of the file
  INTEGER,INTENT(in) :: isw ! input: switch
  !
  INTEGER :: iud3dyn
       !na_j, na_k
  INTEGER :: nu_i, na_i, i, &
             nu_j, na_j, j, &
             nu_k, na_k, k
  COMPLEX (DP), ALLOCATABLE :: aux (:,:,:), d3dyn_shuffled(:,:,:)
  ! auxiliary space
  CHARACTER(len=256) :: order
  INTEGER :: i1,i2,i3, ios

  ! **********
  RETURN   ! *
  ! **********

  CALL d3matrix(d3dyn_x, TRIM(fild3dyn)//"_"//TRIM(filename))


  IF ( .NOT. ionode ) RETURN

  ALLOCATE  (aux( 3 * nat, 3 * nat, 3 * nat))
  ALLOCATE  (d3dyn_shuffled( 3 * nat, 3 * nat, 3 * nat))
  !
  d3dyn_shuffled = d3dyn_x

!   IF(nsymq>1)THEN
  ! Take the D3 matrix to cartesian axes

  CALL d3_pat2cart(d3dyn_shuffled, nat, patq(1)%u, patq(2)%u, patq(3)%u)

!     WRITE(*,"(7x,'< Symmetrizing D3 matrix and rotating to cartesian axis.')")
 CALL d3_symmetrize(d3dyn_shuffled, kplusq(1)%xq, kplusq(2)%xq, kplusq(3)%xq, &
                     s, invs, rtau, irt, irgq, at, bg, nsymq, nat, irotmq, minus_q) !, npert_i, npert_f)
!   ENDIF

  CALL getenv('D3ORDER',order)
  READ(order,*,iostat=ios) i1,i2,i3
  IF(ios==0.and.(i1/=1.and.i2/=2.and.i3/=3)) THEN
    WRITE(*,"(7x,'< Reordering D3 matrix from',3i2,' to 1 2 3')") i1,i2,i3
    CALL d3_shuffle_global( 1,2,3, i1,i2,i3, .false., d3dyn_shuffled )
  ENDIF
  !
  aux = d3dyn_shuffled

  iud3dyn = 57

  OPEN (unit = iud3dyn, file = TRIM(filename), status = 'unknown')
  DO na_i = 1,nat
  DO i    = 1,3
  nu_i = i + (na_i-1)*3
!     print*, na_i,i,nu_i
    !
    DO na_j = 1,nat
    DO j    = 1,3
    nu_j = j + (na_j-1)*3
      !
      DO na_k = 1,nat
      DO k    = 1,3
      nu_k = k + (na_k-1)*3
          write(iud3dyn, '(3(3i2,1x),f14.9,2f30.20)') &
                nu_i, na_i, i, nu_j, na_j, j, nu_k, na_k, k, &
                ABS(aux(nu_i, nu_j, nu_k)), aux(nu_i, nu_j, nu_k)
      ENDDO !k
      ENDDO !na_k
      !
    ENDDO !j
    ENDDO !na_j
    write(iud3dyn,'("------------------------------------")')
    !
  ENDDO !i
  ENDDO !na_i
  CLOSE (iud3dyn)

  DEALLOCATE (aux, d3dyn_shuffled)
  RETURN
END SUBROUTINE writed3dyn_5
