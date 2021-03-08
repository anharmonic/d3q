!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
MODULE d3matrix_module
!-----------------------------------------------------------------------
CONTAINS
!-----------------------------------------------------------------------
SUBROUTINE d3matrix(d3dyn_in, d3dyn_basename, symmetrize)
  !-----------------------------------------------------------------------
  !
  ! This routine is driver which computes the symmetrized derivative
  ! of the dynamical matrix at q and in the star of q.
  ! The result is written on a iud3dyn file
  !
  USE kinds,           ONLY : DP
  USE ions_base,       ONLY : nat, ityp, ntyp => nsp, tau, atm, amass
  USE cell_base,       ONLY : at, bg, ibrav, celldm
  USE symm_base,       ONLY : s, irt, invs, nsym
  USE lr_symm_base,    ONLY : rtau, nsymq, irgq, minus_q, irotmq
  USE kplus3q,         ONLY : kplusq
  USE d3_basis,        ONLY : patq
  USE io_global,       ONLY : stdout, ionode
  USE qstar_d3_module, ONLY : star_3q
  USE d3_basis,        ONLY : d3_pat2cart, patq, d3_3idx_2_6idx, d3_6idx_2_3idx, &
                              d3_cart2crys, d3_crys2cart
  USE rotate_d3,       ONLY : rotate_and_add_d3
  USE d3_symmetry,     ONLY : d3_symmetrize
  USE d3_shuffle,      ONLY : d3perms, d3_shuffle_global, dperms, iperms
  USE d3_control,      ONLY : print_star, print_perm, print_trev
  USE d3matrix_io2,    ONLY : write_d3dyn_xml2

  IMPLICIT NONE
  COMPLEX(DP),INTENT(in) :: d3dyn_in( 3 * nat, 3 * nat, 3 * nat)
  CHARACTER(len=*),INTENT(in) ::d3dyn_basename

  INTEGER :: isym !, tr_factor
  INTEGER :: dperms_print, nst3q_print
  LOGICAL,OPTIONAl,INTENT(in) :: symmetrize
  INTEGER,PARAMETER :: iud3dyn = 1666
  ! degeneracy of the star of q
  ! index of q in the star of a given sym.op.
  ! index of -q in the star of q (0 if not present)
  ! counter on atoms
  ! counter on atomic type
  ! generic counter
  !
  ! Symmetrizes the dynamical matrix w.r.t. the small group of q
  !
  INTEGER :: istq
  !
  ! For the star of the q triplet:
  INTEGER :: nst3q
  ! n3q : degeneracy of the star of q triplets
  INTEGER :: is3q (48), im3q
  ! is3q : index of q in the star for a given sym
  ! im3q  : index of -q in the star (0 IF not present)
  REAL(DP) :: sx3q(3, 3, 48)
  !
  INTEGER :: i,j,k, iperm, npoint
  !
  COMPLEX(DP),ALLOCATABLE :: d3dyn(:,:,:), d3perm(:,:,:)
  COMPLEX(DP),ALLOCATABLE :: p3dyn(:,:,:, :,:,:), p3perm(:,:,:, :,:,:)
  CHARACTER(len=8),PARAMETER :: sub = "d3matrix"
  !
  ! Prepare internal workspace
  ALLOCATE(d3dyn(3*nat, 3*nat, 3*nat))
  !
  CALL start_clock('d3matrix')
  !
  d3dyn = d3dyn_in
  !
  ! Take the D3 matrix to cartesian axis
  CALL d3_pat2cart (d3dyn, nat, patq(1)%u, patq(2)%u, patq(3)%u)
  !
  ! Symmetrize
  CALL d3_symmetrize(d3dyn, kplusq(1)%xq, kplusq(2)%xq, kplusq(3)%xq, &
                     s, invs, rtau, irt, irgq, at, bg, nsymq, nat, irotmq, minus_q) !, npert_i, npert_f)
  ! NOTE: here d3dyn contains the SYMMETRIZED matrix on CARTESIAN axis
  !
  ! Generates the star of the q triplet
  CALL star_3q(kplusq(1)%xq, kplusq(2)%xq, kplusq(3)%xq, &
               at, bg, nsym, s, invs, nst3q, sx3q, is3q, im3q)
  !
  ! Decide if print the entire star or only the original q1,q2,q3
  nst3q_print = 1
  IF (print_star) nst3q_print = nst3q
  !
  ! Decide if print also the inequivalent permutations
  dperms_print = 1
  IF (print_perm) dperms_print = dperms
  !
  ! Check if time-reversal symmetry will be used
!  tr_factor = 1
!  IF(im3q==0) tr_factor = 2
  !
!  WRITE(stdout, '(9x,a,i2,a)') "Triplet of qs produces",dperms," distinct permutations"
!  WRITE(stdout, '(9x,a,i3,a)') "For each permutation we have", nst3q," triplets in the star"
!  IF(im3q==0) WRITE(stdout, '(5x,a,i3)')   "In addition we have the -q1,-q2,-q3 triplets."
!  WRITE(stdout, '(7x,a,i3)')   "TOTAL NUMBER OF TRIPLETS: ", dperms*nst3q*tr_factor
  !
  ! For each triplet of q's in the star, rotate D3 and write it on file
  ALLOCATE(p3dyn(3,3,3, nat,nat,nat), &
           p3perm(3,3,3, nat,nat,nat))
  ALLOCATE(d3perm(3*nat,3*nat,3*nat))
  !
  npoint = 0
  IF(.not.ionode) RETURN
  !
  STAR_OF_3Q : & ! <<<================================================================
  DO istq = 1, nst3q_print
      !
      ! Repack d3 as (3,3,3, nat,nat,nat) this also makes a copy of d3dyn in p3perm
      CALL d3_3idx_2_6idx(nat, d3dyn, p3perm)
      !
      ! Take to cryst axis
      CALL d3_cart2crys(nat, at, bg, p3perm)
      !
      ! Apply sym.ops.that generate this triplet of q's
      p3dyn = (0._dp, 0._dp)
      DO isym = 1, nsym
        IF( is3q(isym) == istq ) &
          CALL rotate_and_add_d3(p3perm, p3dyn, nat, isym, s, invs, irt, rtau, sx3q(:,:,istq))
        !
      ENDDO
      ! Renormalize:
      p3dyn = p3dyn / nsymq
      !
!       p3dyn = p3perm
      ! Take back to cart axis
      CALL d3_crys2cart(nat, at, bg, p3dyn)
      !
      ! NOTE: from here on p3dyn contain the symmetrized D3 matrix of this triplet
      !       it MUST NOT be changed in the permutations loop!
      !
      PERMS_OF_3Q : & ! <----------------------------------------------------------------
      DO iperm = 1,dperms_print
        npoint = npoint+1
        !
        i = d3perms(iperms(iperm))%i
        j = d3perms(iperms(iperm))%j
        k = d3perms(iperms(iperm))%k
        !
        ! Repack the matrix, shuffle it and repack it again
        CALL d3_6idx_2_3idx(nat, p3dyn, d3perm)
        CALL d3_shuffle_global(nat,  1,2,3, i,j,k, .false., d3perm )
        CALL d3_3idx_2_6idx(nat, d3perm, p3perm)
        !
        ! Report
        WRITE(stdout, '(5x,3(2x,a,i3)," (",3i2," )")') &
                      "point #", npoint, ": star #", istq, "perm #", iperm, i,j,k
!        WRITE(stdout, '(9x,a,3f10.6,a)') "q1: (", sx3q(:,i,istq), " )"
!        WRITE(stdout, '(9x,a,3f10.6,a)') "q2: (", sx3q(:,j,istq), " )"
!        WRITE(stdout, '(9x,a,3f10.6,a)') "q3: (", sx3q(:,k,istq), " )"
        FLUSH( stdout )
  !
  ! In the entire code, for historical reasons, we compute the complex-conjugate of 
  ! the D3 matrix, this is most unfortunate, however we change back to the proper
  ! D3 matrix here, before writing it to file.

        ! Write to file
        CALL write_d3dyn_xml2(d3dyn_basename, sx3q(:,i,istq), sx3q(:,j,istq), sx3q(:,k,istq), &
                             CONJG(p3perm), ntyp, nat, ibrav, celldm, at, ityp, tau, atm, amass)
        !
        ! If -q1,-q2,-q3 is not in the star, use time reversal to compute its D3 and write it to file:
        IF( im3q == 0 .and. print_trev) THEN
          npoint = npoint+1
          WRITE(stdout, '(5x,3(2x,a,i3)," (",3i2," ) ",a)') &
                        "point #", npoint, ": star #", istq, "perm #", iperm, i,j,k, "(time reversal)"
!          WRITE(stdout, '(9x,a,3f10.6,a)') "q1: (", -sx3q(:,i,istq), " )"
!          WRITE(stdout, '(9x,a,3f10.6,a)') "q2: (", -sx3q(:,j,istq), " )"
!          WRITE(stdout, '(9x,a,3f10.6,a)') "q3: (", -sx3q(:,k,istq), " )"
          !
          FLUSH( stdout )
          CALL write_d3dyn_xml2(d3dyn_basename, &
                     -sx3q(:,i,istq), -sx3q(:,j,istq), -sx3q(:,k,istq), &
                     p3perm, ntyp, nat, ibrav, celldm, at, ityp, tau, atm, amass)
        ENDIF
        !
      ENDDO &
      PERMS_OF_3Q     ! ---------------------------------------------------------------->
      !
  ENDDO &
  STAR_OF_3Q     ! ================================================================>>>
  !
!   IF (npoint /= dperms*nst3q*tr_factor) CALL errore(sub, 'Something is wrong with q points...', 1)
  !
  DEALLOCATE(d3dyn, d3perm)
  DEALLOCATE(p3dyn, p3perm)
  !
  CALL stop_clock('d3matrix')
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE d3matrix
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE write_d3dyn_raw (xq1,xq2,xq3, iud3dyn, nat, d3dyn)
  !-----------------------------------------------------------------------
  !
  !     writes in a file the third derivative of dynamical matrix
  USE kinds,       ONLY : DP
  USE io_global,   ONLY : ionode
  !
  IMPLICIT NONE
  !
  REAL(DP)               :: xq1(3), xq2(3), xq3(3)
  COMPLEX(DP),INTENT(in) :: d3dyn (3*nat, 3*nat, 3*nat)
  INTEGER,INTENT(in)     :: iud3dyn, nat
  INTEGER :: nu_i, na_i, i, &
             nu_j, na_j, j, &
             nu_k, na_k, k
  IF ( .NOT. ionode ) RETURN

  WRITE(iud3dyn, '(3f12.6)') xq1
  WRITE(iud3dyn, '(3f12.6)') xq2
  WRITE(iud3dyn, '(3f12.6)') xq3

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
          WRITE(iud3dyn, '(3(3i2,1x),f14.9,2f30.20)') &
                nu_i, na_i, i, nu_j, na_j, j, nu_k, na_k, k, &
                ABS(d3dyn(nu_i, nu_j, nu_k)), d3dyn(nu_i, nu_j, nu_k)
      ENDDO !k
      ENDDO !na_k
      !
    ENDDO !j
    ENDDO !na_j
    WRITE(iud3dyn,'("------------------------------------")')
    !
  ENDDO !i
  ENDDO !na_i

  RETURN

END SUBROUTINE write_d3dyn_raw

!
!-----------------------------------------------------------------------
END MODULE d3matrix_module
!-----------------------------------------------------------------------
