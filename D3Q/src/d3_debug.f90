! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE d3_debug

  LOGICAL :: &
    dbg_do_dwfc = .true. ,       &
    dbg_do_dpdvp = .true.
  !
  LOGICAL :: &
    dbg_do_dpdvdp = .true. ,     &
    dbg_do_dpdpdv = .true. ,     &
    dbg_do_drhod2v = .true. ,    &
    dbg_do_rhod3v = .true. ,     &
    dbg_do_ion = .true. ,        &
    dbg_do_smearing = .true. ,   &
      dbg_do_smr_ijk = .true. ,  &
      dbg_do_smr_ij = .true. ,   &
      dbg_do_smr_g = .true. ,    &
    dbg_do_exc = .true. ,        &
      dbg_exc_do_gga = .true.,   &
    dbg_do_nlcc = .true. ,       &
      dbg_do_nlcc_0 = .true. ,   &
      dbg_do_nlcc_123 = .true.
  !
  LOGICAL :: &
    dbg_write_d3_parts = .false., &
    dbg_add_core = .true. 
    
   LOGICAL :: dbg_full_bands = .false.

CONTAINS
!-----------------------------------------------------------------------
SUBROUTINE read_d3_debug(u)
  !-----------------------------------------------------------------------
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  INTEGER,INTENT(in) :: u
  INTEGER :: ios
   CHARACTER(len=512) line
  namelist / d3_debug / &
   dbg_do_dwfc, dbg_do_dpdvp, dbg_do_dpdvdp, dbg_do_dpdpdv, dbg_do_dpdvdp, &
   dbg_do_drhod2v, dbg_do_rhod3v, dbg_do_ion, &
   dbg_do_smearing, dbg_do_smr_ijk, dbg_do_smr_ij, dbg_do_smr_g, &
   dbg_do_exc, dbg_exc_do_gga, &
   dbg_do_nlcc, dbg_do_nlcc_0, dbg_do_nlcc_123, &
   dbg_write_d3_parts, dbg_add_core, &
   dbg_full_bands

    dbg_do_dwfc       = .true.
    dbg_do_dpdvp      = .true.
    dbg_do_dpdvdp     = .true.
    dbg_do_dpdpdv     = .true.
    dbg_do_drhod2v    = .true.
    dbg_do_rhod3v     = .true.
    dbg_do_ion        = .true.
    dbg_do_smearing   = .true.
      dbg_do_smr_ijk  = .true.
      dbg_do_smr_ij   = .true.
      dbg_do_smr_g    = .true.
    dbg_do_exc        = .true.
      dbg_exc_do_gga   = .true.
    dbg_do_nlcc       = .true.
      dbg_do_nlcc_0   = .true.
      dbg_do_nlcc_123 = .true.
   !
   dbg_write_d3_parts = .false.
   dbg_add_core      = .true.
   !
   dbg_full_bands = .false.
   !
   WRITE(stdout,'(5x,a)') "Checking for debug instructions."
   READ (u, d3_debug, iostat = ios)
   !
   IF(ios == 0) THEN
     WRITE(stdout,'(5x,a)') "REMARK: Debug instructions found!"
     WRITE(stdout,'(5x,a)') "_____________ debug flags _____________"
     WRITE(stdout, d3_debug)
     WRITE(stdout,'(5x,a)') "_____________ debug flags _____________"
   ENDIF

   RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE read_d3_debug
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE bcast_d3_debug
  !-----------------------------------------------------------------------
#ifdef __MPI
    USE mp,             ONLY : mp_bcast
    USE io_global,      ONLY : ionode_id
    USE mp_world,   ONLY : world_comm
    !
    IMPLICIT NONE
    !
    CALL mp_bcast( dbg_do_dwfc, ionode_id, world_comm)
    CALL mp_bcast( dbg_do_dpdvp, ionode_id, world_comm)
    !
    CALL mp_bcast( dbg_do_dpdvdp, ionode_id, world_comm)
    CALL mp_bcast( dbg_do_dpdpdv, ionode_id, world_comm)
    CALL mp_bcast( dbg_do_drhod2v, ionode_id, world_comm)
    CALL mp_bcast( dbg_do_rhod3v, ionode_id, world_comm)
    CALL mp_bcast( dbg_do_ion, ionode_id, world_comm)
    CALL mp_bcast( dbg_do_smearing, ionode_id, world_comm)
      CALL mp_bcast( dbg_do_smr_ijk, ionode_id, world_comm)
      CALL mp_bcast( dbg_do_smr_ij, ionode_id, world_comm)
      CALL mp_bcast( dbg_do_smr_g, ionode_id, world_comm)
    CALL mp_bcast( dbg_do_exc, ionode_id, world_comm)
      CALL mp_bcast( dbg_exc_do_gga, ionode_id, world_comm)
    CALL mp_bcast( dbg_do_nlcc, ionode_id, world_comm)
      CALL mp_bcast( dbg_do_nlcc_0, ionode_id, world_comm)
      CALL mp_bcast( dbg_do_nlcc_123, ionode_id, world_comm)
   !
   CALL mp_bcast( dbg_write_d3_parts, ionode_id, world_comm)
   CALL mp_bcast( dbg_add_core, ionode_id, world_comm)
   CALL mp_bcast( dbg_full_bands, ionode_id, world_comm)
#endif
    !   
  !-----------------------------------------------------------------------
END SUBROUTINE bcast_d3_debug
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE dbgwrite_d3dyn (d3dyn_x, filename, isw)
  !-----------------------------------------------------------------------
  !
  !     writes in a file the third derivative of dynamical matrix
  !     isw = +1  :  d3dyn_x is in cartesian axis
  !     isw = -1  :  rotates d3dyn_x from the basis of pattern to
  !                      cartesian axis
  USE kinds,           ONLY : DP
  USE ions_base,       ONLY : nat
  USE io_global,       ONLY : ionode
  USE d3_basis,        ONLY : patq
  USE kplus3q,         ONLY : kplusq
  !
  USE ions_base,       ONLY : nat
  USE cell_base,       ONLY : at, bg
  USE symm_base,       ONLY : s, irt, invs
  USE lr_symm_base,    ONLY : rtau,irgq,irotmq,nsymq,minus_q
  USE d3_basis,        ONLY : d3_pat2cart, patq
  USE d3_symmetry,     ONLY : d3_symmetrize
  USE d3_shuffle,      ONLY : d3_shuffle_global
  USE d3matrix_module, ONLY : d3matrix
  USE d3com,           ONLY : fild3dyn
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
  !RETURN   ! *
  ! **********
  IF(.not.dbg_write_d3_parts) RETURN

  CALL d3matrix(d3dyn_x, TRIM(fild3dyn)//"."//TRIM(filename))

  ! **********
  RETURN   ! *
  ! **********

  !Also write in human easily readable format:
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
    CALL d3_shuffle_global(nat,  1,2,3, i1,i2,i3, .false., d3dyn_shuffled )
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
END SUBROUTINE dbgwrite_d3dyn

END MODULE d3_debug
