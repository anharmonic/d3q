!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE d3_efermi_shift

  USE kinds, ONLY : DP
  REAL(DP), ALLOCATABLE ::  ef_sh(:) ! E_Fermi shift

  PRIVATE
  PUBLIC :: set_efsh, read_efsh, write_efsh !subroutines
  PUBLIC :: ef_sh                       !variables
  !
  LOGICAL,SAVE  :: first = .true. ! used for initializations
  REAL(DP),SAVE :: dos_ef         ! density of states at Ef

CONTAINS
!-----------------------------------------------------------------------
SUBROUTINE set_efsh (drhoscf, imode0, irr, npe)
  !-----------------------------------------------------------------------
  !  This routine calculates the FermiEnergy shift
  !   and stores it in the variable ef_sh
  !
  USE kinds,          ONLY : DP
  USE constants,      ONLY : eps8, rytoev
  USE io_global,      ONLY : stdout
  USE cell_base,      ONLY : omega
  USE pwcom,          ONLY : nbnd, degauss, ngauss, et, ef
  USE fft_interfaces, ONLY : fwfft
  USE fft_base,       ONLY : dfftp
  USE gvect,          ONLY : gg !, nl
  USE qpoint,         ONLY : nksq
  USE modes,          ONLY : npertx
  USE mp_pools,       ONLY : inter_pool_comm, intra_pool_comm
  USE mp,             ONLY : mp_sum
  USE kplus3q,        ONLY : kplusq, nbnd_max
  USE control_lr,     ONLY : nbnd_occ
  IMPLICIT NONE
  INTEGER :: npe, imode0, irr
  ! input: the number of perturbation
  ! input: the position of the current mode
  ! input: index of the current irr. rep.
  COMPLEX(DP) :: drhoscf (dfftp%nnr, npe)
  ! input: variation of the charge density

  INTEGER :: ipert, ik, ikk, ibnd
  ! counters
  COMPLEX(DP) :: delta_n, def (npertx)
  ! the change in electron number
  ! the change of the Fermi energy for each perturbation
  REAL(DP) :: weight, wdelta
  ! kpoint weight
  ! delta function weight
  REAL(DP),EXTERNAL :: w0gauss
  ! Used for initialization
  !
  CALL start_clock('set_efsh')
  ! first call: calculates density of states at Ef
  !
  IF (first) THEN
     first = .false.
     dos_ef = 0.d0
     DO ik = 1, nksq
        ikk = kplusq(0)%ikqs(ik)
        weight = kplusq(0)%wk(ik)
        DO ibnd = 1, nbnd_max !_occ(ikk)
           wdelta = w0gauss((ef - et(ibnd, ikk))/degauss, ngauss) / degauss
           dos_ef = dos_ef + weight * wdelta
!             write(*,'(4i6,99e24.12)') ik, ibnd, ikk, ngauss, &
!               (ef - et(ibnd, ikk))/degauss, w0gauss((ef - et(ibnd, ikk))/degauss, ngauss), weight, ef, et(ibnd,ikk), degauss
        ENDDO
     ENDDO
#ifdef __MPI
     call mp_sum( dos_ef, inter_pool_comm )
#endif
  ENDIF
  !
  IF(dos_ef<1.d-8) THEN
    WRITE(stdout,'(a,1e24.12,/,a,/,a)')&
        "WARNING! very low DOS at Fermi energy:", dos_ef, &
        "probably not enough k-points for this smearing! Or system is actually an",&
         "insulator. Wrong results are likely: setting E_f shift to zero."
      CALL stop_clock('set_efsh')
      def(1:npe) = 0._dp
      RETURN
  ENDIF
  !
  ! determines Fermi energy shift (such that each pertubation is neutral)
  !
  DO ipert = 1, npe
     call  fwfft('Rho', drhoscf(:,ipert), dfftp)
#ifdef __MPI
     delta_n = (0._dp, 0._dp)
     IF (gg(1) < eps8) delta_n = omega * drhoscf (dfftp%nl(1), ipert)
     call mp_sum ( delta_n, intra_pool_comm )
#else
     delta_n = omega * drhoscf(dfftp%l(1), ipert)
#endif
!      print*, ipert, delta_n, dos_ef, drhoscf (nl(1), ipert)
     def(ipert) = - delta_n/dos_ef
!      PRINT*,"fermi stuff:", delta_n, drhoscf(nl(1), ipert), dos_ef
  ENDDO
  !
  ! symmetrizes the Fermi energy shift
  CALL sym_def1 (def, irr)
  !
  ! impose ef_sh to be real
  DO ipert = 1, npe
     ef_sh(imode0 + ipert) = DBLE(def(ipert))
  ENDDO
!  print*, "WARNING: EFSH=0"
!  ef_sh = 0._dp
  !
!   WRITE( stdout, '(11x,"--> E_f shift (prt: Ry/eV) =",4(1i3,":",1f11.6,1f8.3," /",1f11.6,1f8.3))')  &
!        (ipert, def(ipert), rytoev*def(ipert), ipert = 1, npe)
  WRITE( stdout, '(11x,"--> E_f shift (prt: Ry/eV) =",(1i3,":",1f14.4,1f8.3," /",1f14.4,1f8.3))')  &
       1, def(1), rytoev*def(1)
  IF(npe>1) &
  WRITE( stdout, '(39x,1i3,":",1f14.4,1f8.3," /",1f14.4,1f8.3)')  &
       (ipert, def(ipert), rytoev*def(ipert), ipert=2,npe)
  !
  CALL stop_clock('set_efsh')
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE set_efsh
!----------------------------------------------------------------------
!
!---------------------------------------------------------------------
SUBROUTINE sym_def1 (def, irr)
  !---------------------------------------------------------------------
  ! Symmetrizes the first order changes of the Fermi energies of an
  ! irreducible representation. These objects are defined complex because
  ! perturbations may be complex
  !
  ! Used in the q=0 metallic case only.
  !
  USE kinds,        ONLY : DP
  USE d3_symmetry,  ONLY : sg => sym_gamma
  USE lr_symm_base, ONLY : nsymq, minus_q
  !USE io_global,   ONLY : stdout
  IMPLICIT NONE
  INTEGER :: irr ! input: the representation under consideration
  COMPLEX(DP) :: def (sg%npertx) ! inp/out: the fermi energy changes

  INTEGER :: ipert, jpert, isym, irot
  ! counter on perturbations
  ! counter on perturbations
  ! counter on symmetries
  ! the rotation

  COMPLEX(DP),ALLOCATABLE :: w_def (:) ! the fermi energy changes (work array)
  !

  DO ipert = 1, sg%npert(irr)
     def (ipert) = DBLE(def(ipert))
  ENDDO
  !
  ALLOCATE(w_def(sg%npertx))
  IF (minus_q) THEN
     w_def = (0.d0, 0.d0)
     DO ipert = 1, sg%npert(irr)
        DO jpert = 1, sg%npert(irr)
           w_def (ipert) = w_def(ipert) + sg%tmq(jpert, ipert, irr) &
                * def (jpert)
        ENDDO
     ENDDO
     DO ipert = 1, sg%npert(irr)
        def (ipert) = 0.5d0 * (def(ipert) + CONJG(w_def(ipert)) )
     ENDDO
  ENDIF
  !
  IF (nsymq == 1) THEN
    DEALLOCATE(w_def)
    RETURN
  ENDIF
  !
  ! Here we symmetrize with respect to the small group of q
  !
  w_def(:) = (0.d0, 0.d0)
  DO ipert = 1, sg%npert(irr)
     DO isym = 1, nsymq !sg%nsymq <-- k points have the symmetries of the triplet 
        irot = sg%irgq(isym)
!        print*, "Symmetrizing ef with op.", irot
        DO jpert = 1, sg%npert(irr)
           w_def(ipert) = w_def(ipert) &
                         + sg%t(jpert, ipert, irot, irr) * def(jpert)
!            print*, "sym_def1", sg%nsymq, irr, ipert, isym, jpert, sg%t(jpert, ipert, irot, irr)
        ENDDO
     ENDDO
  ENDDO

!   WRITE(stdout,'(15x,a,i3,a)') "efsh symmetrized with", sg%nsymq, " sym. ops."
  !
  ! normalize and exit
  !
  def (:) = w_def(:) / DBLE(sg%nsymq)
  !
  DEALLOCATE(w_def)
  !
  RETURN
  !---------------------------------------------------------------------
END SUBROUTINE sym_def1
!---------------------------------------------------------------------
!
!----------------------------------------------------------------------
SUBROUTINE read_efsh()
  !-----------------------------------------------------------------------
  !
  ! Reads the shift of the Fermi Energy
  !
  USE pwcom,      ONLY : lgauss
  USE d3_iofiles, ONLY : iuef
  USE io_global,  ONLY : ionode, ionode_id
  USE mp,         ONLY : mp_bcast
  USE mp_world,   ONLY : world_comm
  USE kplus3q,    ONLY : kplusq
  !
  IMPLICIT NONE
  !
  INTEGER :: ios
  !
  IF (.not. lgauss) RETURN
  !
  ef_sh = 0._dp
  !
  IF (.not.kplusq(1)%lgamma .and. &
      .not.kplusq(2)%lgamma .and. &
      .not.kplusq(3)%lgamma ) THEN
    RETURN
  ENDIF
  !
  IF ( ionode ) THEN
     !
     REWIND (unit = iuef)
     READ (iuef, iostat = ios) ef_sh
     !
  END IF
  !
  CALL mp_bcast(ios, ionode_id, world_comm)
  CALL errore ('read_ef', 'reading iuef', ABS(ios) )
  !
  CALL mp_bcast( ef_sh, ionode_id, world_comm )
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE read_efsh
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
SUBROUTINE write_efsh()
  !-----------------------------------------------------------------------
  USE pwcom,     ONLY : lgauss
  USE d3_iofiles, ONLY : iuef
  USE io_global, ONLY : ionode
  USE mp,        ONLY : mp_bcast

  IF((.not.lgauss) .or. (.not.ionode)) RETURN
  !
  WRITE (iuef) ef_sh
  FLUSH(iuef)
  !
  RETURN
  !-----------------------------------------------------------------------
END SUBROUTINE write_efsh
!----------------------------------------------------------------------
!
!---------------------------------------------------------------------
END MODULE d3_efermi_shift
!---------------------------------------------------------------------
