!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! Crappy libraries, to be deleted
module d3_tools
  !
  USE kinds,     only : DP
  USE constants, ONLY : tpi
  !
  INTEGER,PARAMETER  :: NTYPX     = 10
  REAL(DP),PARAMETER :: eps8      = 1.d-6
  INTEGER, PARAMETER :: seek      = 6
  !
  TYPE d3_system_info
    ! atoms
    INTEGER              :: ntyp
    REAL(DP), ALLOCATABLE:: amass(:)
    CHARACTER(len=3  ), ALLOCATABLE   :: atm(:)
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
    REAL(DP)             :: q1(3,48), q2(3,48), q3(3,48)
    ! phonon switches (mostly unused here)
!     REAL(DP)             :: epsil(3,3)
!     LOGICAL              :: lrigid
    ! dynamical matrix
    COMPLEX(DP), ALLOCATABLE :: phiq(:,:,:,:,:,:)
    COMPLEX(DP), ALLOCATABLE :: Do(:,:,:), Dn(:,:,:)
  END TYPE
  !
contains
!-----------------------------------------------------------------------
subroutine read_d3(iudyn, nat, ntyp, ityp, atm, amass, tau, ibrav, at, celldm, xp,xq,xpq, phi)
  !-----------------------------------------------------------------------
  implicit none
  !
  ! input variables
  !
  integer,intent(in)  :: iudyn
  integer,intent(out) :: nat, ntyp,  ibrav
  integer,allocatable,intent(out) :: ityp(:)
  real(dp),intent(out):: celldm(6)
  real (DP),intent(out) :: xp(3),xq(3),xpq(3), at(3,3)
  real (DP),allocatable,intent(out) :: tau (:,:), amass(:)
  character(len=3),allocatable :: atm(:)
  ! unit number
  ! number of atom in the unit cell
  complex (DP),allocatable,intent(out) :: phi (:,:,:,:,:,:)
  !  derivative of the dynamical matrix
  character(len=1024) :: dummy
  integer :: dummy1, nt, na, j
  !
  ! All variables read from file that need dynamical allocation
  !
  !
  read (iudyn, *) dummy
   read (iudyn, *)
  read (iudyn, *) ntyp, nat, ibrav, celldm
  if(ibrav==0)then
    read(iudyn,*) at(:,1)
    read(iudyn,*) at(:,2)
    read(iudyn,*) at(:,3)
  endif

  allocate(tau(3,nat))
  allocate(ityp(nat))
  allocate(atm(ntyp))
  allocate(amass(ntyp))
  allocate(phi(3,3,3,nat,nat,nat))

  do nt = 1, ntyp
     read (iudyn, * ) dummy1, atm (nt), amass (nt)
  enddo
  do na = 1, nat
     read (iudyn, *) dummy1, ityp (na) , (tau (j, na) , j = 1, 3)
  enddo

  call read_d3dyn(xp,xq,xpq, phi, nat, iudyn)

end subroutine read_d3
!
!-----------------------------------------------------------------------
subroutine read_d3dyn (xp,xq,xpq, phi, nat, iudyn, wrmode_)
  !-----------------------------------------------------------------------
  !
  USE kinds, only : DP
  implicit none
  !
  ! input variables
  !
  integer,intent(in) :: iudyn, nat
  ! unit number
  ! number of atom in the unit cell
  complex (DP),intent(out) :: phi (3, 3, 3, nat, nat, nat)
  !  derivative of the dynamical matrix
  real (DP),intent(out) :: xp(3), xq(3), xpq(3)
  ! the q vector
  logical,optional,intent(in) :: wrmode_ (3 * nat)
  logical :: wrmode (3 * nat)
  ! if .true. this mode is to be read
  !
  ! local variables
  !
  character(len=1024) :: dummyc
  integer :: na, nb, nc, icar, jcar, kcar, i, dummy1, dummy2
  ! counters on atoms
  ! cartesian coordinate counters
  ! generic counter
  if(present(wrmode_)) then
    wrmode = wrmode_
  else
    wrmode = .true.
  endif

  read(iudyn, *)
  read(iudyn, *) dummyc
  read(iudyn, *)
!   read(iudyn, *)
  read(iudyn, '(a12,3f14.9)') dummyc, xp(1:3)
  read(iudyn, '(a12,3f14.9)') dummyc, xq(1:3)
  read(iudyn, '(a12,3f14.9)') dummyc, xpq(1:3)
!   read (iudyn, ) (xq (icar), icar = 1, 3)
!   read(iudyn, *)
  do i = 1, 3 * nat
     if (wrmode (i) ) then
        read(iudyn, *)
        read (iudyn, *) dummyc,dummy1
        read(iudyn, *)
          if(dummy1 .ne. i) &
            call errore('read_d3dyn', 'something wrong', 2)
        nc = (i - 1) / 3 + 1
        kcar = i - 3 * (nc - 1)
        do na = 1, nat
           do nb = 1, nat
!               read (iudyn, '(a1024)') dummyc
!               print*, dummyc
              read (iudyn, *) dummy1, dummy2
              if(dummy1 .ne. na .or. dummy2 .ne. nb) &
                call errore('read_d3dyn', 'something wrong', 1)
              do icar = 1, 3
!                  read (iudyn, '(3e24.12)') (phi (kcar, icar, jcar, nc, na, nb) &
                 read (iudyn, '(6e24.12)') (phi (kcar, icar, jcar, nc, na, nb), &
                                            jcar = 1, 3)
              enddo
           enddo
        enddo
     endif

  enddo

  return
! 9000 format(/,5x,'Third derivative in cartesian axes', &
!        &       //,5x,'q = ( ',3f14.9,' ) ',/)
end subroutine read_d3dyn

!----------------------------------------------------------------------------

FUNCTION check_int_linearcombination(test_vector, inverted_basis)
  IMPLICIT NONE
  REAL(DP),INTENT(IN) :: test_vector(3)
  REAL(DP),INTENT(IN) :: inverted_basis(3,3)
  REAL(DP) :: coefficients(3)
  LOGICAL :: check_int_linearcombination
  !
  coefficients = MATMUL(transpose(inverted_basis), test_vector)
  !
  ! wat20100312 : okay, I had to change this because it didn't work for negative translations.
  !               Feel free to change it back if you encounter problems.
  check_int_linearcombination = ALL(ABS(ABS(coefficients)-NINT(ABS(coefficients)))<eps8 ) 
  !check_int_linearcombination = ALL(ABS(coefficients-INT(coefficients+.5_dp))<eps8 )
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

! Function get_delta()
! 
! END Function
!
end module d3_tools
