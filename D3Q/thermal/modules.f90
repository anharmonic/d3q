!
! CONTAINS FOUR MODULES:
! MPI_BASE    : used instead of including mpif.h and to store mpi-related variables
! cdiagh_work : workig space to diagonalize the dynamical matrix (via LAPACK)
! constants   : physical constants
! common_variables : code variables
!
!----------------------------------------------------------------------------
Module MPI_BASE
  Include 'mpif.h'
  Integer :: nrank, nsize
!
 Contains
  !
  Subroutine start_parallel ()
     Call mpi_init (ierr)
     Call mpi_comm_rank (mpi_comm_world, nrank, ierr)
     Call mpi_comm_size (mpi_comm_world, nsize, ierr)
     nrank = nrank + 1
     Write (*,*) "MPI start: ", ierr, "     ", nrank, "/", nsize
  End Subroutine
  !
  Subroutine stop_parallel ()
     Call mpi_finalize (ierr)
  End Subroutine
  !
End Module
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
Module constants
  !----------------------------------------------------------------------------
  !
  !
  ! .. The constants needed everywhere
  !
  Implicit None
  !
  Save
  !
  ! ... Mathematical constants
  Integer, Parameter :: DP = Selected_Real_Kind (14, 200)
  !
  Real (DP), Parameter :: pi = 3.14159265358979323846_DP
  Real (DP), Parameter :: tpi = 2.0_DP * pi
  Real (DP), Parameter :: fpi = 4.0_DP * pi
  Real (DP), Parameter :: sqrtpi = 1.77245385090551602729_DP
  Real (DP), Parameter :: sqrttpi = 2.506628275
  Real (DP), Parameter :: sqrtpm1 = 1.0_DP / sqrtpi
  Real (DP), Parameter :: sqrt2 = 1.41421356237309504880_DP
  !
  ! ... Physical constants, SI (NIST CODATA 2006), Web Version 5.1
  !     http://physics.nist.gov/constants
  Real (DP), Parameter :: H_PLANCK_SI = 6.62606896E-34_DP ! J s
  Real (DP), Parameter :: H_BAR_ERG = 1.05457172647E-27_DP ! erg s
  Real (DP), Parameter :: H_BAR_SI = 1.05457172647E-34_DP ! J s
  Real (DP), Parameter :: H_BAR_RY = 4.837764376E-17_DP ! Ry s
  Real (DP), Parameter :: K_BOLTZMANN_SI = 1.3806504E-23_DP ! J K^-1
  Real (DP), Parameter :: ELECTRON_SI = 1.602176487E-19_DP ! C
  Real (DP), Parameter :: ELECTRONVOLT_SI = 1.602176487E-19_DP ! J
  Real (DP), Parameter :: ELECTRONMASS_SI = 9.10938215E-31_DP ! Kg
  Real (DP), Parameter :: HARTREE_SI = 4.35974394E-18_DP ! J
  Real (DP), Parameter :: RYDBERG_SI = HARTREE_SI / 2.0_DP ! J
  Real (DP), Parameter :: BOHR_RADIUS_SI = 0.52917720859E-10_DP ! m
  Real (DP), Parameter :: AMU_SI = 1.660538782E-27_DP ! Kg
  Real (DP), Parameter :: C_SI = 2.99792458E+8_DP ! m sec^-1
  !
  ! ... Physical constants, atomic units:
  ! ... AU for "Hartree" atomic units (e = m = hbar = 1)
  ! ... RY for "Rydberg" atomic units (e^2=2, m=1/2, hbar=1)
 !
  Real (DP), Parameter :: K_BOLTZMANN_AU = K_BOLTZMANN_SI / HARTREE_SI
  Real (DP), Parameter :: K_BOLTZMANN_RY = K_BOLTZMANN_SI / RYDBERG_SI
  !
  ! ... Unit conversion factors: energy and masses
  !
  Real (DP), Parameter :: AUTOEV = HARTREE_SI / ELECTRONVOLT_SI
  Real (DP), Parameter :: RYTOEV = AUTOEV / 2.0_DP
  Real (DP), Parameter :: AMU_AU = AMU_SI / ELECTRONMASS_SI
  Real (DP), Parameter :: AMU_RY = AMU_AU / 2.0_DP
  !
  ! ... Unit conversion factors: atomic unit of time, in s and ps
  !
  Real (DP), Parameter :: AU_SEC = H_PLANCK_SI / tpi / HARTREE_SI
  Real (DP), Parameter :: AU_PS = AU_SEC * 1.0E+12_DP
  !
  ! ... Unit conversion factors: pressure (1 Pa = 1 J/m^3, 1GPa = 10 Kbar )
  !
  Real (DP), Parameter :: AU_GPA = HARTREE_SI / BOHR_RADIUS_SI ** 3 / 1.0E+9_DP
  Real (DP), Parameter :: RY_KBAR = 10.0_DP * AU_GPA / 2.0_DP
  !
  ! ... Unit conversion factors: 1 debye = 10^-18 esu*cm
  ! ...                                  = 3.3356409519*10^-30 C*m
  ! ...                                  = 0.208194346 e*A
  ! ... ( 1 esu = (0.1/c) Am, c=299792458 m/s)
  !
  Real (DP), Parameter :: DEBYE_SI = 3.3356409519_DP * 1.0E-30_DP ! C*m
  Real (DP), Parameter :: AU_DEBYE = ELECTRON_SI * BOHR_RADIUS_SI / DEBYE_SI
  !
  Real (DP), Parameter :: eV_to_kelvin = ELECTRONVOLT_SI / K_BOLTZMANN_SI
  Real (DP), Parameter :: ry_to_kelvin = RYDBERG_SI / K_BOLTZMANN_SI
  Real (DP), Parameter :: ry_to_cm1 = 13.6058d0 * 8065.5d0
  !
  ! .. Unit conversion factors: Energy to wavelength
  !
  Real (DP), Parameter :: EVTONM = 1E+9_DP * H_PLANCK_SI * C_SI / ELECTRONVOLT_SI
  Real (DP), Parameter :: RYTONM = 1E+9_DP * H_PLANCK_SI * C_SI / RYDBERG_SI
  !
  !  Speed of light in atomic units
  !
  Real (DP), Parameter :: C_AU = C_SI / BOHR_RADIUS_SI * AU_SEC
  !
  ! ... zero up to a given accuracy
  !
  Real (DP), Parameter :: eps4 = 1.0E-4_DP
  Real (DP), Parameter :: eps6 = 1.0E-6_DP
  Real (DP), Parameter :: eps8 = 1.0E-8_DP
  Real (DP), Parameter :: eps12 = 1.0E-12_DP
  Real (DP), Parameter :: eps14 = 1.0E-14_DP
  Real (DP), Parameter :: eps16 = 1.0E-16_DP
  Real (DP), Parameter :: eps24 = 1.0E-24_DP
  Real (DP), Parameter :: eps32 = 1.0E-32_DP
  !
  Real (DP), Parameter :: gsmall = 1.0E-12_DP
  !
  Real (DP), Parameter :: e2 = 2.0_DP ! the square of the electron charge
  Real (DP), Parameter :: degspin = 2.0_DP ! the number of spins per level
  !
  Real (DP), Parameter :: amconv = AMU_RY
  Real (DP), Parameter :: uakbar = RY_KBAR
  Real (DP), Parameter :: bohr_radius_cm = BOHR_RADIUS_SI * 100.0_DP
  Real (DP), Parameter :: BOHR_RADIUS_ANGS = bohr_radius_cm * 1.0E8_DP
  Real (DP), Parameter :: BOHR_TO_M = 0.52917721092E-10_DP
  Real (DP), Parameter :: ANGSTROM_AU = 1.0_DP / BOHR_RADIUS_ANGS
  Real (DP), Parameter :: DIP_DEBYE = AU_DEBYE
  Real (DP), Parameter :: AU_TERAHERTZ = AU_PS
  Real (DP), Parameter :: AU_TO_OHMCMM1 = 46000.0_DP ! (ohm cm)^-1
  Real (DP), Parameter :: RY_TO_THZ = 1.0_DP / AU_TERAHERTZ / fpi
  Real (DP), Parameter :: RY_TO_CMM1 = 1.E+10_DP * RY_TO_THZ / C_SI
  Real (DP), Parameter :: RY_TO_WATT = 45.059494E-3_DP
!
  Complex (DP), Parameter :: cone = cmplx (1.d0, 0.d0), czero = cmplx (0.d0, 0.d0)
  !
End Module constants
!
!
Module cdiagh_work
  Use constants, Only: DP
  Implicit None
  Save
  Integer :: lwork, n_first
  Real (Kind=DP), Allocatable :: rwork (:)
  Complex (Kind=DP), Allocatable :: work (:)
End Module cdiagh_work
!
!
Module common_variables
  Use constants, Only: DP
  Implicit None
  !------------------------------------------------------------------------------------
  ! Variables read in input
  Character (Len=80) :: filemat2 ! name of file containg 2nd order force constants
  Character (Len=80) :: filemat3 ! name of file containg 3rd order force constants
  Integer :: nq1 (3)! dimensions of the external q-point grid (crystall coords.)
  Real (DP) :: deltaq1 (3)!      shift of the external q-point grid (crystall coords.)
  Integer :: nq2 (3)! dimansions of the internal q-point grid (crystall coords.)
  Real (DP) :: deltaq2 (3)!      shift of the internal q-point grid (crystall coords.)
  Logical :: asr2, asr3
  Integer :: ispeed2, ispeed3, asr3iters
  Real (DP) :: alpha ! esponent of Q_0 for the preconditioning case CONJUGATE GRADIENT
  Logical :: ltest
  Real (DP) :: x_iso ! percentage of isotopes
!
  Integer, Parameter :: ntempx = 20
  Integer :: ntemp, nxcryst
  Integer :: xcryst (3)
  Real (DP) :: temp (ntempx), sigma (ntempx), sigmam1 (ntempx), tempm1 (ntempx)
  Logical, Allocatable :: ltobedone (:)
  !------------------------------------------------------------------------------------
  ! variables for parallel part
  Integer :: ierr
  Integer :: nq1loc
  Integer, Parameter :: nq1loc_max = 300000
  Integer :: iq1_loc (nq1loc_max)
  !------------------------------------------------------------------------------------
  ! output file number
  Integer :: iun1
  !------------------------------------------------------------------------------------
  ! Parameters describing the system
  Integer :: ntyp, nat, ibrav
  Integer, Parameter :: ntypx = 10
  Real (DP) :: celldm (6), at (3, 3), bg (3, 3)
  Real (DP), Allocatable :: tau (:,:)
  Character (Len=3) :: atm (ntypx)
  Real (DP) :: amass (ntypx)
  Real (DP) :: omega
  Integer, Allocatable :: ityp (:)
  Real (DP) :: fracmass2 ! fracmass = delta_m / amass_b  ...for isotopic distribution
  Real (DP) :: const_iso
!
  Integer :: nq1tot, nq2tot
  Integer :: nat3, nat32, nat33
  Real (DP), Allocatable :: sqrtm1 (:)
  Integer :: dim1, dimrun
  Integer :: dimG
  Real (DP) :: Lcas ! Casimir Length
  Real (DP) :: Lcasm1 ! Casimir Length-1
  Real (DP) :: Fcas !correction factor for border scattering
  Logical :: lbordo
 ! -----------------------------------------------------------------------------------
 ! Indeces
  Integer :: iq1_, ipq2_, imq2_, ipq3_, imq3_, ipq3b_
!
  !------------------------------------------------------------------------------------
  ! 2nd order force constants
  Integer :: nRbig2_in ! # R-points, as read from file
  Integer, Allocatable :: iRbig2_in (:,:)!   R-points, as read from file
  Complex (DP), Allocatable :: mat2_in (:,:,:)!   Force constants
!
  Integer :: nRbig2 (3)! number or R-points,      after reordering
  Integer :: nRbig2t ! total number or R-poins, after reordering
  Integer, Allocatable :: iRbig2 (:,:)! R-points,                after reordering (crystal coords)
  Real (DP), Allocatable :: Rbig2 (:,:)! same as iRbig2
  Complex (DP), Allocatable :: mat2 (:,:,:)! Force constants,         after reordering
!
!
  !------------------------------------------------------------------------------------
  ! 3rd order force constants (they depends on two R vectors, associated to the 2nd and 3rd indexes)
  ! The R-vector associated to the 1st index is zero
  Integer :: nRbig3_in ! # R-points, as read from file
  Integer, Allocatable :: iRbig3_in (:,:,:)!   R-points, as read from file
  Complex (DP), Allocatable :: mat3_in (:,:,:,:)!   Force constants
  Integer :: nRbig31 (3)! # R-points of 2nd index, after rerdering
  Integer :: nRbig31t
  Integer, Allocatable :: iRbig31 (:,:)!   R-points of 2nd index, after rerdering (crystal coords)
  Real (DP), Allocatable :: Rbig31 (:,:)! same as iRbig31
  Integer :: nRbig32 ! # R-points of 3rd index, after rerdering
  Integer, Allocatable :: iRbig32 (:,:)!   R-points of 3rd index, after rerdering (crystal coords)
  Real (DP), Allocatable :: Rbig32 (:,:)! same as iRbig32
  Complex (DP), Allocatable :: mat3_rr (:,:,:)! Force constants,         after reordering
!
  !------------------------------------------------------------------------------------
  !
  Integer, Parameter :: ic_out = 3, ic_med = 2, ic_inn = 1
!  integer, parameter, dimension=3 :: ic_corr(:) = 1,2,3
!
  !------------------------------------------------------------------------------------
  ! Auxiliary quantities which are used during the Fourier interpolation
  Complex (DP), Allocatable :: mat3_r (:,:)
  Complex (DP), Allocatable :: mat3_rzy (:,:), mat3_rz (:,:)! ispeed=3
  Complex (DP), Allocatable :: mat3_rzyB (:,:), mat3_rzB (:,:)! ispeed=3
  Complex (DP), Allocatable :: mat3q (:,:,:), mat3qb (:,:,:)
!
  Complex (DP), Allocatable :: phasem3_q1t (:)! ispeed=2/3
  Complex (DP), Allocatable :: phasem3_q2t (:), phasem3_q2tc (:), phasem3_q2t_ (:,:), phasem3_q2t_init (:)! ispeed=2
  Complex (DP), Allocatable :: phasem3_q2 (:,:), phasem3_q2_ (:,:), phasem3_q2_init (:,:)! ispeed=3
  Complex (DP), Allocatable :: phasem3_q2c (:,:)! ispeed=3
!
!  complex(DP), allocatable :: d3mm_aux(:,:,:), d3mm_aux1(:,:,:), d3mm_aux2(:,:,:)
  Complex (DP), Allocatable :: d3mm_aux (:,:,:,:), d3mm_aux1 (:,:,:,:), d3mm_aux2 (:,:,:,:)
!
  !------------------------------------------------------------------------------------
  ! After the Fourier interpolation:
  Real (DP) :: q1d (3), q2d (3), mq2d (3), q3d (3), q3db (3), mq3d (3) ! The three q-points
  Complex (DP), Allocatable :: zz1_x (:,:,:), zz2_x (:,:,:), zz3_x (:,:,:), &
                               zz3b_x (:,:,:), zz2b_c_x (:,:,:)! The corresponding eigenvectors
  Real (DP), Allocatable :: freq1 (:), freq2 (:), freq3 (:), freq3b (:)! The corresponding frequencies
  Real (DP), Allocatable :: bos1 (:), bos2 (:), bos3 (:), bos3b (:)!  Bose-Einst populations
  Real (DP), Allocatable :: d3mmx (:,:,:,:), d3mmx_1 (:,:,:,:)!  The scattering coefficients
  ! 
  Real (DP), Allocatable :: lambd (:,:,:,:), lambdiso (:,:,:,:) ! thermal transport
!
  !--------------------------------------------------------------------------------
  !-  dyndiag
  Logical :: on ! true: evaluate new freq and zz; false: take previously evaluated values
  !----------------------------------------------------------------------------------------
  !umklapp variables
  Real (DP), Allocatable :: q1ws (:,:), q2ws (:,:), q3ws (:,:), q3bws (:,:)
  Integer, Allocatable :: neqq1 (:), neqq2 (:), neqq3 (:), neqq3b (:)
  Real (DP) :: const_kuk, const_k
  !
  !----------------------------------------------------------------------------------------
  ! group velocities
  Real (DP), Allocatable :: velph (:,:,:)
  Real (DP), Allocatable :: freq1_ (:,:), freq3_ (:,:), freq3b_ (:,:), freq2_ (:,:)
  Complex (DP), Allocatable :: zz1n (:,:,:,:), zz3n (:,:,:,:), zz3bn (:,:,:,:), zz2n (:,:,:,:)
! NEW: Eigenvactors evaluated in agreement with velocity diagonalization
  !----------------------------------------------------------------------------------------
  ! thermal conductivity SMA
!
  Real (DP) :: const_cond (ntempx), tcondc_ (3, 3, ntempx), tcondc_SMA (3, 3)
  Real (DP), Allocatable :: velph_ (:,:,:)
!
  !----------------------------------------------------------------------------------------
  ! thermal conductivity FIS & CG
  Real (DP), Allocatable :: Q_0 (:,:,:,:), Q_0rad (:,:,:,:), Q_0radm1 (:,:,:,:)
  Real (DP), Allocatable :: sum_IT (:,:,:,:), sum_ITq (:,:,:,:), sum_ITq2 (:,:,:,:), sum_IT2 (:,:,:,:,:,:)
  Real (DP), Allocatable :: F_0 (:,:,:,:), F_old (:,:,:,:), F_new (:,:,:,:), F_rs (:,:,:,:), matrix (:,:,:,:)
  Real (DP) :: const_condFI, tcondFI_ (3, 3, ntempx), tcond (3, 3, ntempx),&
               tcond2 (3, 3, ntempx), tcond3 (3, 3, ntempx), tcondrs (3, 3, ntempx)
!
  Real (DP), Allocatable :: ggrd (:,:), gp (:,:), ggt (:,:)
  Real (DP), Allocatable :: hh (:,:), hh_old (:,:)
  Real (DP) :: treshk
!
  !------------------------------------------------------------------------------------
  ! Variables used for the timing
  Real (DP) :: T0, T1, ttot_0, ttot_1, tinit_0, tinit_1
  Real (DP) :: tphase, tsetup2, tsetup3, tdindiag, tbose, td3mm, tmodes
  Real (DP), External :: nanosec
  !------------------------------------------------------------------------------------
  !
  CONTAINS
  !-----------------------------------------------------------------------------------------
  Subroutine alloc_dealloc_common (isw)
  !-----------------------------------------------------------------------------------------
     Use cdiagh_work
     Implicit None
     Integer :: isw
!
     Integer :: nRbig2x
!
     If (isw == 1) Then
        Allocate (mat3_r(nat33, nRbig31t))
        Allocate (mat3q(nat3, nat3, nat3))
        Allocate (mat3qb(nat3, nat3, nat3))

        Allocate (zz1_x(nat3, nat3, 3), zz2_x(nat3, nat3, 3), &
                  zz3_x(nat3, nat3, 3), zz3b_x(nat3, nat3, 3),& 
                  zz2b_c_x(nat3, nat3, 3) )
        Allocate (zz1n(nat3, nat3, nq1tot, 3))
        Allocate (zz2n(nat3, nat3, nq2tot, 3))
        Allocate (zz3n(nat3, nat3, nq2tot, 3))
        Allocate (zz3bn(nat3, nat3, nq2tot, 3))
        Allocate (freq1(nat3), freq2(nat3), freq3(nat3), freq3b(nat3))
!
        dim1 = nat3 * nq1tot
        dimrun = 3 * ntemp
        dimG = nat3 * nq1tot * 3 * ntemp
        Allocate (ggrd(dim1, dimrun), gp(dim1, dimrun), ggt(dim1, dimrun))
        Allocate (hh(dim1, dimrun), hh_old(dim1, dimrun))
!
        If (ispeed3 == 3) Then
           Allocate (phasem3_q1t(nRbig32))
           Allocate (mat3_rzy(nat33, nRbig31(ic_inn)*nRbig31(ic_med)), mat3_rz(nat33, nRbig31(ic_inn)))
           Allocate (mat3_rzyB(nat33, nRbig31(ic_inn)*nRbig31(ic_med)), mat3_rzB(nat33, nRbig31(ic_inn)))
           nRbig2x = Max (Max(nRbig31(1), nRbig31(2)), nRbig31(3))
           Allocate (phasem3_q2(nRbig2x, 3), phasem3_q2_(nRbig2x, 3), phasem3_q2_init(nRbig2x, 3))
           Allocate (phasem3_q2c(nRbig2x, 3))
        End If
        Allocate (d3mm_aux(nat3, nat3, nat3, 3), d3mm_aux1(nat3, nat3, nat3, 3), d3mm_aux2(nat3, nat3, nat3, 3))

        Allocate (d3mmx(nat3, nat3, nat3, 3), d3mmx_1(nat3, nat3, nat3, 3))
        Allocate (bos1(nat3), bos2(nat3), bos3(nat3), bos3b(nat3))
        Allocate (lambd(nat3, ntemp, 3, nq1tot))
        Allocate (lambdiso(nat3, ntemp, 3, nq1tot))
        Allocate (velph(3, nat3, nq1tot))
        Allocate (freq1_(nat3, nq1tot))
        Allocate (freq2_(nat3, nq2tot))
        Allocate (freq3_(nat3, nq2tot))
        Allocate (freq3b_(nat3, nq2tot))
        Allocate (velph_(3, 3, nq1tot))
        Allocate (Q_0(nat3, nq1tot, 3, ntemp))

        Allocate (Q_0rad(nat3, nq1tot, 3, ntemp))
        Allocate (Q_0radm1(nat3, nq1tot, 3, ntemp))

        Allocate (F_0(nat3, nq1tot, 3, ntemp))
        Allocate (F_old(nat3, nq1tot, 3, ntemp))
        Allocate (F_new(nat3, nq1tot, 3, ntemp))
        Allocate (sum_IT(nat3, nq1tot, 3, ntemp))

        Allocate (sum_ITq(nat3, nq1tot, 3, ntemp))
        Allocate (sum_ITq2(nat3, nq1tot, 3, ntemp))
        Allocate (matrix(nat3, nq1tot, 3, ntemp))
!
        Allocate (neqq1(nq1tot), neqq2(nq2tot), neqq3(nq2tot), neqq3b(nq2tot))
!
        Allocate (q1ws(3, nq1tot), q2ws(3, nq2tot), q3ws(3, nq2tot), q3bws(3, nq2tot))
!
!______________________________________________________________________________
     Else If (isw ==-1) Then  ! DEALLOCATE
!______________________________________________________________________________
    
 ! These were allocated in read_mat2R
        Deallocate (ityp)
        Deallocate (tau)
        Deallocate (sqrtm1)
 ! These were allcoated n change_m2r_index
        Deallocate (iRbig2, Rbig2, mat2)
 ! These were allcoated n change_m3r_index
        Deallocate (iRbig31, Rbig31, iRbig32, Rbig32, mat3_rr)
 ! These were allocated by cdiagh
        Deallocate (work)
        Deallocate (rwork)
!
        Deallocate (mat3_r)
        Deallocate (mat3q)
        Deallocate (mat3qb)
        Deallocate (zz1_x, zz2_x, zz3_x, zz3b_x, zz2b_c_x)
        Deallocate (zz1n, zz2n, zz3n, zz3bn)
        Deallocate (freq1, freq2, freq3, freq3b)
        Deallocate (freq1_)
        Deallocate (freq2_)
        Deallocate (freq3_)
        Deallocate (freq3b_)
        Deallocate (velph)
        Deallocate (velph_)
        Deallocate (Q_0)
        Deallocate (Q_0rad)
        Deallocate (Q_0radm1)
        Deallocate (F_0)
        Deallocate (F_new)
        Deallocate (F_old)
        Deallocate (sum_IT)
        Deallocate (sum_ITq)
        Deallocate (sum_ITq2)
        Deallocate (matrix)
        Deallocate (ltobedone)
!
        If (ALLOCATED(neqq1)) DEALLOCATE (neqq1)
        If (ALLOCATED(neqq2)) DEALLOCATE (neqq2)
        If (ALLOCATED(neqq3)) DEALLOCATE (neqq3)
        If (ALLOCATED(neqq3b)) DEALLOCATE (neqq3b)
        If (ALLOCATED(q1ws)) DEALLOCATE (q1ws)
        If (ALLOCATED(q2ws)) DEALLOCATE (q2ws)
        If (ALLOCATED(q3ws)) DEALLOCATE (q3ws)
        If (ALLOCATED(q3bws)) DEALLOCATE (q3bws)
!
        If (ispeed3 == 3) Then
           Deallocate (phasem3_q1t)
           Deallocate (phasem3_q2t, phasem3_q2tc, phasem3_q2t_, phasem3_q2t_init)
           Deallocate (phasem3_q2, phasem3_q2_, phasem3_q2_init)
           Deallocate (mat3_rzy, mat3_rz)
           Deallocate (mat3_rzyB, mat3_rzB)
           Deallocate (phasem3_q2, phasem3_q2c, phasem3_q2_, phasem3_q2_init)
        End If
        Deallocate (d3mm_aux)
        Deallocate (d3mm_aux1)
        Deallocate (d3mm_aux2)
        Deallocate (d3mmx)
        Deallocate (d3mmx_1)
        Deallocate (bos1, bos2, bos3, bos3b)
        Deallocate (lambd)
        Deallocate (lambdiso)
!
        If (ALLOCATED(ggrd)) DEALLOCATE (ggrd)
        If (ALLOCATED(gp)) DEALLOCATE (gp)
        If (ALLOCATED(hh)) DEALLOCATE (hh)
        If (ALLOCATED(hh_old)) DEALLOCATE (hh_old)
        If (ALLOCATED(ggt)) DEALLOCATE (ggt)
!
!
     End If
!
!
!
!
     Return
  End Subroutine alloc_dealloc_common
!----------------------------------------------------------------------------------
!
End Module common_variables
