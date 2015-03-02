!
! Written by Lorenzo Paulatto (2013-2015) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
! <<^V^\\=========================================//-//-//========//O\\//
MODULE timers
  USE nanoclock, ONLY : nanotimer
  IMPLICIT NONE
  
  TYPE(nanotimer) :: t_freq =   nanotimer("ph interp & diag"), &
                     t_bose =   nanotimer("bose distrib"), &
                     t_sum  =   nanotimer("sum modes"), &
                     t_fc3int = nanotimer("fc3 interpolate"), &
                     t_fc3m2  = nanotimer("fc3 modulus sq"), &
                     t_fc3rot = nanotimer("fc3 rotate"), &
                     t_lwphph = nanotimer("lw ph-ph"), &
                     t_lwisot = nanotimer("lw isotopes"), &
                     t_lwcasi = nanotimer("lw casimir")
!    *      lw isotopes *       48.455432 s     *   272.221529 ms/call *    7.365 % wtime *          178 calls *
!    *       lw casimir *        0.092465 s     *     0.519465 ms/call *    0.014 % wtime *          178 calls *
!    *         lw ph-ph *      608.567777 s     *  3418.920093 ms/call *   92.498 % wtime *          178 calls *
!   * Contributions to ph-ph linewidth time:
!    *      freq & bose *      205.748674 s     *     0.070546 ms/call *   31.272 % wtime *      2916530 calls *
!    *        sum modes *       17.402060 s     *     0.007956 ms/call *    2.645 % wtime *      2187264 calls *
!    *  fc3 interpolate *      361.900451 s     *     0.496374 ms/call *   55.006 % wtime *       729088 calls *
!    *   fc3 modulus sq *        0.967292 s     *     0.001327 ms/call *    0.147 % wtime *       729088 calls *
!    *       fc3 rotate *       17.824916 s     *     0.024448 ms/call *    2.709 % wtime *       729088 calls *

  TYPE(nanotimer) :: t_tksma = nanotimer("SMA thermalk"), &
                     t_tksum = nanotimer("sum of TK terms")
END MODULE