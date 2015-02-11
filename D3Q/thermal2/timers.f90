!
! Written by Lorenzo Paulatto (2013-2015) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
!
! <<^V^\\=========================================//-//-//========//O\\//
MODULE timers
  USE nanoclock, ONLY : nanotimer
  IMPLICIT NONE
  
  TYPE(nanotimer) :: t_freq =   nanotimer("freq & bose"), &
                     t_sum  =   nanotimer("sum modes"), &
                     t_fc3int = nanotimer("fc3 interpolate"), &
                     t_fc3m2  = nanotimer("fc3 modulus sq"), &
                     t_fc3rot = nanotimer("fc3 rotate")
          
END MODULE