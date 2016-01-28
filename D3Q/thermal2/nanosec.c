/*
! Written by Lorenzo Paulatto (2013-2015) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
*/

#include <time.h>

double c_nanosec()
{
        double static new_epoch = -1.;
        struct timespec T;
        clock_gettime(CLOCK_REALTIME, &T);
        if ( new_epoch < 0.) new_epoch = (double) T.tv_sec;
        return (((double)T.tv_sec-new_epoch) + ((double)T.tv_nsec)*1.e-9);
}
