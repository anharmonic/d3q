/*
! Written by Lorenzo Paulatto (2013-2015) IMPMC @ UPMC / CNRS UMR7590
!  released under the CeCILL licence v 2.1
!  <http://www.cecill.info/licences/Licence_CeCILL_V2.1-fr.txt>
*/
/*#if defined(__APPLE__)
double c_nanosec()
{ return 0.;  }
#else*/

#include <time.h>
#include <sys/time.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

double c_nanosec()
{
        double static new_epoch = -1.;
        struct timespec ts;

#if defined(__MACH__)
        clock_serv_t cclock;
        mach_timespec_t mts;
        host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
        clock_get_time(cclock, &mts);
        mach_port_deallocate(mach_task_self(), cclock);
        ts.tv_sec = mts.tv_sec;
        ts.tv_nsec = mts.tv_nsec;
#else
        clock_gettime(CLOCK_REALTIME, &ts);
#endif
        if ( new_epoch < 0.) new_epoch = (double) ts.tv_sec;
        return (((double)ts.tv_sec-new_epoch) + ((double)ts.tv_nsec)*1.e-9);
}

