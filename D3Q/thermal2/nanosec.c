
#include <time.h>

double c_nanosec()
{
        double static new_epoch = -1.;
        struct timespec T;
        clock_gettime(CLOCK_REALTIME, &T);
        if ( new_epoch < 0.) new_epoch = (double) T.tv_sec;
        return (((double)T.tv_sec-new_epoch) + ((double)T.tv_nsec)*1.e-9);
}
