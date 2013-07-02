/*
  Copyright (C) 2002-2006 Quantum ESPRESSO group
  This file is distributed under the terms of the
  GNU General Public License. See the file `License'
  in the root directory of the present distribution,
  or http://www.gnu.org/copyleft/gpl.txt .
*/

#include <time.h>
//#include <sys/resource.h>
//#include <unistd.h>



double nanosec_(double *new_epoch )
{
        struct timespec T;
        clock_gettime(CLOCK_REALTIME, &T);
        return (((double)T.tv_sec-new_epoch[0]) + ((double)T.tv_nsec)*1.e-9);
}

double NANOSEC_(double *new_epoch )
{
        struct timespec T;
        clock_gettime(CLOCK_REALTIME, &T);
        return (((double)T.tv_sec-new_epoch[0]) + ((double)T.tv_nsec)*1.e-9);
}

