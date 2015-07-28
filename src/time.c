//
//  time.c
//  
//
//  Created by Johannes Dietschreit on 29.05.15.
//
//

#include <sys/time.h> /* for time of day */
#include <stdlib.h> 

#include "time.h"


/* get the time of day from the system clock, and store it (in seconds) */
double time_of_day(void) {
#if defined(_MSC_VER)
    double t;
    
    t = GetTickCount();
    t = t / 1000.0;
    
    return t;
#else
    struct timeval tm;
    //    struct timezone tz;
    
    //    gettimeofday(&tm, &tz);
    gettimeofday(&tm, NULL);
    return((double)(tm.tv_sec) + (double)(tm.tv_usec)/1000000.0);
#endif
}
