#ifndef MACTIME_H
#define MACTIME_H

#ifdef __MACH__

#include <sys/time.h>

#define CLOCK_REALTIME 1234

int clock_gettime(int clk_id, struct timespec* t);

#endif /* __MACH__ */

#endif /* MACTIME_H */
