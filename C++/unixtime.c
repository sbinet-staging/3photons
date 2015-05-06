#include "unixtime.h"
#include <stdio.h>
#include <sys/types.h>
#include <sys/times.h>
#include <unistd.h>

     /* code for including or porting to Ada
        returns time from program start to call time */

     /* this works with SunOS and HP/UX. Not all Unixes have the same
        CPU time calls; your mileage might vary.     */

int unixtime ()
{
  struct tms tmsGot;

  times (&tmsGot);
  return (tmsGot.tms_utime);
}

int sysconf_sc_clk_tck ()
{
  //printf ("SYSCONF %d", sysconf (_SC_CLK_TCK));
  return (sysconf (_SC_CLK_TCK));
}

int usertime ()
{
  struct tms tmsGot;

  times (&tmsGot);
  return (tmsGot.tms_utime);
}

int systtime ()
{
  struct tms tmsGot;

  times (&tmsGot);
  return (tmsGot.tms_stime);
}
//            struct tms {
//                clock_t tms_utime;  /* dur�e utilisateur          */
//                clock_t tms_stime;  /* dur�e syst�me              */
//                clock_t tms_cutime; /* dur�e utilisateur des fils */
//                clock_t tms_cstime; /* dur�e syst�me des fils     */
// sysconf(_SC_CLK_TCK);
