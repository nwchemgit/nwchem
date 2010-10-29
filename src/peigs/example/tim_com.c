/*
 $Id$
 *======================================================================
 *
 * DISCLAIMER
 *
 * This material was prepared as an account of work sponsored by an
 * agency of the United States Government.  Neither the United States
 * Government nor the United States Department of Energy, nor Battelle,
 * nor any of their employees, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
 * ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY,
 * COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT,
 * SOFTWARE, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT
 * INFRINGE PRIVATELY OWNED RIGHTS.
 *
 * ACKNOWLEDGMENT
 *
 * This software and its documentation were produced with Government
 * support under Contract Number DE-AC06-76RLO-1830 awarded by the United
 * States Department of Energy.  The Government retains a paid-up
 * non-exclusive, irrevocable worldwide license to reproduce, prepare
 * derivative works, perform publicly and display publicly by or for the
 * Government, including the right to distribute to other Government
 * contractors.
 *
 *======================================================================
 *
 *  -- PEIGS  routine (version 2.1) --
 *     Pacific Northwest Laboratory
 *     July 28, 1995
 *
 *======================================================================
 */
#include <stdio.h>
#include <math.h>

#include "globalp.c.h"

void tim_com(t_total, tcom_me, iwork, work )

   Integer         *iwork;
   DoublePrecision t_total, tcom_me;
   DoublePrecision *work;
/*
 * Compute and print min, max, and mean time spend in mxread/mxwrit
 * and fraction of total cpu time spent in comm
 *
 * t_total = total time for routine being timed (e.g.,for pdspevx or pdspgvx)
 * tcom_me = communication time for this processor
 * iwork   = workspace of length mxnprc_()
 * work    = workspace big enough for mxcombv1 ( i.e., length bufsiz bytes )
 *
 *
 */
{
   Integer                naproc, me, i;
   DoublePrecision        t0, t_min, t_max, t_mean;

   extern void            gmax00(), gsum00();

   extern Integer         mxnprc_(), mxmynd_();

   naproc = mxnprc_();
   me     = mxmynd_();

   for( i=0; i < naproc; i++ )
     iwork[i] = i;

   t0 = tcom_me;

   if( t0 < 0.e0 )
     t0 = 0.e0;

   t_max = t0;
   gmax00( (char *) &t_max, 1, 1, 997, 0, naproc, iwork, work );

   t_min = -t0;
   if( t_min == 0 )
      t_min = -t_max;

   gmax00( (char *) &t_min, 1, 1, 998, 0, naproc, iwork, work );
   t_min = -t_min;

   t_mean = t0;
   gsum00( (char *) &t_mean, 1, 5, 999, 0, naproc, iwork, work );
   t_mean /= naproc;

  if (me == 0 ){
    fprintf(stderr, " Min  communication time = %f \n", t_min);
    fprintf(stderr, " Max  communication time = %f \n", t_max);
    fprintf(stderr, " Mean communication time = %f \n", t_mean);
    fprintf(stderr, " (Max comm)/(total time) = %f \n", t_max/t_total);
  }

   return;
}
