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
#include <memory.h>
#include <math.h>
#include <stdlib.h>

#define ffabs(a) ((a) >= (0.) ? (a) : (-a))

#include "globalp.c.h"

void pgexit( info, msg, proclist, nprocs, work )
     char            *msg;
     Integer         *info, *proclist, *nprocs;
     DoublePrecision *work;
{
  /*
    This routine does a parallel, global exit if integer *info is non-zero
    on any processor in proclist.
    
    a global combine is performed on abs(*info),
    if info is non-zero upon return, an exit is called with info = -51.
    
    info     = integer to check. On return info = 1 if info
               was non-zero on any processor in proclist
    nprocs   = number of processor ids in proclist
    proclist = array of processor ids on which to check n,
               there must be no repeated processor ids in proclist.
               proclist must be identical on all processors in proclist.

    msg =  message to print before exiting if n <> 0.
    work = workspace or length at least bufsiz bytes (see cmbbrf1.h)

    
    */
  
  static Integer TYPE = 10;

  Integer m;
  
  extern void    gi_sum();
  extern void    xstop_();
  extern Integer mxmynd_();

  m = ffabs( *info );
  
  *info = 0;
  
  if ( *nprocs > 1 )
    gi_sum( &m, 1, TYPE, proclist[0], *nprocs, proclist, work);
  
  if ( m != 0 ) {
    *info = -51;
    fprintf( stderr, " %s  me = %d exiting via pgexit. \n", msg, mxmynd_() );
    xstop_( info );
  }

  return;
}
