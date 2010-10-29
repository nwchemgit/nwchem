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

/*
  PeIGS communication utility; not supported
  Performs global max on vectors
  */

void gmax00(buf, items, datatype, msgtype, root, snumprocs, plist, work)
     /*
       this implementation of gmax00
       leaves the result in all nodes participating nodes
       (0..snumprocs-1) instead of just the root.
       
       Performs global max on DoublePrecision precision scalars or vectors
       
       */
     char *buf;    /* scalar or vector to be combined */
     Integer items;    /* number of such in bytes */
     Integer datatype; /* not used */
     Integer msgtype;  /* message type for communication */
     Integer root;     /* the root for the b-tree global combine */
     Integer *plist, snumprocs;    /* the list of processors, # of processors */
     DoublePrecision *work;  /* workspace containing at least bufsiz bytes (see cmbbrf.h) */
{
  Integer isize;
  extern Integer maxdv_();
  extern void mxcombv1_();
  
  isize = sizeof(DoublePrecision);
  mxcombv1_ ( buf, maxdv_, &isize, &items, &snumprocs, plist, &msgtype, (char *)work);
  
  return;
}

