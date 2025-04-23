*
* $Id$
*
c
c     Communications buffer.  bufsiz is in bytes.  As with 'data',
c     (see mxcomb.cpp) 'bufdat' really ought to be integer*4.
c
         integer bufsiz
#ifndef NCUBE_NODE
         parameter (bufsiz=10000)
#else
c      parameter (bufsiz=512)
      parameter (bufsiz=8*1024)
#endif
c     Place datbuf in common to force alignment.  The djunk variable
c     is left over from experiments on how alignment changes
c     performance.  (On an iPSC/860, making djunk type integer*4
c     degrades performance of the combine by a factor of 5X.)
C*
C*    Pass datbuf as an INTEGER array argument.
C*
C*      common /cmbbfr/ djunk, datbuf
C*      double precision djunk
C*      integer datbuf(bufsiz/4)
