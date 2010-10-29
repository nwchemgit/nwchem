      SUBROUTINE YLAGTS( JOB, N, A, B, C, D, IN, Y, TOL, INFO )
      implicit none
#include "y64.fh"
c
* $Id$
c
      INTEGER            INFO, JOB, N
      DOUBLE PRECISION   TOL
*     ..
*     .. Array Arguments ..
      INTEGER            IN( * )
      DOUBLE PRECISION   A( * ), B( * ), C( * ), D( * ), Y( * )
c
      INTGR4            job4, INFO4, N4,IN4(N)
      integer i
c
      n4=n
      job4=job
      do i=1,n
         in4(i)=in(i)
      enddo
      call DLAGTS( JOB4, N4, A, B, C, D, IN4, Y, TOL, INFO4 )
      info=info4
      return
      end
