      SUBROUTINE yLAGTF( N, A, LAMBDA, B, C, TOL, D, IN, INFO )
c
* $Id$
c
      implicit none
#include "y64.fh"
      INTEGER            INFO, N
      DOUBLE PRECISION   LAMBDA, TOL
*     ..
*     .. Array Arguments ..
      INTEGER            IN( * )
      DOUBLE PRECISION   A( * ), B( * ), C( * ), D( * )
c
      INTGR4            INFO4, N4,IN4(N)
      integer i
c
      n4=n
      call DLAGTF( N4, A, LAMBDA, B, C, TOL, D, IN4, INFO4 )
      do i=1,n
         in(i)=in4(i)
      enddo
      info=info4
      return
      end
