      subroutine yLARNV( IDIST, ISEED, N, X )
c
* $Id$
c
      implicit none
#include "y64.fh"
      INTEGER            IDIST, N
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      DOUBLE PRECISION   X( * )
c
      INTGR4 idist4,n4,iseed4(4)
c
      idist4=idist
      n4=n
      iseed4(1)=iseed(1)
      iseed4(2)=iseed(2)
      iseed4(3)=iseed(3)
      iseed4(4)=iseed(4)
c     
      call DLARNV( IDIST4, ISEED4, N4, X )
      iseed(1)=iseed4(1)
      iseed(2)=iseed4(2)
      iseed(3)=iseed4(3)
      iseed(4)=iseed4(4)
      return
      end
