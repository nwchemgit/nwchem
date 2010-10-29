      double precision function yasum(n,dx,incx)
*
* $Id$
*
      implicit none
#include "y64.fh"
c
      double precision dx(*)
      integer incx,n
c
      INTGR4 incx4,n4
      external dasum
      double precision dasum
c
      yasum = 0.0d0
      if( n.le.0 .or. incx.le.0 )return
      incx4=incx
      n4=n
      yasum=dasum(n4,dx,incx4)
      return
      end
