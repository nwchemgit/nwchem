      subroutine  xscal(n,za,zx,incx)
*
* $Id$
*
      implicit none
#include "y64.fh"
      double complex za,zx(*)
      integer incx,n
c
      INTGR4 incx4,n4
c
      if( n.le.0 .or. incx.le.0 )return
c
      n4=n
      incx4=incx
c
      call zscal(n4,za,zx,incx4)
      return
      end
