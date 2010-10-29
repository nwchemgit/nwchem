      subroutine  yscal(n,da,dx,incx)
*
* $Id$
*
c
      implicit none
#include "y64.fh"
      double precision da,dx(*)
      integer n,incx
c
      INTGR4 n4,incx4
c
      n4=n
      incx4=incx
c
      call dscal(n4,da,dx,incx4)
      return
      end
