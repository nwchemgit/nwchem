      subroutine yswap (n,dx,incx,dy,incy)
*
* $Id$
*
      implicit none
#include "y64.fh"
      double precision dx(*),dy(*)
      integer incx,incy,n
c
      INTGR4 incx4,incy4,n4
c
      if(n.le.0)return
c
      incx4=incx
      incy4=incy
      n4=n
c
      call dswap (n4,dx,incx4,dy,incy4)
      return
      end
