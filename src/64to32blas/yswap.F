      subroutine yswap (n,dx,incx,dy,incy)
*
* $Id: yswap.F,v 1.1 2005-11-02 17:14:02 edo Exp $
*
      implicit none
#include "y64.fh"
      double precision dx(*),dy(*),dtemp
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
