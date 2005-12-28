      double complex function xdotc(n,dx,incx,dy,incy)
*
* $Id: xdotc.F,v 1.1 2005-12-28 18:09:49 edo Exp $
*
      implicit none
#include "y64.fh"
      double complex dx(*),dy(*)
      integer incx,incy,n
c
      INTGR4 incx4,incy4,n4
      external zdotc
      double complex zdotc
c
      incx4=incx
      incy4=incy
      n4=n
c
      xdotc= zdotc (n4,dx,incx4,dy,incy4)
      return
      end
