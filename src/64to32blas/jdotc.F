      double complex function jdotc(n,dx,incx,dy,incy)
*
* $Id: xdotc.F 19697 2010-10-29 16:57:34Z d3y133 $
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
      jdotc= zdotc (n4,dx,incx4,dy,incy4)
      return
      end
