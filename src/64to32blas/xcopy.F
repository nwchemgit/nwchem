      subroutine  xcopy(n,dx,incx,dy,incy)
*
* $Id: xcopy.F,v 1.1 2005-12-28 18:09:49 edo Exp $
*
      implicit none
#include "y64.fh"
      double complex dx(*),dy(*)
      integer incx,incy,n
c
      INTGR4  n4,incx4,incy4
c
      n4=n
      incx4=incx
      incy4=incy
c
      call  zcopy(n4,dx,incx4,dy,incy4)
      return
      end
