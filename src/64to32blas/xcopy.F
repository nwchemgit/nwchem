      subroutine  xcopy(n,dx,incx,dy,incy)
*
* $Id$
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
