      subroutine  ycopy(n,dx,incx,dy,incy)
*
* $Id: ycopy.F,v 1.1 2005-11-02 17:13:59 edo Exp $
*
      implicit none
#include "y64.fh"
      double precision dx(*),dy(*)
      integer incx,incy,n
c
      INTGR4  n4,incx4,incy4
c
      n4=n
      incx4=incx
      incy4=incy
c
      call  dcopy(n4,dx,incx4,dy,incy4)
      return
      end
