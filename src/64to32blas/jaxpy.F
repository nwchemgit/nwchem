      subroutine jaxpy(n,za,zx,incx,zy,incy)
*
* $Id: xaxpy.F 21373 2011-11-02 20:52:29Z d3p852 $
*
      implicit none
#include "y64.fh"
c
      double complex zx(*),zy(*),za
      integer incx,incy,n
c
      INTGR4 incx4,incy4,n4
c
      if(n.le.0)return
      incx4 = incx
      incy4 = incy
      n4 = n
      call zaxpy(n4,za,zx,incx4,zy,incy4)
      return
      end
