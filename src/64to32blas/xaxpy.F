      subroutine xaxpy(n,za,zx,incx,zy,incy)
*
* $Id: xaxpy.F,v 1.1 2005-11-02 17:13:58 edo Exp $
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
      call zaxpy(n4,za,zx,incx4,zy,incy4)
      return
      end
