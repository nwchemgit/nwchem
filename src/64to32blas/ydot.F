      double precision function ydot(n,dx,incx,dy,incy)
*
* $Id: ydot.F,v 1.1 2005-11-02 17:13:59 edo Exp $
*
      implicit none
#include "y64.fh"
      double precision dx(*),dy(*)
      integer incx,incy,n
c
      INTGR4 incx4,incy4,n4
      external ddot
      double precision ddot
c
      incx4=incx
      incy4=incy
      n4=n
c
      ydot= ddot(n4,dx,incx4,dy,incy4)
      return
      end
