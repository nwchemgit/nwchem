      integer function ijamax(n,dx,incx)
*
* $Id: ixamax.F 19697 2010-10-29 16:57:34Z d3y133 $
*
      implicit none
#include "y64.fh"
      double complex dx(*)
      integer incx,n
c
      INTGR4 incx4,n4
      external izamax
      INTGR4 izamax
c
      incx4=incx
      n4=n
c
      ijamax=izamax(n4,dx,incx4)
      return
      end
