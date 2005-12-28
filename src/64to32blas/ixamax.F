      integer function ixamax(n,dx,incx)
*
* $Id: ixamax.F,v 1.1 2005-12-28 18:09:49 edo Exp $
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
      ixamax=izamax(n4,dx,incx4)
      return
      end
