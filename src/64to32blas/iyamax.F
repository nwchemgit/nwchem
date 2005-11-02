      integer function iyamax(n,dx,incx)
*
* $Id: iyamax.F,v 1.1 2005-11-02 17:13:58 edo Exp $
*
      implicit none
#include "y64.fh"
      double precision dx(*)
      integer incx,n
c
      INTGR4 incx4,n4
      external idamax
      INTGR4 idamax
c
      incx4=incx
      n4=n
c
      iyamax=idamax(n4,dx,incx4)
      return
      end
