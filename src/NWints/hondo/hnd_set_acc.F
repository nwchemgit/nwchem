c $Id$
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     hondo : Set HONDO integral tolerance based on that of NWChem
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hnd_set_acc(pnl_tol)
      implicit none
#include "hnd_tol.fh"
      double precision pnl_tol, temp
      integer i
c
c     pnl_tol comes in as a double precision
c     and we need to get the exponent 10**(-itol)
c     temp way of getting the necessary info
c
      i = 0
      temp = pnl_tol
      itol = 15
  100 continue
      if ((temp.le.10.0).and.(temp.ge.1.0)) then
        itol = i
      else
        i = i + 1
        temp = temp * 10.0
        go to 100
      endif
      return
      end
