c
c     $Id: pspw_gradient.F,v 1.2 2002-02-27 19:05:21 bylaska Exp $
c
      logical function pspw_gradient(rtdb)
      implicit none
      integer rtdb
      logical  cgsdv5
      external cgsdv5

      pspw_gradient = cgsdv5(rtdb,1)
      return
      end

