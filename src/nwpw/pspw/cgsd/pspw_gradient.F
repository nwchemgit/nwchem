c     $Id: pspw_gradient.F,v 1.1 2001-08-30 00:25:11 edo Exp $
      logical function pspw_gradient(rtdb)
      implicit none
      integer rtdb
      logical  cgsdv5
      external cgsdv5

      pspw_gradient = cgsdv5(rtdb,1)
      return
      end

