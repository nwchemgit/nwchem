c     $Id: pspw_noit_energy.F,v 1.1 2002-04-25 00:34:39 bylaska Exp $
      logical function pspw_noit_energy(rtdb)
      implicit none
      integer rtdb

      logical  cgsdv5
      external cgsdv5

      pspw_noit_energy = cgsdv5(rtdb,-1)
      return
      end

