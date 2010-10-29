c     $Id$
      logical function pspw_noit_energy(rtdb)
      implicit none
      integer rtdb

      logical  cgsdv5
      external cgsdv5

      pspw_noit_energy = cgsdv5(rtdb,-1)
      return
      end

