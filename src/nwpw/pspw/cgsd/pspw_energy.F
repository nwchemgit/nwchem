c     $Id$
      logical function pspw_energy(rtdb)
      implicit none
      integer rtdb

      logical  cgsdv5
      external cgsdv5

      pspw_energy = cgsdv5(rtdb,0)
      return
      end

