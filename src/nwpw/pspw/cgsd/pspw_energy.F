c     $Id: pspw_energy.F,v 1.1 2001-08-30 00:25:11 edo Exp $
      logical function pspw_energy(rtdb)
      implicit none
      integer rtdb

      logical  cgsdv5
      external cgsdv5

      pspw_energy = cgsdv5(rtdb,0)
      return
      end

