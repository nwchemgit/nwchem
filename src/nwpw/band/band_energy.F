c     $Id: band_energy.F,v 1.2 2001-11-17 01:01:16 bylaska Exp $
      logical function band_energy(rtdb)
      implicit none
      integer rtdb
      logical  band_minimizer
      external band_minimizer

      band_energy = band_minimizer(rtdb,1)
      return
      end

