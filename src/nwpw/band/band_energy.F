c     $Id: band_energy.F,v 1.3 2003-12-02 19:16:46 bylaska Exp $
      logical function band_energy(rtdb)
      implicit none
      integer rtdb
      logical  band_minimizer
      external band_minimizer

      band_energy = band_minimizer(rtdb,0)
      return
      end

