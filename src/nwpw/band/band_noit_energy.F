c     $Id: band_noit_energy.F,v 1.1 2004-11-13 02:25:29 bylaska Exp $

      logical function band_noit_energy(rtdb)
      implicit none
      integer rtdb
      logical  band_minimizer
      external band_minimizer

      band_noit_energy = band_minimizer(rtdb,-1)
      return
      end

