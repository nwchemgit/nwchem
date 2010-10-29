c     $Id$

      logical function band_noit_energy(rtdb)
      implicit none
      integer rtdb
      logical  band_minimizer
      external band_minimizer

      band_noit_energy = band_minimizer(rtdb,-1)
      return
      end

