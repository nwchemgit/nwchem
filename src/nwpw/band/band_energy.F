c     $Id$
      logical function band_energy(rtdb)
      implicit none
      integer rtdb
      logical  band_minimizer
      external band_minimizer

      band_energy = band_minimizer(rtdb,0)
      return
      end

