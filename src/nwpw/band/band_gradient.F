c
c     $Id: band_gradient.F,v 1.2 2002-01-20 23:49:20 bylaska Exp $
c

      logical function band_gradient(rtdb)
      implicit none
      integer rtdb
      logical  band_minimizer
      external band_minimizer

      band_gradient = band_minimizer(rtdb,1)
      return
      end

