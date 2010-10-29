c     $Id$
      logical function paw_energy(rtdb)
      implicit none
      integer rtdb

      logical  paw_minimizer
      external paw_minimizer

      paw_energy = paw_minimizer(rtdb,0)
      return
      end
