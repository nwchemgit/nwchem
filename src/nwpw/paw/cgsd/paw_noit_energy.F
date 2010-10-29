c     $Id$
      logical function paw_noit_energy(rtdb)
      implicit none
      integer rtdb

      logical  paw_minimizer
      external paw_minimizer

      paw_noit_energy = paw_minimizer(rtdb,-1)
      return
      end
