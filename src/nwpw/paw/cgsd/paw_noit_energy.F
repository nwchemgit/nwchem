c     $Id: paw_noit_energy.F,v 1.1 2004-11-29 16:05:19 bylaska Exp $
      logical function paw_noit_energy(rtdb)
      implicit none
      integer rtdb

      logical  paw_minimizer
      external paw_minimizer

      paw_noit_energy = paw_minimizer(rtdb,-1)
      return
      end
