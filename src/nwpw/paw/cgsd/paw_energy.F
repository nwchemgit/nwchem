c     $Id: paw_energy.F,v 1.2 2004-08-01 02:02:44 bylaska Exp $
      logical function paw_energy(rtdb)
      implicit none
      integer rtdb

      logical  paw_minimizer
      external paw_minimizer

      paw_energy = paw_minimizer(rtdb,0)
      return
      end
