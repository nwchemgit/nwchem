c     $Id: paw_energy.F,v 1.1 2004-07-31 21:54:30 bylaska Exp $
      logical function paw_energy(rtdb)
      implicit none
      integer rtdb

      !logical  paw_cgsd
      !external paw_cgsd

      !paw_energy = paw_cgsd(rtdb,0)
      paw_energy = .true.
      return
      end
