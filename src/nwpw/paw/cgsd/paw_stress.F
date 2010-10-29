c     $Id$

      logical function paw_stress(rtdb)
      implicit none
      integer rtdb

      logical  paw_stress_numerical
      external paw_stress_numerical

      paw_stress = paw_stress_numerical(rtdb)
      return
      end
