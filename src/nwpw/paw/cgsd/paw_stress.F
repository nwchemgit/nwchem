c     $Id: paw_stress.F,v 1.1 2004-11-29 16:05:19 bylaska Exp $

      logical function paw_stress(rtdb)
      implicit none
      integer rtdb

      logical  paw_stress_numerical
      external paw_stress_numerical

      paw_stress = paw_stress_numerical(rtdb)
      return
      end
