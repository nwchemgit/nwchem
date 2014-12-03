c     $Id$

      logical function band_stress(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"

*     **** local variables ****
      logical numerical

*     **** external functions ****
      logical  band_minimizer,band_stress_numerical
      external band_minimizer,band_stress_numerical

      if (.not.
     >  btdb_get(rtdb,'nwpw:stress_numerical',mt_log,1,numerical))
     >  numerical = .false.

      if (numerical) then
        band_stress = band_stress_numerical(rtdb)
      else
        band_stress = band_minimizer(rtdb,3)
      end if

      return
      end

