c
c     $Id$
c

      logical function band_gradient(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"

*     *** local variables ****
      logical includestress,value

*     **** external functions ****
      logical  band_minimizer,band_stress
      external band_minimizer,band_stress


      if (.not.btdb_get(rtdb,'includestress',mt_log,1,includestress))
     >   includestress = .false.

*     **** only calculate atom gradients ****
      if (.not.includestress) then
         value = band_minimizer(rtdb,1)
      else
        value = band_stress(rtdb)
      end if

      band_gradient = value
      return
      end

