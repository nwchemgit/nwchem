c
c     $Id: band_gradient.F,v 1.3 2004-11-29 16:05:16 bylaska Exp $
c

      logical function band_gradient(rtdb)
      implicit none
      integer rtdb

#include "mafdecls.fh"
#include "rtdb.fh"

*     *** local variables ****
      logical includestress,value

*     **** external functions ****
      logical  band_minimizer,band_stress
      external band_minimizer,band_stress


      if (.not.rtdb_get(rtdb,'includestress',mt_log,1,includestress))
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

