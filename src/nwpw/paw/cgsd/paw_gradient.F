c
c     $Id: paw_gradient.F,v 1.3 2004-11-29 16:05:19 bylaska Exp $
c
      logical function paw_gradient(rtdb)
      implicit none
      integer rtdb
      logical  paw_minimizer,paw_stress
      external paw_minimizer,paw_stress

#include "mafdecls.fh"
#include "rtdb.fh"

*     *** local variables ****
      logical includestress,value

      if (.not.rtdb_get(rtdb,'includestress',mt_log,1,includestress))
     >   includestress = .false.

*     **** only calculate atom gradients ****
      if (.not.includestress) then
         value = paw_minimizer(rtdb,1)

*     **** include stress calculation ****
      else
         value = paw_stress(rtdb)
      end if

      paw_gradient = value
      return
      end

