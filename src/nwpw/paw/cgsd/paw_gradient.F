c
c     $Id: paw_gradient.F,v 1.2 2004-08-01 02:02:44 bylaska Exp $
c
      logical function paw_gradient(rtdb)
      implicit none
      integer rtdb
      logical  paw_minimizer
      external paw_minimizer

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
         value = paw_minimizer(rtdb,3)
      end if

      paw_gradient = value
      return
      end

