c
c     $Id: pspw_gradient.F,v 1.3 2002-05-04 19:28:30 bylaska Exp $
c
      logical function pspw_gradient(rtdb)
      implicit none
      integer rtdb
      logical  cgsdv5
      external cgsdv5

#include "mafdecls.fh"
#include "rtdb.fh"

*     *** local variables ****
      logical includestress,value

      if (.not.rtdb_get(rtdb,'includestress',mt_log,1,includestress)) 
     >   includestress = .false.

*     **** only calculate atom gradients ****
      if (.not.includestress) then
         value = cgsdv5(rtdb,1)

*     **** include stress calculation ****
      else
         value = cgsdv5(rtdb,3)
      end if

      pspw_gradient = value
      return
      end

