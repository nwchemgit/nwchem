c
c     $Id: pspw_gradient.F,v 1.4 2004-11-29 16:05:21 bylaska Exp $
c
      logical function pspw_gradient(rtdb)
      implicit none
      integer rtdb

#include "mafdecls.fh"
#include "rtdb.fh"

*     *** local variables ****
      logical includestress,value

*     **** external functions ****
      logical  cgsdv5,pspw_stress
      external cgsdv5,pspw_stress

      if (.not.rtdb_get(rtdb,'includestress',mt_log,1,includestress)) 
     >   includestress = .false.

*     **** only calculate atom gradients ****
      if (.not.includestress) then
         value = cgsdv5(rtdb,1)

*     **** include stress calculation ****
      else
         value = pspw_stress(rtdb)
      end if

      pspw_gradient = value
      return
      end

