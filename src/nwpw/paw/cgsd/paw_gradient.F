c
c     $Id: paw_gradient.F,v 1.1 2004-07-31 21:54:30 bylaska Exp $
c
      logical function paw_gradient(rtdb)
      implicit none
      integer rtdb
      !logical  paw_cgsd
      !external paw_cgsd

#include "mafdecls.fh"
#include "rtdb.fh"

*     *** local variables ****
      logical includestress,value

      if (.not.rtdb_get(rtdb,'includestress',mt_log,1,includestress))
     >   includestress = .false.

*     **** only calculate atom gradients ****
      if (.not.includestress) then
         !value = paw_cgsd(rtdb,1)

*     **** include stress calculation ****
      else
         !value = paw_cgsd(rtdb,3)
      end if

      value = .true.
      paw_gradient = value
      return
      end

