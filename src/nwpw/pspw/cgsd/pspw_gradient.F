c
c     $Id$
c
      logical function pspw_gradient(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"

*     *** local variables ****
      logical includestress,value,numerical

*     **** external functions ****
      logical  cgsdv5,pspw_stress,pspw_force_numerical
      external cgsdv5,pspw_stress,pspw_force_numerical

      call pspw_cutoff_boot_psi(rtdb)

      if (.not.btdb_get(rtdb,'includestress',mt_log,1,includestress)) 
     >   includestress = .false.

*     **** only calculate atom gradients ****
      if (.not.includestress) then
          if (.not.
     >      btdb_get(rtdb,'nwpw:force_numerical',mt_log,1,numerical))
     >      numerical = .false.
         if (numerical) then
            value = pspw_force_numerical(rtdb)
         else
            value = cgsdv5(rtdb,1)
         end if

*     **** include stress calculation ****
      else
         value = pspw_stress(rtdb)
      end if

      pspw_gradient = value
      return
      end

