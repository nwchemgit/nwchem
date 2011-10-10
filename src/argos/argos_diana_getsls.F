      subroutine argos_diana_getsls
c
c     get atom selection
c     ------------------
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      if(.not.lsels) then
      call argos_diana_all(int_mb(i_ssel),0,dbl_mb(i_swt),nsel)
      lsels=.true.
      endif
      call argos_diana_select(card,byte_mb(i_snam),int_mb(i_ssel),
     + dbl_mb(i_swt),
     + dbl_mb(i_xref),nsel)
      if(me.eq.0) then
      write(*,3533) nsel,nsa
 3533 format(' Selected',i5,' out of',i5,' atoms for superposition')
      endif
c
      return
      end
c $Id$
