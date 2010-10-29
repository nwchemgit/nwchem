      subroutine argos_diana_getsel
c
c $Id$
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
      if(.not.lsel) then
      call argos_diana_all(int_mb(i_isel),0,dbl_mb(i_wt),nsel)
      lsel=.true.
      endif
      call argos_diana_select(byte_mb(i_snam),int_mb(i_isel),
     + dbl_mb(i_wt),
     + dbl_mb(i_xref),nsel)
      if(me.eq.0) then
      write(*,3003) nsel,nsa
 3003 format(' Selected',i5,' out of',i5,' atoms')
      endif
      return
      end
