      subroutine argos_diana_getslo
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
      if(.not.lselo) then
      call argos_diana_all(int_mb(i_osel),0,dbl_mb(i_owt),nselo)
      lselo=.true.
      endif
      call argos_diana_select(card,byte_mb(i_snam),int_mb(i_osel),
     + dbl_mb(i_owt),
     + dbl_mb(i_xref),nselo)
      if(me.eq.0) then
      write(*,3633) nselo,nsa
 3633 format(' Selected',i5,' out of',i5,' atoms for output')
      endif
c
      return
      end
c $Id$
