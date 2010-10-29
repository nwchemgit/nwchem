      subroutine dia_getsel
c
c $Id$
c
c     get atom selection
c     ------------------
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      if(.not.lsel) then
      call dia_all(int_mb(i_isel),0,dbl_mb(i_wt),nsel)
      lsel=.true.
      endif
      call dia_select(byte_mb(i_snam),int_mb(i_isel),dbl_mb(i_wt),
     + dbl_mb(i_xref),nsel)
      if(me.eq.0) then
      write(*,3003) nsel,nsa
 3003 format(' Selected',i5,' out of',i5,' atoms')
      endif
      return
      end
      subroutine dia_getsls
c
c     get atom selection
c     ------------------
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      if(.not.lsels) then
      call dia_all(int_mb(i_ssel),0,dbl_mb(i_swt),nsel)
      lsels=.true.
      endif
      call dia_select(card,byte_mb(i_snam),int_mb(i_ssel),dbl_mb(i_swt),
     + dbl_mb(i_xref),nsel)
      if(me.eq.0) then
      write(*,3533) nsel,nsa
 3533 format(' Selected',i5,' out of',i5,' atoms for superposition')
      endif
c
      return
      end
      subroutine dia_getslo
c
c     get atom selection
c     ------------------
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      if(.not.lselo) then
      call dia_all(int_mb(i_osel),0,dbl_mb(i_owt),nselo)
      lselo=.true.
      endif
      call dia_select(card,byte_mb(i_snam),int_mb(i_osel),dbl_mb(i_owt),
     + dbl_mb(i_xref),nselo)
      if(me.eq.0) then
      write(*,3633) nselo,nsa
 3633 format(' Selected',i5,' out of',i5,' atoms for output')
      endif
c
      return
      end
