      subroutine ana_rtdbin
c
c $Id: ana_rtdbin.F,v 1.1 1998-12-29 18:01:04 d3j191 Exp $
c
      implicit none
c
#include "ana_common.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
c
c     input options from rtdb
c
      if(me.eq.0) then
      if(.not.rtdb_get(irtdb,'ana:irmsd',mt_int,1,irmsd))
     + call errquit('ana_rtdbin: rtdb_get failed',1)
      if(irmsd.eq.1) then
      if(.not.rtdb_cget(irtdb,'ana:rmsfil',1,rmsfil))
     + call errquit('ana_rtdbin: rtdb_get failed',2)
      endif
      endif
c
      if(np.gt.1) then
      call ga_brdcst(mag_d00,irmsd,ma_sizeof(mt_int,1,mt_byte),0)
      endif
c
      return
      end
