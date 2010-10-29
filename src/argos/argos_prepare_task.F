      logical function argos_prepare_task(irtdb,itask,lfnout)
c
c $Id$
c
      implicit none
c
#include "rtdb.fh"
#include "global.fh"
#include "msgids.fh"
#include "mafdecls.fh"
#include "util.fh"
c
      logical scf,esp
      external scf,esp
c
      integer irtdb,itask,lfnout
      logical lstate
c
      call ga_sync()
      call ga_brdcst(pre_i00,itask,ma_sizeof(mt_int,1,mt_byte),0)
c
      argos_prepare_task=.true.
c
      if(itask.eq.1) then
c
      lstate=rtdb_parallel(.true.)
c
      argos_prepare_task=.false.
c
c      if(.not.rtdb_cput(irtdb,'scf:print',1,'none'))
c     + call md_abort('argos_prepare_task1: rtdb_cput failed',9999)
c
      if(scf(irtdb)) then
c
c      if(.not.rtdb_cput(irtdb,'esp:print',1,'none'))
c     + call md_abort('argos_prepare_task1: rtdb_cput failed',9999)
c
      if(esp(irtdb)) then
      argos_prepare_task=.true.
      endif
      endif
c
      lstate=rtdb_parallel(.false.)
c
      endif
c
      return
      end
