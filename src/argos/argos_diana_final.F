      subroutine argos_diana_final()
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "util.fh"
#include "rtdb.fh"
c
      logical lopen
c
      if(me.eq.0) close(unit=lfncmd,status='delete')
c      if(me.eq.0) close(unit=lfncmd,status='keep')
c
      if(active) call argos_diana_edfinal()
c
      if(.not.ma_pop_stack(l_sbnd))
     + call md_abort('Could not deallocate sbnd',0)
      if(.not.ma_pop_stack(l_iequi))
     + call md_abort('Could not deallocate iequi',0)
      if(.not.ma_pop_stack(l_imol))
     + call md_abort('Could not deallocate imol',0)
      if(.not.ma_pop_stack(l_iram))
     + call md_abort('Could not deallocate iram',0)
c
      call argos_diana_finish(byte_mb(i_snam),.true.)
c
      if(.not.ma_pop_stack(l_wrk))
     + call md_abort('Could not deallocate iwrk',0)
c
      inquire(unit=lfnpca,opened=lopen)
      if(lopen) close(lfnpca,status='keep')
      inquire(unit=lfnord,opened=lopen)
      if(lopen) close(lfnord,status='keep')
      inquire(unit=lfndat,opened=lopen)
      if(lopen) close(lfndat,status='keep')
c
      ignore=rtdb_parallel(.true.)
c
      return
      end
c $Id$
