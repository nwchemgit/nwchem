      subroutine cf_finish()
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer id
c
c     deallocate memory
c
      if(ipme.gt.0) call pme_finish()
c
      call cf_free()
c
      do 1 id=nbs,1,-1
c
      if(id.eq.2.and.msa.gt.0) then
      if(.not.ma_pop_stack(l_wsm))
     + call errquit('Failed to deallocate wsm',id)
      endif
c
      if(id.eq.1.and.mwa.gt.0) then
      if(.not.ma_pop_stack(l_iwq))
     + call errquit('Failed to deallocate iwq',id)
      if(.not.ma_pop_stack(l_iwa))
     + call errquit('Failed to deallocate iwa',id)
      endif
c
      if(.not.ma_pop_stack(l_ixcl(id)))
     + call errquit('Failed to deallocate ixcl',id)
      if(.not.ma_pop_stack(l_itrd(id)))
     + call errquit('Failed to deallocate itrd',id)
c
      if(.not.ma_pop_stack(l_rimp(id)))
     + call errquit('Failed to deallocate rimp',id)
      if(.not.ma_pop_stack(l_iimp(id)))
     + call errquit('Failed to deallocate iimp',id)
      if(.not.ma_pop_stack(l_imp(id)))
     + call errquit('Failed to deallocate imp',id)
c
      if(.not.ma_pop_stack(l_rdih(id)))
     + call errquit('Failed to deallocate rdih',id)
      if(.not.ma_pop_stack(l_idih(id)))
     + call errquit('Failed to deallocate idih',id)
      if(.not.ma_pop_stack(l_dih(id)))
     + call errquit('Failed to deallocate dih',id)
c
      if(.not.ma_pop_stack(l_rang(id)))
     + call errquit('Failed to deallocate rang',id)
      if(.not.ma_pop_stack(l_iang(id)))
     + call errquit('Failed to deallocate iang',id)
      if(.not.ma_pop_stack(l_ang(id)))
     + call errquit('Failed to deallocate ang',id)
c
      if(.not.ma_pop_stack(l_rbnd(id)))
     + call errquit('Failed to deallocate rbnd',id)
      if(.not.ma_pop_stack(l_ibnd(id)))
     + call errquit('Failed to deallocate ibnd',id)
      if(.not.ma_pop_stack(l_bnd(id)))
     + call errquit('Failed to deallocate bnd',id)
c
    1 continue
c
      if(.not.ma_pop_stack(l_chg))
     + call errquit('Failed to deallocate chg',0)
c
      if(.not.ma_pop_stack(l_ias))
     + call errquit('Failed to deallocate ias',0)
c
      if(.not.ma_pop_stack(l_num))
     + call errquit('Failed to deallocate num',0)
c
      if(.not.ma_pop_stack(l_nam))
     + call errquit('Failed to deallocate nam',0)
c
      if(.not.ma_pop_stack(l_mas))
     + call errquit('Failed to deallocate mas',0)
c
      if(.not.ma_pop_stack(l_vdw))
     + call errquit('Failed to deallocate vdw',0)
c
      return
      end
