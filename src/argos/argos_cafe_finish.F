      subroutine argos_cafe_finish()
c $Id$
      implicit none
c
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      integer id
c
c     deallocate memory
c
      if(nhop.gt.0) call qhop_finish()
c
      if(ipme.gt.0) call argos_pme_finish()
c
      if(numpmf.gt.0) then
      if(.not.ma_pop_stack(l_upmf))
     + call md_abort('Failed to deallocate upmf',0)
      if(.not.ma_pop_stack(l_wpmf))
     + call md_abort('Failed to deallocate wpmf',0)
      if(.not.ma_pop_stack(l_ypmf))
     + call md_abort('Failed to deallocate ypmf',0)
      if(.not.ma_pop_stack(l_xpmf))
     + call md_abort('Failed to deallocate xpmf',0)
      if(.not.ma_pop_stack(l_rpmf))
     + call md_abort('Failed to deallocate rpmf',0)
      if(.not.ma_pop_stack(l_jpmf))
     + call md_abort('Failed to deallocate jpmf',0)
      if(.not.ma_pop_stack(l_ipmf))
     + call md_abort('Failed to deallocate ipmf',0)
      endif
c
      if(nhop.gt.0) then
      if(.not.ma_pop_stack(l_timhop))
     + call md_abort('Failed to deallocate timhop',0)
      if(.not.ma_pop_stack(l_lsthop))
     + call md_abort('Failed to deallocate lsthop',0)
      if(.not.ma_pop_stack(l_pda))
     + call md_abort('Failed to deallocate pda',0)
      if(.not.ma_pop_stack(l_uda))
     + call md_abort('Failed to deallocate uda',0)
      if(.not.ma_pop_stack(l_rda))
     + call md_abort('Failed to deallocate rda',0)
      if(.not.ma_pop_stack(l_lda))
     + call md_abort('Failed to deallocate lda',0)
      endif
      if(ndrs.gt.0) then
      if(.not.ma_pop_stack(l_xdrs))
     + call md_abort('Failed to deallocate xdrs',0)
      if(.not.ma_pop_stack(l_rdrs))
     + call md_abort('Failed to deallocate rdrs',0)
      if(.not.ma_pop_stack(l_idrs))
     + call md_abort('Failed to deallocate idrs',0)
      endif
c
      call argos_cafe_free()
c
      if(mmult.gt.0) then
      if(.not.ma_pop_stack(l_fmul))
     + call md_abort('Failed to deallocate fmul',0)
      if(.not.ma_pop_stack(l_xmul))
     + call md_abort('Failed to deallocate xmul',0)
      if(.not.ma_pop_stack(l_imul))
     + call md_abort('Failed to deallocate imul',0)
      if(.not.ma_pop_stack(l_ixmul))
     + call md_abort('Failed to deallocate ixmul',0)
      endif
c
      do 1 id=nbs,1,-1
c
      if(id.eq.2.and.msa.gt.0) then
      if(.not.ma_pop_stack(l_wsm))
     + call md_abort('Failed to deallocate wsm',id)
      endif
c
      if(id.eq.1.and.mwa.gt.0) then
      if(.not.ma_pop_stack(l_iwq))
     + call md_abort('Failed to deallocate iwq',id)
      if(.not.ma_pop_stack(l_iwa))
     + call md_abort('Failed to deallocate iwa',id)
      endif
c
      if(.not.ma_pop_stack(l_ixcl(id)))
     + call md_abort('Failed to deallocate ixcl',id)
      if(.not.ma_pop_stack(l_itrd(id)))
     + call md_abort('Failed to deallocate itrd',id)
c
      if(.not.ma_pop_stack(l_rimp(id)))
     + call md_abort('Failed to deallocate rimp',id)
      if(.not.ma_pop_stack(l_iimp(id)))
     + call md_abort('Failed to deallocate iimp',id)
      if(.not.ma_pop_stack(l_imp(id)))
     + call md_abort('Failed to deallocate imp',id)
c
      if(.not.ma_pop_stack(l_rdih(id)))
     + call md_abort('Failed to deallocate rdih',id)
      if(.not.ma_pop_stack(l_idih(id)))
     + call md_abort('Failed to deallocate idih',id)
      if(.not.ma_pop_stack(l_dih(id)))
     + call md_abort('Failed to deallocate dih',id)
c
      if(.not.ma_pop_stack(l_rub(id)))
     + call md_abort('Failed to deallocate rub',id)
c
      if(.not.ma_pop_stack(l_rang(id)))
     + call md_abort('Failed to deallocate rang',id)
      if(.not.ma_pop_stack(l_iang(id)))
     + call md_abort('Failed to deallocate iang',id)
      if(.not.ma_pop_stack(l_ang(id)))
     + call md_abort('Failed to deallocate ang',id)
c
      if(.not.ma_pop_stack(l_rbnd(id)))
     + call md_abort('Failed to deallocate rbnd',id)
      if(.not.ma_pop_stack(l_ibnd(id)))
     + call md_abort('Failed to deallocate ibnd',id)
      if(.not.ma_pop_stack(l_bnd(id)))
     + call md_abort('Failed to deallocate bnd',id)
c
    1 continue
c
      if(.not.ma_pop_stack(l_chg))
     + call md_abort('Failed to deallocate chg',0)
c
      if(.not.ma_pop_stack(l_ias))
     + call md_abort('Failed to deallocate ias',0)
c
      if(.not.ma_pop_stack(l_typ))
     + call md_abort('Failed to deallocate typ',0)
c
      if(.not.ma_pop_stack(l_num))
     + call md_abort('Failed to deallocate num',0)
c
      if(.not.ma_pop_stack(l_nam))
     + call md_abort('Failed to deallocate nam',0)
c
      if(.not.ma_pop_stack(l_mas))
     + call md_abort('Failed to deallocate mas',0)
c
      if(.not.ma_pop_stack(l_vdw))
     + call md_abort('Failed to deallocate vdw',0)
c
      if(.not.ma_pop_stack(l_sti))
     + call md_abort('Failed to deallocate sti',0)
c
      if(lqhop) then
      if(.not.ma_pop_stack(l_mprot))
     + call md_abort('Failed to deallocate mprot',0)
      endif
c
      if(.not.ma_pop_stack(l_snam))
     + call md_abort('Failed to deallocate snam',0)
c
c      if(lanal) call ana_finish(byte_mb(i_snam),.false.)
c
      return
      end
c $Id$
