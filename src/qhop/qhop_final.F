      subroutine qhop_alloc
c
c $Id: qhop_final.F,v 1.6 2003-08-18 17:35:32 d3j191 Exp $
c
      implicit none
c
#include "qhop_common.fh"
#include "mafdecls.fh"
c
c     original dimensions
c
      maxseq=5000
      mxseq=30
      maxpar=10000
      mxar=mxseq*30
      maxhv=4
c
c     allocate memory for parameter arrays
c
      if(.not.ma_push_get(mt_int,2*maxseq,'iptseq',l_ptseq,i_ptseq))
     + call errquit('Failed to allocate iptseq',0)
      if(.not.ma_push_get(mt_int,mxseq,'iptarat',l_ptarat,i_ptarat))
     + call errquit('Failed to allocate iptarat',0)
      if(.not.ma_push_get(mt_int,4*mxar,'iarat',l_iarat,i_iarat))
     + call errquit('Failed to allocate iarat',0)
      if(.not.ma_push_get(mt_dbl,2*mxar,'racs',l_racs,i_racs))
     + call errquit('Failed to allocate racs',0)
      if(.not.ma_push_get(mt_int,maxpar,'iptpar',l_ptpar,i_ptpar))
     + call errquit('Failed to allocate iptpar',0)
      if(.not.ma_push_get(mt_dbl,31*maxpar,'par',l_par,i_par))
     + call errquit('Failed to allocate par',0)
c
      return
      end
      subroutine qhop_finish
c
      implicit none
c
#include "qhop_common.fh"
#include "mafdecls.fh"
c
      if(.not.ma_pop_stack(l_par))
     + call errquit('qhop_final: Failed to deallocate par',0)
      if(.not.ma_pop_stack(l_ptpar))
     + call errquit('qhop_final: Failed to deallocate ptpar',0)
      if(.not.ma_pop_stack(l_racs))
     + call errquit('qhop_final: Failed to deallocate racs',0)
      if(.not.ma_pop_stack(l_iarat))
     + call errquit('qhop_final: Failed to deallocate arat',0)
      if(.not.ma_pop_stack(l_ptarat))
     + call errquit('qhop_final: Failed to deallocate ptarat',0)
      if(.not.ma_pop_stack(l_ptseq))
     + call errquit('qhop_final: Failed to deallocate ptseq',0)
c
      return
      end
