      subroutine qhop_final
c
c $Id: qhop_final.F,v 1.1 2003-08-11 19:43:18 d3j191 Exp $
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
