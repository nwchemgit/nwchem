      subroutine xc_exit_index(l_cntoce, l_cntobfr, l_cetobfr,
     &                         l_rdens_atom)
c
c     $Id$
c
c     Cleans up indexing and reduced density arrays used by several
c     XC routines
c
c     BGJ - 8/98
c
      implicit none
#include "errquit.fh"
      integer l_cntoce, l_cntobfr, l_cetobfr, l_rdens_atom
c
#include "mafdecls.fh"
c
      if (.not.ma_pop_stack(l_rdens_atom))
     &   call errquit('xc_exit_index: cannot pop stack',0, MA_ERR)
      if (.not.ma_pop_stack(l_cetobfr))
     &   call errquit('xc_exit_index: cannot pop stack',0, MA_ERR)
      if (.not.ma_pop_stack(l_cntobfr))
     &   call errquit('xc_exit_index: cannot pop stack',0, MA_ERR)
      if (.not.ma_pop_stack(l_cntoce))
     &   call errquit('xc_exit_index: cannot pop stack',0, MA_ERR)
c
      return
      end
