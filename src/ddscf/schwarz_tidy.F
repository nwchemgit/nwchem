      subroutine schwarz_tidy()
C$Id: schwarz_tidy.F,v 1.3 1995-02-02 23:12:15 d3g681 Exp $
      implicit none
#include "cschwarz.fh"
#include "mafdecls.fh"
c
      if (.not. MA_free_heap(l_cat))
     $     call errquit('schwarz_tidy: ma_free_heap failed', 0)
      if (.not. MA_free_heap(l_csh))
     $     call errquit('schwarz_tidy: ma_free_heap failed', 0)
c
      end
