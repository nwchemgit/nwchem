      subroutine schwarz_tidy()
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
