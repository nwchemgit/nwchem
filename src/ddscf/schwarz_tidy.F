      subroutine schwarz_tidy()
C$Id: schwarz_tidy.F,v 1.4 2003-10-17 22:54:33 carlfahl Exp $
      implicit none
#include "errquit.fh"
#include "cschwarz.fh"
#include "mafdecls.fh"
c
      if (.not. MA_free_heap(l_cat))
     $     call errquit('schwarz_tidy: ma_free_heap failed', 0, MA_ERR)
      if (.not. MA_free_heap(l_csh))
     $     call errquit('schwarz_tidy: ma_free_heap failed', 0, MA_ERR)
c
      end
