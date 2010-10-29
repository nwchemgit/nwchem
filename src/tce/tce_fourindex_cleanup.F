      subroutine tce_fourindex_cleanup(d_f1,l_f1_offset,k_f1_offset)
c
c $Id$
c 
c Written by Jeff Hammond, January 2008.
c
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
c#include "geom.fh"
c#include "rtdb.fh"
c#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "sf.fh"
#include "inp.fh"
#include "errquit.fh"
c#include "tce.fh"
c#include "tce_main.fh"
#include "tce_ints.fh"
      call deletefile(d_f1)
      if (.not.ma_pop_stack(l_f1_offset))
     1  call errquit("tce_fourindex_cleanup: MA problem",
     2               l_f1_offset,MA_ERR)
      return
      end
