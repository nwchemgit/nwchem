      subroutine tce_tidy(rtdb)
c
c $Id: tce_tidy.F,v 1.1 2002-10-18 23:42:09 sohirata Exp $
c
c     Terminate TCE modulem.
c
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "tce.fh"
c
      integer rtdb             ! RTDB handle
      integer i
c
      do i=ipol,1,-1
        if (.not.ma_pop_stack(l_irs(i))) call errquit
     1    ('tce_tidy: MA problem',0)
      enddo
      do i=ipol,1,-1
        if (.not.ma_pop_stack(l_occ(i))) call errquit
     1    ('tce_tidy: MA problem',0)
        if (.not.ma_pop_stack(l_evl(i))) call errquit
     1    ('tce_tidy: MA problem',0)
      enddo
      if (.not.geom_destroy(geom))
     1  call errquit('tce_tidy: geom_destroy failed',0)
      if (.not.bas_destroy(ao_bas_han))
     1  call errquit('tce_tidy: bas_destroy failed',0)
      call fock_2e_tidy(rtdb)
      call int_terminate()
      call schwarz_tidy()
      return
      end
