      subroutine tce_tidy(rtdb)
c
c $Id: tce_tidy.F,v 1.4 2002-11-11 18:25:55 sohirata Exp $
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
      logical nodezero
c
      nodezero = (ga_nodeid().eq.0)
      if (nodezero.and.util_print('time',print_debug)) then
        write(LuOut,9000) 'createfile',cpusecs(1)
        write(LuOut,9000) 'deletefile',cpusecs(2)
        write(LuOut,9000) 'put_block', cpusecs(3)
        write(LuOut,9000) 'add_block', cpusecs(4)
        write(LuOut,9000) 'get_block', cpusecs(5)
      endif
      do i=ipol,1,-1
        if (.not.ma_pop_stack(l_irs(i))) call errquit
     1    ('tce_tidy: MA problem',0)
        if (.not.ma_pop_stack(l_occ(i))) call errquit
     1    ('tce_tidy: MA problem',1)
        if (.not.ma_pop_stack(l_evl(i))) call errquit
     1    ('tce_tidy: MA problem',2)
      enddo
      if (.not.geom_destroy(geom))
     1  call errquit('tce_tidy: geom_destroy failed',0)
      if (.not.bas_destroy(ao_bas_han))
     1  call errquit('tce_tidy: bas_destroy failed',0)
      call fock_2e_tidy(rtdb)
      call int_terminate()
      call schwarz_tidy()
 9000 format(1x,'CPU time for ',A,f20.1)
      return
      end
