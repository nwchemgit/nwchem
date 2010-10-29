      subroutine tce_grad_tidy(rtdb)
c
c $Id$
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
#include "tce_main.fh"
#include "errquit.fh"
#include "dra.fh"
c
      integer rtdb             ! RTDB handle
      integer i
      logical nodezero
c
      nodezero = (ga_nodeid().eq.0)
      if (nodezero.and.util_print('time',print_debug)) then
        write(LuOut,9000) 'createfile   ', cpusecs(1), cpusecs(51)
        write(LuOut,9000) 'deletefile   ', cpusecs(2), cpusecs(52)
        write(LuOut,9000) 'put_block    ', cpusecs(3), cpusecs(53)
        write(LuOut,9000) 'add_block    ', cpusecs(4), cpusecs(54)
        write(LuOut,9000) 'get_block    ', cpusecs(5), cpusecs(55)
        write(LuOut,9000) 'reconcilefile', cpusecs(6), cpusecs(56)
        write(LuOut,9000) 'tce_sort     ', cpusecs(7), cpusecs(57)
      endif
      do i=ipol,1,-1
        if (.not.ma_pop_stack(l_irs(i))) call errquit
     1    ('tce_tidy: MA problem',1,MA_ERR)
        if (.not.ma_pop_stack(l_occ(i))) call errquit
     1    ('tce_tidy: MA problem',2,MA_ERR)
        if (.not.ma_pop_stack(l_evl(i))) call errquit
     1    ('tce_tidy: MA problem',3,MA_ERR)
      enddo
      if (.not.geom_destroy(geom))
     1  call errquit('tce_tidy: geom_destroy failed',0,GEOM_ERR)
      if (.not.bas_destroy(ao_bas_han))
     1  call errquit('tce_tidy: bas_destroy failed',0,BASIS_ERR)
      if (ioalg.eq.5) then
        if (dra_terminate().ne.0)
     1  call errquit('tce_tidy: dra termination failed',0,DISK_ERR)
      endif
      call fock_2e_tidy(rtdb)
c      call int_terminate()
c      call schwarz_tidy()
 9000 format(1x,'cpu & wall time in seconds for ',A,2f20.1)
      return
      end
