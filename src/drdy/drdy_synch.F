      subroutine drdy_synch(opcode,calling_routine)
      implicit none
#include "drdyP.fh"
#include "nwc_drdyP.fh"
#include "tcgmsg.fh"
#include "msgids.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "printlevels.fh"
*::functions::
      logical  task_energy, task_gradient, task_hessian
      external task_energy, task_gradient, task_hessian
*::passed::
      integer opcode
      character*(*)calling_routine
*::local::
      character*255 drdy_theory
      integer current_print_level
      logical status
      logical status_rtdb, status_ignore
*
*      write(6,*)' drdy_synch called with code ',opcode,
*     &    ' from ',calling_routine,' node =',me_drdy
*      call util_flush(6)
      if (me_drdy.eq.master) then
        if (.not.
     &        rtdb_cput(my_rtdb,'scf:input vectors',1,'atomic')) then
          write(luout,*)' could not reset input vectors to atomic '
          call errquit('drdy_synch: fatal error ',911)          
        endif
        if (.not.
     &        rtdb_cput(my_rtdb,'dft:input vectors',1,'atomic')) then
          write(luout,*)' could not reset input vectors to atomic '
          call errquit('drdy_synch: fatal error ',911)          
        endif
        if (.not.
     &        rtdb_cput(my_rtdb,'mp2:input vectors',1,'atomic')) then
          write(luout,*)' could not reset input vectors to atomic '
          call errquit('drdy_synch: fatal error ',911)          
        endif
        if (.not.
     &        rtdb_cput(my_rtdb,'mcscf:input vectors',1,'atomic')) then
          write(luout,*)' could not reset input vectors to atomic '
          call errquit('drdy_synch: fatal error ',911)          
        endif
        if (.not.rtdb_cget(my_rtdb,'drdy:theory',1,drdy_theory)) then
          write(luout,*)' could not determine drdy:theory'
          call errquit('drdy_synch: fatal error ',911)
        endif
        if (.not.rtdb_cput(my_rtdb,'task:theory',1,drdy_theory)) then
          write(luout,*)' could not write drdy:theory -> theory'
          call errquit('drdy_synch: fatal error ',911)
        endif
      endif
00001 continue
*      write(6,*)' B4 msg_drdy_opcode = ',msg_drdy_opcode,
*     &    ' opcode = ',opcode,' node=',me_drdy
*      call util_flush(6)
      call ga_brdcst(msg_DRDY_OPCODE, opcode, mitob(1), master)
*      write(6,*)' AF msg_drdy_opcode = ',msg_drdy_opcode,
*     &    ' opcode = ',opcode,' node=',me_drdy
*      call util_flush(6)
*     call util_print_common_data
      call util_print_get_level(current_print_level)
      call util_print_set_level(print_none)
*     call util_print_common_data
      if (opcode .eq. DRDY_CODE_ENERGY) then
        status_rtdb   = rtdb_parallel(.true.)
*        call setdbg(0)
        status        = task_energy(my_rtdb)
*        call setdbg(1)
        status_ignore = rtdb_parallel(status_rtdb)
      elseif (opcode .eq. DRDY_CODE_GRADIENT) then
        status_rtdb   = rtdb_parallel(.true.)
*        call setdbg(0)
        status        = task_gradient(my_rtdb)
*        call setdbg(1)
        status_ignore = rtdb_parallel(status_rtdb)
      elseif (opcode .eq. DRDY_CODE_HESSIAN) then
        status_rtdb   = rtdb_parallel(.true.)
*        call setdbg(0)
        status        = task_hessian(my_rtdb)
        call ga_sync()
        status        = status.and.task_gradient(my_rtdb)
*        call setdbg(1)
        status_ignore = rtdb_parallel(status_rtdb)
      elseif (opcode .eq. DRDY_CODE_DONE) then
        call ga_sync()
        return
      else
        write(luout,*)' unknown opcode ',opcode
        call errquit('drdy_synch: fatal error ',911)
      endif
      call util_print_set_level(current_print_level)
      call ga_sync()
      if (me_drdy.eq.master) then
        return
      endif
      opcode = 0
      goto 00001
      end
