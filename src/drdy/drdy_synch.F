C $Id$
      subroutine drdy_synch(opcode,calling_routine)
      implicit none
#include "errquit.fh"
#include "drdyP.fh"
#include "nwc_drdyP.fh"
#include "msgids.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "printlevels.fh"
*::functions::
      logical  task_energy, task_gradient, task_hessian
      logical nw_inp_from_character
      external nw_inp_from_character
      external task_energy, task_gradient, task_hessian
*::passed::
      integer opcode
      character*(*)calling_routine
*::local::
      character*32 drdy_theory, bname, ename, basis
      character*255 input
      integer current_print_level
      logical status
      logical status_rtdb, status_ignore
      logical extra_sync
      integer found_task
*
*      write(fu6,*)' drdy_synch called with code ',opcode,
*     &    ' from ',calling_routine,' node =',me_drdy
*      call util_flush(6)
      if (me_drdy.eq.master) then
        extra_sync = .true.
        if (.not.
     &        rtdb_cput(my_rtdb,'scf:input vectors',1,'atomic')) then
          write(fu6,*)' could not reset input vectors to atomic '
          call errquit('drdy_synch: fatal error ',911, RTDB_ERR)
        endif
        if (.not.
     &        rtdb_cput(my_rtdb,'dft:input vectors',1,'atomic')) then
          write(fu6,*)' could not reset input vectors to atomic '
          call errquit('drdy_synch: fatal error ',911, RTDB_ERR)
        endif
*        if (.not.
*     &        rtdb_cput(my_rtdb,'mp2:input vectors',1,'atomic')) then
*          write(fu6,*)' could not reset input vectors to atomic '
*          call errquit('drdy_synch: fatal error ',911, RTDB_ERR)
*        endif
        if (.not.
     &        rtdb_cput(my_rtdb,'mcscf:input vectors',1,'atomic')) then
          write(fu6,*)' could not reset input vectors to atomic '
          call errquit('drdy_synch: fatal error ',911, RTDB_ERR)
        endif
c
c If we are not doing a single point calculation at a higher level
c of theory, just set the theory in the rtdb.
c
        if (opcode .ne. DRDY_CODE_SPENERGY) then
          if (.not.rtdb_cget(my_rtdb,'drdy:theory',1,drdy_theory)) then
            write(fu6,*)' could not determine drdy:theory'
            call errquit('drdy_synch: fatal error ',911, RTDB_ERR)
          endif
          call inp_lcase(drdy_theory)
          if (.not.rtdb_cput(my_rtdb,'task:theory',1,drdy_theory)) then
            write(fu6,*)' could not write drdy:theory -> theory'
            call errquit('drdy_synch: fatal error ',911, RTDB_ERR)
          endif
          if (rtdb_cget(my_rtdb,'drdy:input',1,input)) then
            if (.not.nw_inp_from_character(my_rtdb,input)) then
              write(fu6,*) 'could not process drdy:theory input'
              call errquit('drdy_synch: fatal error ',555, RTDB_ERR)
            endif
            extra_sync=.false.
          endif
c
c Otherwise we need to set the information for the single point
c information (theory,basis,ecp,input).
c
        else
          if (.not.rtdb_cget(my_rtdb,'drdy:sptheory',1,
     &                       drdy_theory)) then
            write(fu6,*)' could not determine drdy:sptheory'
            call errquit('drdy_synch: fatal error ',911, RTDB_ERR)
          endif
          if (.not.rtdb_cput(my_rtdb,'task:theory',1,drdy_theory)) then
            write(fu6,*)' could not write drdy:sptheory -> theory'
            call errquit('drdy_synch: fatal error ',911, RTDB_ERR)
          endif
c
c Save original information and set new info if available.  Note that
c input overwrites the old information and is essentially impossible to
c get back.  The user is responsible for input on the "theory" line of
c the drdy input to get back the "defaults".
c
          if (.not. rtdb_cget(my_rtdb,'ao basis',1,bname)) bname = ' '
          if (.not. rtdb_cget(my_rtdb,'ecp basis',1,ename)) ename = ' '
          if (rtdb_cget(my_rtdb,'drdy:spbasis',1,basis)) then
            if (.not.rtdb_cput(my_rtdb,'ao basis',1,basis)) then
              write(fu6,*) 'could not write drdy:spbasis->ao basis'
              call errquit('drdy_synch: fatal error ',555, RTDB_ERR)
            endif
          endif
          if (rtdb_cget(my_rtdb,'drdy:specp',1,basis)) then
            if (.not.rtdb_cput(my_rtdb,'ecp basis',1,basis)) then
              write(fu6,*) 'could not write drdy:specp->ecp basis'
              call errquit('drdy_synch: fatal error ',555, RTDB_ERR)
            endif
          endif
          if (rtdb_cget(my_rtdb,'drdy:spinput',1,input)) then
            if (.not.nw_inp_from_character(my_rtdb,input)) then
              write(fu6,*) 'could not process drdy:theorysp input'
              call errquit('drdy_synch: fatal error ',555, RTDB_ERR)
            endif
            extra_sync = .false.
          endif
        endif
* input routines have ga_sync() calls in them, so we have to do a ga_sync()
* call to match that, but if we don't do input reading, then we fake it
* on the master
        if (extra_sync)  then
*         So, we didn't parse any input, so fake it
          call ga_sync()
          found_task = 0
          call ga_igop(33,found_task,1,'+')
        endif
      endif

00001 continue
      if (me_drdy.ne.master) then
        call ga_sync()
        found_task = 0
        call ga_igop(33,found_task,1,'+')
      endif

*      write(fu6,*)' B4 msg_drdy_opcode = ',msg_drdy_opcode,
*     &    ' opcode = ',opcode,' node=',me_drdy
*      call util_flush(6)
      call ga_brdcst(msg_DRDY_OPCODE, opcode,
     &     MA_sizeof(MT_INT,1,MT_BYTE), master)
      call ga_sync()
*      write(fu6,*)' AF msg_drdy_opcode = ',msg_drdy_opcode,
*     &    ' opcode = ',opcode,' node=',me_drdy
*      call util_flush(6)
*     call util_print_common_data
      call util_print_get_level(current_print_level)
      call util_print_set_level(print_none)
*     call util_print_common_data
      if ((opcode .eq. DRDY_CODE_ENERGY).or.
     &    (opcode .eq. DRDY_CODE_SPENERGY)) then
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
        write(fu6,*)' unknown opcode ',opcode
        call errquit('drdy_synch: fatal error ',911, RTDB_ERR)
      endif
      call util_print_set_level(current_print_level)
      call ga_sync()
      if (me_drdy.eq.master) then
c
c If doing the single point energy, reset info to original
c
        if (opcode .eq. DRDY_CODE_SPENERGY) then
          status_ignore = rtdb_delete(my_rtdb, 'ao basis')
          status_ignore = rtdb_delete(my_rtdb, 'ecp basis')
          if (bname .ne. ' ') then
            if (.not. rtdb_cput(my_rtdb,'ao basis',1,bname)) 
     &        call errquit
     &        ('drdy_synch:failed resetting user basis name',0,
     &       RTDB_ERR)
          end if
          if (ename .ne. ' ') then
            if (.not. rtdb_cput(my_rtdb,'ecp basis',1,ename)) 
     &         call errquit
     &         ('drdy_synch:failed resetting user ecp name',0,
     &       RTDB_ERR)
          end if
        end if
        return
      endif
      opcode = 0
      goto 00001
      end
