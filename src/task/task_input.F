      subroutine task_input(rtdb)
*
* $Id: task_input.F,v 1.28 2002-04-24 16:04:59 bylaska Exp $
*
      implicit none
#include "stdio.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "util.fh"
#include "mafdecls.fh"
      integer rtdb              ! [input] data base handle
c
c     Called by process 0 only.  Does nothing but read the input
c     and set corresponding the database parameters
c
c     Called whenever the input parser encounters a task directive.
c     Terminates by calling errquit if there was any error.
c
c     task [qmmm] <string theory> [<string operation = energy>] [numerical || analytic] 
c          [ignore]
c
      logical numerical, ignore, oqmmm, analytic
      character*30 test, operation
      character*32 theory
      character*128 err_string
c
      theory    = ' '
      operation = ' '
      analytic  = .false.
      numerical = .false.
      ignore    = .false.
      oqmmm     = .false.
c
      call rtdb_clean(rtdb, 'task:') ! Delete ALL previous task parameters
      
c
      if (.not. inp_a(test)) goto 1000 ! task keyword
 10   if (inp_a(test)) then
         if (inp_compare(.false.,test,'energy')) then
            operation = 'energy'
         else if (inp_compare(.false.,test,'optimize')) then
            operation = 'optimize'
         else if (inp_compare(.false.,test,'saddle')) then
            operation = 'saddle'
         else if (inp_compare(.false.,test,'lst')) then
            operation = 'lst'
         else if (inp_compare(.false.,test,'gradient')) then
            operation = 'gradient'
         else if (inp_compare(.false.,test,'hessian')) then
            operation = 'hessian'
         else if (inp_compare(.false.,test,'dynamics')) then
            operation = 'dynamics'
         else if (inp_compare(.false.,test,'thermodynamics')) then
            operation = 'thermodynamics'
         else if (inp_compare(.false.,test,'et')) then
            operation = 'et'
         else if (inp_compare(.false.,test,'freq') .or.
     $            inp_compare(.false.,test,'frequencies')) then
            operation = 'freq'
         else if (inp_compare(.false.,test,'numerical')) then
            numerical = .true.
         else if (inp_compare(.false.,test,'analytic')) then
            analytic  = .true.
         else if (inp_compare(.false.,test,'vscf')) then
            operation  = 'vscf'
         else if (inp_compare(.false.,test,'shell')) then
            theory = 'shell'
            call task_shell_input(rtdb)
         else if (inp_compare(.false.,test,'debug')) then
            if (.not. rtdb_put(rtdb, 'task:debug', MT_LOG, 1, .true.))
     $           call errquit('task_input: failed putting debug flag',0)
         else if (inp_compare(.false.,test,'ignore')) then
            ignore = .true.
         else if (inp_compare(.false.,test,'qmmm')) then
            oqmmm = .true.
         else if (inp_compare(.false.,test,
     >                        'steepest_descent')) then
            operation='steepest_descent'
         else if (inp_compare(.false.,test,
     >                        'Car-Parrinello')) then
            operation='Car-Parrinello'
         else if (inp_compare(.false.,test,
     >                        'conjugate_gradient')) then
            operation='conjugate_gradient'
         else if (inp_compare(.false.,test,
     >                        'wavefunction_initializer')) then
            operation='wavefunction_initializer'
         else if (inp_compare(.false.,test,
     >                        'v_wavefunction_initializer')) then
            operation='v_wavefunction_initializer'
         else if (inp_compare(.false.,test,
     >                        'wavefunction_expander')) then
            operation='wavefunction_expander'
         else if (inp_compare(.false.,test,
     >                        'psp_formatter')) then
            operation='psp_formatter'
         else if (inp_compare(.false.,test,
     >                        'psp_generator')) then
            operation='psp_generator'
         else if (inp_compare(.false.,test,
     >                        'cell_optimize')) then
            operation='cell_optimize'
         else if (inp_compare(.false.,test,
     >                        'stress')) then
            operation='stress'
         else if (inp_compare(.false.,test,
     >                        'pspw_dplot')) then
            operation='pspw_dplot'
         else if (theory .eq. ' ') then
            theory = test
            call inp_lcase(theory)
         else
            goto 1000
         endif
         goto 10
      endif
      call inp_clear_err()
c
      if (numerical.and.analytic) then
        write(luout,*)' ********** fatal input error ************'
        write(luout,*)' both numerical and analytic specified '
        goto 1000
      endif
      if (numerical) then
         if (.not. rtdb_put(rtdb, 'task:numerical', mt_log, 1, 
     $        numerical)) call errquit('task_input: rtdb is corrupt',0)
      endif
      if (analytic) then
         if (.not. rtdb_put(rtdb, 'task:analytic', mt_log, 1, 
     $        analytic)) call errquit('task_input: rtdb is corrupt',0)
      endif
      if (oqmmm) then
        if (.not.rtdb_put(rtdb, 'task:QMMM', mt_log, 1, .true.))
     &       call errquit('task_input:task:QMMM rtdb is corrupt',911)
        if (theory.eq.'md') call errquit
     &      ('task_input:theory must be a QM method for QMMM',911)
      endif
      if (ignore) then
         if (.not. rtdb_put(rtdb, 'task:ignore', mt_log, 1, 
     $        ignore)) call errquit('task_input: rtdb is corrupt',0)
      endif
c
      if (theory .ne. ' ') then
         if (.not. rtdb_cput(rtdb, 'task:theory', 1, theory))
     $        call errquit('task_input: rtdb is corrupt', 0)
      endif
      if (operation .ne. ' ') then
         if (.not. rtdb_cput(rtdb, 'task:operation',1,operation))
     $        call errquit('task_input: rtdb is corrupt', 0)
      endif
c
      return
c
 1000 err_string =             'usage: task <string theory>'
      err_string = err_string//' [<string operation = energy>] '
      err_string = err_string//'[numerical || analytic] [ignore]'
      call errquit(err_string,911)
c
c     THIS IS NEVER EXECUTED.  IT IS HERE BEACUSE -ltask IS AT
c     THE FRONT OF THE LINK PATH AND THE STATE ROUTINES ARE NOT
c     REFERENCED UNTIL LATER.
c
      call task_save_state(0,'a')
      call task_restore_state(0,'a')
      call task_delete_state(0,'a')
c
      end
      subroutine task_shell_input(rtdb)
      implicit none
#include "rtdb.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer rtdb
c
c     Read input of the shell task directive
c
c     task shell [(<integer-range processor = 0>|all)] <string command>
c
      character*1023 command
      character*8 test
      integer range(3)
c
      if (inp_n_field().lt.3 .or. inp_n_field().gt.4)
     $     call errquit('task shell: invalid syntax',0)
c
      if (inp_n_field() .eq. 4) then
         if (.not. inp_irange(range(1),range(2),range(3))) then
            if (inp_a(test)) then
               if (inp_compare(.false.,test,'all')) then
                  range(1) = 0
                  range(2) = ga_nnodes()-1
                  range(3) = 1
                  goto 10
               endif
            endif
         else
            goto 10
         endif
         call errquit('task shell: bad format for processor range',0)
c
 10      call inp_clear_err()
         if (range(1).lt.0 .or. range(1).ge.ga_nnodes() .or.
     $        range(2).lt.0 .or. range(2).ge.ga_nnodes() .or.
     $        range(3).lt.1)
     $        call errquit('task shell: invalid process range',0)
         if (.not. rtdb_put(rtdb, 'task:shell range', mt_int, 3, range))
     $        call errquit('task shell: rtdb corrupt',0)
      endif
c
      if (.not. inp_a(command))
     $     call errquit('task shell: failed reading command',0)
c     
      if (.not. rtdb_cput(rtdb,'task:shell command',1,command))
     $     call errquit('task shell: rtdb corrupt',0)
c
      end





