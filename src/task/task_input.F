      subroutine task_input(rtdb)
      implicit none
#include "rtdb.fh"
#include "inp.fh"
#include "util.fh"
      integer rtdb              ! [input] data base handle
c
c     Called by process 0 only.  Does nothing but read the input
c     and set corresponding the database parameters
c
c     Called whenever the input parser encounters a task directive.
c     Terminates by calling errquit if there was any error.
c
c     task <string theory> [<string operation = energy>]
c
      character*20 test, theory, operation
c
      logical  rjhtest, raktest, atwtest, task_check_sum, mdmain, cphf
      logical  nwargos, property, gapss
      logical  task_energy, task_gradient, task_optimize, task_freq
c
      external rjhtest, raktest, atwtest, task_check_sum, mdmain, cphf
      external nwargos, property, gapss
      external task_energy, task_gradient, task_optimize, task_freq
c      
      theory = ' '
      operation = ' '
c
      call rtdb_clean(rtdb, 'task:') ! Delete ALL previous task parameters
c
      if (.not. inp_a(test)) goto 1000 ! task keyword
 10   if (inp_a(test)) then
         if (inp_compare(.false.,test,'energy')) then
            operation = 'energy'
         else if (inp_compare(.false.,test,'optimize')) then
            operation = 'optimize'
         else if (inp_compare(.false.,test,'freq') .or.
     $            inp_compare(.false.,test,'frequencies')) then
            operation = 'freq'
         else if (theory .eq. ' ') then
            theory = test
            call inp_lcase(theory)
         else
            goto 1000
         endif
         goto 10
      endif
c
      if (theory .eq. ' ') goto 1000 ! No default theory
      if (.not. rtdb_cput(rtdb, 'task:theory', 1, theory))
     $     call errquit('task_input: rtdb is corrupt', 0)
      if (operation .ne. ' ') then
         if (.not. rtdb_cput(rtdb, 'task:operation',1,operation))
     $        call errquit('task_input: rtdb is corrupt', 0)
      endif
c
      return
c
 1000 call errquit('usage: task <string theory>'//
     $     ' [<string operation = energy>]',0)
c
c     THIS IS NEVER EXECUTED.  IT IS HERE BEACUSE -ltask IS AT
c     THE FRONT OF THE LINK PATH AND THE STATE ROUTINES ARE NOT
c     REFERENCED UNTIL LATER.
c
      call task_save_state()
      call task_restore_state()
      call task_delete_state()
c
      end

         
