      subroutine task(rtdb)
      implicit none
#include "rtdb.fh"
#include "inp.fh"
#include "util.fh"
      integer rtdb              ! [input] data base handle
c
c     Called by ALL processes.  task_input has read the
c     task directive and put stuff into the database.  Get the
c     data out and invoke the desired action.
c
c     If the operation is in the list of those supported by generic 
c     routines (currently energy, gradient, hessian, optimize), then 
c     the generic routine is called.  Otherwise, a match is attemped 
c     for a specialized routine.  If no operation is specified 
c     and no specialized routine located, then it is assumed that
c     a generic energy calculation is required.
c
c     This needs extending to accomodate QM/MM and other mixed methods
c     by having both MM and QM pieces specified (e.g., task md dft).
c
      character*20 theory, operation
      logical status
      character*80 errmsg
c
      logical  rjhtest, raktest, atwtest, task_check_sum, mdmain, cphf
      logical  nwargos, property, gapss
      logical  task_energy, task_gradient, task_optimize, task_freq
c
      external rjhtest, raktest, atwtest, task_check_sum, mdmain, cphf
      external nwargos, property, gapss
      external task_energy, task_gradient, task_optimize, task_freq
c      
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task: no task input for theory?', 0)
      if (.not. rtdb_cget(rtdb, 'task:operation', 1, operation))
     $     operation = ' '
c
      if (operation .eq. ' ') then ! Look for specialized routine
         if (inp_compare(.false.,'rjhtest', theory)) then
            status = rjhtest(rtdb)
         else if (inp_compare(.false.,'raktest', theory)) then
            status = raktest(rtdb)
         else if (inp_compare(.false.,'atwtest', theory)) then
            status = atwtest(rtdb)
         else if (inp_compare(.false.,'rtdbprint', theory)) then
            status = rtdb_print(rtdb,.true.)
         else if (inp_compare(.false.,'checksum', theory)) then
            status = task_check_sum(rtdb)
         else if (inp_compare(.false.,'md', theory)) then
            status = mdmain(rtdb)
         else if (inp_compare(.false.,'cphf', theory)) then
            status = cphf(rtdb)
         else if (inp_compare(.false.,'nwargos', theory)) then
            status = nwargos(rtdb)
         else if (inp_compare(.false.,'property', theory)) then
            status = property(rtdb)
         else if (inp_compare(.false.,'gapss', theory)) then
            status = gapss(rtdb) ! This could/should be in task_energy
         else
            operation = 'energy' ! No matching specialized operation
            goto 100            !  execute the default generic operation
         endif
         if (.not. status) then
            errmsg = ' '
            write(errmsg,'(a,'' failed '')')theory(1:inp_strlen(theory))
            call errquit(errmsg, 0)
         endif
         return
      endif
c
 100  if (inp_compare(.false., 'energy', operation)) then
         status = task_energy(rtdb)
      else if (inp_compare(.false., 'gradient', operation)) then
         status = task_gradient(rtdb)
      else if (inp_compare(.false., 'optimize', operation)) then
         status = task_optimize(rtdb)
      else if (inp_compare(.false., 'freq', operation)) then
         status = task_freq(rtdb)
      else
         call errquit('task: unknown generic operation',0)
      endif
c
      end

         
