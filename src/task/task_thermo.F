      logical function task_thermo(rtdb)
*
* $Id: task_thermo.F,v 1.3 2000-10-23 16:24:06 windus Exp $
*
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
      integer rtdb
c
c     Generic NWChem interface to perform free energy simulations
c
c     RTDB input parameters
c     ---------------------
c     task:theory (string) - name of level of theory to use
c     
c     RTDB output parameters
c     ----------------------
c     task:status (logical)- T/F for success/failure
c     if (status) then
c     .  task:cputime (real)  - cpu time to execute the task
c     .  task:walltime (real) - wall time to execute the task
c
c     Also returns status through the function value
c
      logical nwargos,nwmd
      external nwargos,nwmd
c
      logical status
      double precision cpu,wall
      character*80 prefix
      character*32 theory
c
      call ecce_print_module_entry('task thermodynamics')
c
      task_thermo = .false.
c
      cpu  = util_cpusec()
      wall = util_wallsec()
c
c     Right now only have a QM component
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task:energy: theory not specified',0)
c
c     Set output parameters as if failed
c
      status = rtdb_delete(rtdb, 'task:thermodynamics')
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_thermo: failed to invalidate status',0)
c
c     Actually do the deed
c
      prefix=theory           ! Most common scenario
      if(theory.eq.'md') then
      status=nwmd(rtdb)
      elseif(theory.eq.'nwargos') then
      status=nwargos(rtdb)
      else
      call errquit('task_thermo: unknown theory',0)
      endif
c
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - wall
c
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, status))
     $   call errquit('task_thermo: failed to set status',0)
c
      if (status) then
c
         if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $        call errquit('task_dynamics: failed storing cputime',0)
         if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $        call errquit('task_dynamics: failed storing walltime',0)
c
      endif
c
      call ecce_print1('cpu time', mt_dbl, cpu, 1)
      call ecce_print1('wall time', mt_dbl, wall, 1)
      if (status) then
         call ecce_print_module_exit('task thermodynamics', 'ok')
      else
         call ecce_print_module_exit('task thermodynamics', 'failed')
      endif
c         
      task_thermo = status
c
      end

