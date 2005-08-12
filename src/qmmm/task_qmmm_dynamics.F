      logical function task_qmmm_dynamics(rtdb)
*
* $Id: task_qmmm_dynamics.F,v 1.3 2005-08-12 01:12:11 marat Exp $
*
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
      integer rtdb
c
      logical status
      double precision cpu,wall
c
c      call util_print_push() 
c      call util_print_rtdb_load(rtdb, 'qmmm')
c
      call ecce_print_module_entry('task dynamics')
c
      task_qmmm_dynamics = .false.
c
      cpu  = util_cpusec()
      wall = util_wallsec()

c     Set output parameters as if failed
c
      status = rtdb_delete(rtdb, 'task:dynamics')
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_dynamics: failed to invalidate status',0,
     &       RTDB_ERR)
c
c     Actually do the deed
c
      call md_md()
      status = .true.
c
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - wall
c
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, status))
     $   call errquit('task_dynamics: failed to set status',0, RTDB_ERR)
c
      if (status) then
c
         if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $        call errquit('task_dynamics: failed storing cputime',0,
     &       RTDB_ERR)
         if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $        call errquit('task_dynamics: failed storing walltime',0,
     &       RTDB_ERR)
c
      endif
c
      call ecce_print1('cpu time', mt_dbl, cpu, 1)
      call ecce_print1('wall time', mt_dbl, wall, 1)
      if (status) then
         call ecce_print_module_exit('task dynamics', 'ok')
      else
         call ecce_print_module_exit('task dynamics', 'failed')
      endif
c         
      task_qmmm_dynamics = status
c
c      call util_print_pop()
c
      end

