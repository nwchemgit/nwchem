      logical function task_nadtrc(rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
      integer rtdb
c
      logical status, ignore
      logical ocons
      character*32 theory
      double precision cpu, wall
      logical task_aoints_doit
      external task_aoints_doit
      logical oecce
cc
cgk debug
*      write(*,*)'gk: entered task_nadtrc'
cgk end
      oecce = .false.
      cpu  = util_cpusec()
      wall = util_wallsec()
c
      ignore = rtdb_delete(rtdb, 'task:nadtrc')
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_nadtrc: failed to invalidate status',0,
     $       RTDB_ERR)
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task:nadtrc: theory not specified',0, 
     $       RTDB_ERR)

      if (theory .ne. 'sodft') then
        call errquit('task_nadtrc: theory not supported for nadtrc'
     &    ,0,RTDB_ERR)
      else
        status=task_aoints_doit(rtdb,theory) 
      end if
c
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - wall
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, status))
     $     call errquit('task_energy: failed to store status',0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $     call errquit('task_energy: failed storing cputime',0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $     call errquit('task_energy: failed storing walltime',0,
     &       RTDB_ERR)
c
      call ecce_print1('cpu time', mt_dbl, cpu, 1)
      call ecce_print1('wall time', mt_dbl, wall, 1)
c         
      task_nadtrc = status
cgk debug
*      write(*,*)'gk: leaving task_nadtrc'
cgk end
c
      end
