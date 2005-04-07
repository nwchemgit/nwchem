      logical function task_qmmm_optimize(rtdb)
c     $Id: task_qmmm_optimize.F,v 1.10 2005-04-07 01:36:35 marat Exp $
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
#include "qmmm.fh"
      integer rtdb
      double precision cpu, wall
c     
      logical  status
c     
      character*32 theory
      character*32 optimization
      character*32 pname
c     
      logical   driver
      external  driver

      call util_print_push() 
      call util_print_rtdb_load(rtdb, 'qmmm')
c
      pname = "task_qmmm_optimize"
c
      if(qmmm_print_debug())
     >  write(*,*) "in " //pname, ga_nodeid()
c
      if (.not. rtdb_cget(rtdb, 'qmmm:optimization', 1, optimization))
     $     call errquit('qmmm_optimize: failed ',0,0)
c
      cpu  = util_cpusec()
      wall = util_wallsec()
c     
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_optimize: failed to invalidate status',0,
     &       RTDB_ERR)
     
      if(.not.rtdb_cget(rtdb,'task:theory',1,theory))
     +     call errquit('task_optimize: failed rtdb_cget task:theory',0,
     &       RTDB_ERR)

      status = .true.

      if (optimization.eq.'all') then
         if(qmmm_print_debug())
     >      write(*,*) "going into md_em"
         call md_em()
      else if (optimization.eq.'mm') then
         if(qmmm_print_debug())
     >      write(*,*) "going into md_em"
         call md_em()
      else if (optimization.eq.'qm') then
         if (.not. rtdb_put(rtdb,'driver:modsad',mt_int, 1, 0))
     $        call errquit('task_optimize: rtdb?',0, RTDB_ERR)
         if(qmmm_print_debug())
     >      write(*,*) "going into driver"
         status = driver(rtdb)
         call mm_print_info()
      endif
c     
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - wall
      if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $     call errquit('task_optimize: failed storing cputime',0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $     call errquit('task_optimize: failed storing walltime',0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, status))
     $     call errquit('task_optimize: failed to set status',0,
     &       RTDB_ERR)
c     
      call mm_write_restart()
      task_qmmm_optimize = status
c
      call util_print_pop()
c     
      end

