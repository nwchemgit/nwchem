      logical function task_optimize(rtdb)
c     $Id: task_optimize.F,v 1.21 2000-02-29 19:18:03 windus Exp $
      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
      integer rtdb
c     
c     Optimize a geometry using stepper and the generic 
c     task energy/gradient interface.  Eventually will need another
c     layer below here to handle the selection of other optimizers.
c     
c     Since this routine can be directly invoked by application modules
c     no input is processed in this routine.
c     
c     RTDB input parameters
c     ---------------------
c     task:theory (string) - must be set for task_gradient to work
c     
c     RTDB output parameters
c     ----------------------
c     task:energy (real)   - final energy from optimization
c     task:gradient (real) - final gradient from optimization
c     task:status (real)   - T/F on success/failure
c     task:cputime
c     task:walltime
c     geometry             - final geometry from optimization
c     
      double precision cpu, wall
c     
      logical ousedriver, status, ts, omin
c     
      character*32 theory
c     
      logical  task_gradient, stpr_walk_task, driver, nwargos,nwmd
      external task_gradient, stpr_walk_task, driver, nwargos,nwmd
c     
      call ecce_print_module_entry('task optimize')
c     
      if (.not. rtdb_get(rtdb, 'opt:driver', mt_log, 1, ousedriver))
     $     ousedriver = .true.
c     
      cpu  = util_cpusec()
      wall = util_wallsec()
c     
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_optimize: failed to invalidate status',0)
c     
      if (ga_nodeid().eq.0 .and. (.not. ousedriver) .and.
     $     util_print('task_optimize', print_low)) then
         write(6,*)
         write(6,*)
         call util_print_centered(6,
     $        'NWChem Geometry Optimization',
     $        40,.true.)
         write(6,*)
         write(6,*)
      endif
c     
      if(.not.rtdb_cget(rtdb,'task:theory',1,theory))
     +     call errquit('task_optimize: failed rtdb_cget task:theory',0)
      if (theory.eq.'md') then
      status=nwmd(rtdb)
      elseif(theory.eq.'nwargos') then
      status=nwargos(rtdb)
      else if (theory.eq.'qmmm') then
         if (ga_nodeid().eq.0) then
            write(6,*)'task_optimize: qmmm should never be the theory '
            write(6,*)'task_optimize: you probably reversed the order ',
     &           'of the task input '
            write(6,*)'task [qmmm] <string theory> [<string operation ',
     &           '= energy>] [numerical] [ignore]'
         endif
         status = .false.
      else
         if (ousedriver) then
            if (.not. rtdb_put(rtdb,'driver:modsad',mt_int, 1, 0))
     $           call errquit('task_optimize: rtdb?',0)
            status = driver(rtdb)
         else
            ts   = .false.
            omin = .true.
            if (.not.rtdb_put(rtdb,'stepper:ts',MT_LOG,1,ts))
     &           call errquit('stpr_store_rtdb: rtdb_put failed',911)
            if (.not.rtdb_put(rtdb,'stepper:omin',MT_LOG,1,omin))
     &           call errquit('stpr_store_rtdb: rtdb_put failed',911)
            status = stpr_walk_task(rtdb)
*rak:            call stpr_walk_reset()
*rak: 10         if (task_gradient(rtdb)) then
*rak:               if (.not. stpr_walk(rtdb)) goto 10 ! Take next step
*rak:               status = .true.  ! Finished
*rak:            else if (ga_nodeid().eq.0) then
*rak:               write(luout,*)' task_optimize: energy+gradients failed ',
*rak:     $              ' ... optimization incomplete'
*rak:               call util_flush(6)
*rak:               status = .false.
*rak:            endif
         endif
      endif
c     
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - wall
      if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $     call errquit('task_optimize: failed storing cputime',0)
      if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $     call errquit('task_optimize: failed storing walltime',0)
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, status))
     $     call errquit('task_optimize: failed to set status',0)
c     
      call ecce_print1('cpu time', mt_dbl, cpu, 1)
      call ecce_print1('wall time', mt_dbl, wall, 1)
      if (status) then
         call ecce_print_module_exit('task optimize', 'ok')
      else
         call ecce_print_module_exit('task optimize', 'failed')
      endif
c     
      task_optimize = status
c     
      end
      logical function task_saddle(rtdb)
c     $Id: task_optimize.F,v 1.21 2000-02-29 19:18:03 windus Exp $
      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
      integer rtdb
c     
c     Evil twin of task_optimize.  Created with cut/paste, global
c     replace and tweaking of modsad, ts, and omin.
c     
c     Optimize a geometry using stepper and the generic 
c     task energy/gradient interface.  Eventually will need another
c     layer below here to handle the selection of other optimizers.
c     
c     Since this routine can be directly invoked by application modules
c     no input is processed in this routine.
c     
c     RTDB input parameters
c     ---------------------
c     task:theory (string) - must be set for task_gradient to work
c     
c     RTDB output parameters
c     ----------------------
c     task:energy (real)   - final energy from optimization
c     task:gradient (real) - final gradient from optimization
c     task:status (real)   - T/F on success/failure
c     task:cputime
c     task:walltime
c     geometry             - final geometry from optimization
c     
      double precision cpu, wall
c     
      logical ousedriver, status, ts, omin
      integer modsad
c     
      character*32 theory
c     
      logical  task_gradient, stpr_walk_task, driver, nwargos
      external task_gradient, stpr_walk_task, driver, nwargos
c     
      call ecce_print_module_entry('task saddle')
c     
      if (.not. rtdb_get(rtdb, 'opt:driver', mt_log, 1, ousedriver))
     $     ousedriver = .true.
c     
      cpu  = util_cpusec()
      wall = util_wallsec()
c     
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_saddle: failed to invalidate status',0)
c     
      if (ga_nodeid().eq.0 .and. 
     $     util_print('task_saddle', print_low) .and. 
     $     (.not. ousedriver) ) then
         write(6,*)
         write(6,*)
         call util_print_centered(6,
     $        'NWChem Transition State Search',
     $        40,.true.)
         write(6,*)
         write(6,*)
      endif
c     
      if(.not.rtdb_cget(rtdb,'task:theory',1,theory))
     +     call errquit('task_saddle: failed rtdb_cget task:theory',0)
      if (theory.eq.'md'.or.theory.eq.'nwargos') then
         status=nwargos(rtdb)
      else if (theory.eq.'qmmm') then
         if (ga_nodeid().eq.0) then
            write(6,*)'task_saddle: qmmm should never be the theory '
            write(6,*)'task_saddle: you probably reversed the order ',
     &           'of the task input '
            write(6,*)'task [qmmm] <string theory> [<string operation ',
     &           '= energy>] [numerical] [ignore]'
         endif
         status = .false.
      else
         if (ousedriver) then
            modsad = 1
            if (.not. rtdb_put(rtdb,'driver:modsad',mt_int,1,modsad))
     $           call errquit('task_saddle: rtdb corrupt?', 0)
            status = driver(rtdb)
         else
            ts   = .true.
            omin = .false.
            if (.not.rtdb_put(rtdb,'stepper:ts',MT_LOG,1,ts))
     &           call errquit('stpr_store_rtdb: rtdb_put failed',911)
            if (.not.rtdb_put(rtdb,'stepper:omin',MT_LOG,1,omin))
     &           call errquit('stpr_store_rtdb: rtdb_put failed',911)
            status = stpr_walk_task(rtdb)
*rak:            call stpr_walk_reset()
*rak: 10         if (task_gradient(rtdb)) then
*rak:               if (.not. stpr_walk(rtdb)) goto 10 ! Take next step
*rak:               status = .true.  ! Finished
*rak:            else if (ga_nodeid().eq.0) then
*rak:               write(luout,*)' task_saddle: energy+gradients failed ',
*rak:     $              ' ... optimization incomplete'
*rak:               call util_flush(6)
*rak:               status = .false.
*rak:            endif
         endif
      endif
c     
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - wall
      if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $     call errquit('task_saddle: failed storing cputime',0)
      if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $     call errquit('task_saddle: failed storing walltime',0)
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, status))
     $     call errquit('task_saddle: failed to set status',0)
c     
      call ecce_print1('cpu time', mt_dbl, cpu, 1)
      call ecce_print1('wall time', mt_dbl, wall, 1)
      if (status) then
         call ecce_print_module_exit('task saddle', 'ok')
      else
         call ecce_print_module_exit('task saddle', 'failed')
      endif
c     
      task_saddle = status
c     
      end
      logical function task_lst(rtdb)
c     $Id: task_optimize.F,v 1.21 2000-02-29 19:18:03 windus Exp $
      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
      integer rtdb
c
c     Calculate energies ( and gradients ) on a -lst- path   
c
c     Since this routine can be directly invoked by application modules
c     no input is processed in this routine.
c
c     RTDB input parameters
c     ---------------------
c     task:theory (string) - must be set for task_gradient to work
c
c     RTDB output parameters
c     ----------------------
c     task:energy (real)   - final energy from optimization
c     task:gradient (real) - final gradient from optimization
c     task:status (real)   - T/F on success/failure
c     task:cputime
c     task:walltime
c     geometry             - final geometry from optimization
c
      double precision cpu, wall
c
      logical status 
c
      character*32 theory
c
      call errquit('lst broken by rjh ',0)
      task_lst = .false.
c$$$      logical  drv_lst 
c$$$      external drv_lst
c$$$c
c$$$      call ecce_print_module_entry('task lst')
c$$$c
c$$$      cpu  = util_cpusec()
c$$$      wall = util_wallsec()
c$$$c
c$$$      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
c$$$     $     call errquit('task_lst: failed to invalidate status',0)
c$$$c
c$$$      if (ga_nodeid().eq.0 .and.
c$$$     $     util_print('task_lst', print_low)) then
c$$$         write(6,*)
c$$$         write(6,*)
c$$$         call util_print_centered(6,
c$$$     $        'NWChem LST Path Calculation',
c$$$     $        40,.true.)
c$$$         write(6,*)
c$$$         write(6,*)
c$$$      endif
c$$$c
c$$$      if(.not.rtdb_cget(rtdb,'task:theory',1,theory))
c$$$     +     call errquit('task_lst: failed rtdb_cget task:theory',0)
c$$$c
c$$$      status = drv_lst(rtdb)
c$$$c
c$$$      cpu  = util_cpusec() - cpu
c$$$      wall = util_wallsec() - wall
c$$$      if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
c$$$     $     call errquit('task_lst: failed storing cputime',0)
c$$$      if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
c$$$     $     call errquit('task_lst: failed storing walltime',0)
c$$$      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, status))
c$$$     $     call errquit('task_lst: failed to set status',0)
c$$$c
c$$$      call ecce_print1('cpu time', mt_dbl, cpu, 1)
c$$$      call ecce_print1('wall time', mt_dbl, wall, 1)
c$$$      if (status) then
c$$$         call ecce_print_module_exit('task lst', 'ok')
c$$$      else
c$$$         call ecce_print_module_exit('task lst', 'failed')
c$$$      endif
c$$$c
c$$$      task_lst = status
c$$$c
      end
