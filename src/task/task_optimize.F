      logical function task_optimize(rtdb)
c     $Id: task_optimize.F,v 1.4 1996-11-02 07:05:14 d3g681 Exp $
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
      logical ousedriver, status
c
      character*20 theory
c
      logical  task_gradient, stpr_walk, drv_opt, nwargos
      external task_gradient, stpr_walk, drv_opt, nwargos
c
      call ecce_print_module_entry('task optimize')
c
      if (.not. rtdb_get(rtdb, 'opt:driver', mt_log, 1, ousedriver))
     $     ousedriver = .false.
c     
      cpu  = util_cpusec()
      wall = util_wallsec()
c
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_optimize: failed to invalidate status',0)
c
      if (ga_nodeid().eq.0 .and. 
     $     util_print('task_optimize', print_low)) then
         write(6,*)
         write(6,*)
         call util_print_centered(6,'NWChem Geometry Optimization',
     $        40,.true.)
         write(6,*)
         write(6,*)
      endif
c
      if(.not.rtdb_cget(rtdb,'task:theory',1,theory))
     + call errquit('task_optimize: failed rtdb_cget task:theory',0)
      if(theory.eq.'md'.or.theory.eq.'nwargos') then
      status=nwargos(rtdb)
      else
      if (ousedriver) then
         status = drv_opt(rtdb)
      else
 10      if (task_gradient(rtdb)) then
            if (.not. stpr_walk(rtdb)) goto 10 ! Take next step
            status = .true. ! Finished
         else if (ga_nodeid().eq.0) then
            write(luout,*)' task_optimize: energy+gradients failed ',
     $           ' ... optimization incomplete'
            call util_flush(6)
            status = .false.
         endif
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
      if (status) then
         call ecce_print_module_exit('task optimize', 'ok')
      else
         call ecce_print_module_exit('task optimize', 'failed')
      endif
c
      task_optimize = status
c
      end
