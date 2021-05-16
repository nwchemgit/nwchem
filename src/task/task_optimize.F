C> \ingroup task
C> @{
C>
C> \brief The main driver for geometry optimizations 
C>
C> This routine optimizes geometries selecting one of a number of
C> optimizers. The available optimizers are: Stepper, Driver, LBFGS,
C> TROPT, and DELOC.
C>
C> \return Return .true. if an minimum was successfully found, and
C> .false. otherwise.
C>
      logical function task_optimize(rtdb)
c     $Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
      integer rtdb !< [Input] The RTDB handle
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
      logical status, ts, omin
      logical ousedriver ! use driver (default)
      logical ouselbfgs  ! use lbfgs
      logical ousesocket ! use socket
      logical ousedeloc  ! use delocalized internals
CJMC
      logical ousetropt  ! use Trust Region OPTimization
CJMC
c     
      character*32 theory
c     
      logical  task_gradient, stpr_walk_task, nwmd
      external task_gradient, stpr_walk_task, nwmd
c
c     default optimizer (bfgs) driver
      logical  driver
      external driver
c
c     lbfgs optimizer optimizer driver
      logical lbfgsb_driver
      external lbfgsb_driver
c
c     socket optimizer optimizer driver
      logical socket_driver
      external socket_driver
c
c     delocalized internals optimizer driver
      logical deloc_driver
      external deloc_driver
CJMC
c     Trust Region OPTimization
      logical tropt
      external tropt
CJMC
c     
      call ecce_print_module_entry('task optimize')
c     
c     default optimizer
      if (.not. rtdb_get(rtdb, 'opt:driver', mt_log, 1, ousedriver))
     $     ousedriver = .true.
c
c     lbfgs optimizer
      if (.not. rtdb_get(rtdb, 'driver:lbfgs', mt_log, 1, ouselbfgs))
     $     ouselbfgs = .false.

c     socket optimizer
      if (.not. rtdb_get(rtdb, 'driver:socket', mt_log, 1, ousesocket))
     $     ousesocket = .false.
c
c
c     delocalized internals optimizer
      if (.not. rtdb_get(rtdb, 'driver:deloc', mt_log, 1, ousedeloc))
     $     ousedeloc = .false.
CJMC
      if (.not. rtdb_get(rtdb, 'opt:tropt', mt_log, 1, ousetropt))
     $     ousetropt = .false.
CJMC
c
c     turn off the other options if one is set
      if(ouselbfgs) then
          ousedriver = .false.
          ousedeloc = .false.
          ousesocket = .false.
      else if (ousedeloc) then
          ousedriver = .false.
          ouselbfgs = .false.
          ousesocket = .false.
      else if (ousesocket) then
          ousedriver = .false.
          ousedeloc = .false.
          ouselbfgs = .false.
      end if
c     
      cpu  = util_cpusec()
      wall = util_wallsec()
c     
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_optimize: failed to invalidate status',0,
     &       RTDB_ERR)
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
     +     call errquit('task_optimize: failed rtdb_cget task:theory',0,
     &       RTDB_ERR)
      if (theory.eq.'md') then
      status=nwmd(rtdb)
      else if (theory.eq.'qmmm') then
         if (ga_nodeid().eq.0) then
            write(6,*)'task_optimize: qmmm should never be the theory '
            write(6,*)'task_optimize: you probably reversed the order ',
     &           'of the task input '
            write(6,*)'task [qmmm] <string theory> [<string operation ',
     &           '= energy>] [numerical] [ignore]'
         endif
         status = .false.
c      else if (theory.eq.'embed') then
c         if (ga_nodeid().eq.0) then
c           write(6,*)'task [embed] <string theory> [<string operation ',
c     &           '= optimize>] [numerical] [ignore]'
c         endif
c         status = .false.
      else
         if (ousedriver) then
            if (.not. rtdb_put(rtdb,'driver:modsad',mt_int, 1, 0))
     $           call errquit('task_optimize: rtdb?',0, RTDB_ERR)
            status = driver(rtdb)           ! default optimizer
         else if(ouselbfgs) then
            status =  lbfgsb_driver(rtdb)   ! lbfgs optimizer
         else if(ousesocket) then
            status =  socket_driver(rtdb)   ! lbfgs optimizer
         else if(ousedeloc) then
            status =  deloc_driver(rtdb)    ! delocalized internals optimizer
CJMC
         else if(ousetropt) then
            if (.not. rtdb_put(rtdb,'tropt:modsad',mt_int, 1, 0))
     $           call errquit('task_optimize: rtdb?',0, RTDB_ERR)
            status = tropt(rtdb)
CJMC
         else
            ts   = .false.
            omin = .true.
            if (.not.rtdb_put(rtdb,'stepper:ts',MT_LOG,1,ts))
     &           call errquit('stpr_store_rtdb: rtdb_put failed',911,
     &       RTDB_ERR)
            if (.not.rtdb_put(rtdb,'stepper:omin',MT_LOG,1,omin))
     &           call errquit('stpr_store_rtdb: rtdb_put failed',911,
     &       RTDB_ERR)
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
     $     call errquit('task_optimize: failed storing cputime',0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $     call errquit('task_optimize: failed storing walltime',0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, status))
     $     call errquit('task_optimize: failed to set status',0,
     &       RTDB_ERR)
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
c
C>
C> \brief The main driver for saddle point (transition state)
C> optimizations 
C>
C> This routine optimizes saddle points or transition states. In 
C> essence it is a modified version of the geometry optimizer. 
C>
C> \return Return .true. if a saddle point was successfully found, and
C> .false. otherwise.
C>
      logical function task_saddle(rtdb)
c     $Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
      integer rtdb !< [Input] The RTDB handle
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
CJMC
      logical ousetropt
CJMC
      integer modsad
c     
      character*32 theory
c     
      logical  task_gradient, stpr_walk_task, driver
      external task_gradient, stpr_walk_task, driver
CJMC
      logical  tropt
      external tropt
CJMC
c     
      call ecce_print_module_entry('task saddle')
c     
      if (.not. rtdb_get(rtdb, 'opt:driver', mt_log, 1, ousedriver))
     $     ousedriver = .true.
CJMC
      if (.not. rtdb_get(rtdb, 'opt:tropt', mt_log, 1, ousetropt))
     $     ousetropt = .false.

CJMC
c     
      cpu  = util_cpusec()
      wall = util_wallsec()
c     
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_saddle: failed to invalidate status',0,
     &       RTDB_ERR)
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
     +     call errquit('task_saddle: failed rtdb_cget task:theory',0,
     &       RTDB_ERR)
      if (theory.eq.'qmmm') then
         if (ga_nodeid().eq.0) then
            write(6,*)'task_saddle: qmmm should never be the theory '
            write(6,*)'task_saddle: you probably reversed the order ',
     &           'of the task input '
            write(6,*)'task [qmmm] <string theory> [<string operation ',
     &           '= energy>] [numerical] [ignore]'
         endif
         status = .false.
c      else if (theory.eq.'embed') then
c         if (ga_nodeid().eq.0) then
c           write(6,*)'task [embed] <string theory> [<string operation ',
c     &           '= saddle>] [numerical] [ignore]'
c         endif
c         status = .false.
      else
         if (ousedriver) then
            modsad = 1
            if (.not. rtdb_put(rtdb,'driver:modsad',mt_int,1,modsad))
     $           call errquit('task_saddle: rtdb corrupt?', 0, RTDB_ERR)
            status = driver(rtdb)
CJMC
         else if(ousetropt) then
            modsad = 1
            if (.not. rtdb_put(rtdb,'tropt:modsad',mt_int,1,modsad))
     $           call errquit('task_saddle: rtdb corrupt?', 0, RTDB_ERR)
            status = tropt(rtdb)
CJMC
         else
            ts   = .true.
            omin = .false.
            if (.not.rtdb_put(rtdb,'stepper:ts',MT_LOG,1,ts))
     &           call errquit('stpr_store_rtdb: rtdb_put failed',911,
     &       RTDB_ERR)
            if (.not.rtdb_put(rtdb,'stepper:omin',MT_LOG,1,omin))
     &           call errquit('stpr_store_rtdb: rtdb_put failed',911,
     &       RTDB_ERR)
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
     $     call errquit('task_saddle: failed storing cputime',0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $     call errquit('task_saddle: failed storing walltime',0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, status))
     $     call errquit('task_saddle: failed to set status',0,
     &       RTDB_ERR)
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
c
C> \brief The main driver to calculate energies and gradients for
C> geometries on a path (defunct)
C>
C> This routine is current defunct.
C>
      logical function task_lst(rtdb)
c     $Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
      integer rtdb !< [Input] The RTDB handle
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
c$$$      double precision cpu, wall
c
c$$$      logical status 
c
c$$$      character*32 theory
c
      call errquit('lst broken by rjh ',0, UNKNOWN_ERR)
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
c
C> \brief Calculates the minimum energy pathway between two equilibrium
C> configurations
C>
C> Calculates the minimum energy pathway connecting two equilibrium
C> geometries using the Nudged Elastic Band (NEB) method.
C>
C> \return Return .true. if a minimum pathway was found, and .false.
C> otherwise
C>
      logical function task_neb(rtdb)
      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
      integer rtdb !< [Input] The RTDB handle
c
c     Calculate the minimum energy pathway between to ionic configurations
c
c     Since this routine can be directly invoked by application modules
c     no input is processed in this routine.
c
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
      logical status

*     **** external functions ****
      logical  neb
      external neb

      status = neb(rtdb)
      task_neb = status
      return
      end
C> @}
