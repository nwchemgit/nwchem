C> \ingroup task
C> @{
C>
C> \brief Generic interface for Molecular Dynamics (MD) simulations
C>
C> \return Return .true. if successfull, and .false. otherwise.
C>
      logical function task_dynamics(rtdb)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
      integer rtdb !< [Input] The RTDB handle
c
c     Generic NWChem interface to perform molecular dynamics simulations
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
      logical argos,nwmd
      external argos,nwmd
c
      logical task_smd_dynamics
      external task_smd_dynamics
c
      logical  md_driver
      external md_driver
c
      logical  task_freq, leps_gradient, mc_main
      external task_freq, leps_gradient, mc_main
c
      logical status
      double precision cpu,wall
c
      character*32 theory
c
      call ecce_print_module_entry('task dynamics')
c
      task_dynamics = .false.
c
      cpu  = util_cpusec()
      wall = util_wallsec()
c
c     Right now only have a QM component
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task:energy: theory not specified',0, RTDB_ERR)
c
c     Set output parameters as if failed
c
      status = rtdb_delete(rtdb, 'task:dynamics')
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_dynamics: failed to invalidate status',0,
     &       RTDB_ERR)
c
c     Actually do the deed
c
      if(theory.eq."smd") then
c         write(*,*) "smd dynamics is disabled for now"
        status = task_smd_dynamics(rtdb)
c      else if(theory.eq.'leps') then
c        status = task_freq(rtdb)
c        status = mc_main(rtdb,leps_gradient,theory)
      elseif(theory.eq.'argos') then
        status=argos(rtdb)
      else
        status=nwmd(rtdb)
      end if
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
      task_dynamics = status
c
      end
