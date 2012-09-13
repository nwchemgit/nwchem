C> \ingroup task
C> @{
C>
C> \brief Generic interface to perform Reactive Flux Correlation Function
C> (RFCF) calculations
C>
C> Reactive flux correlation function calculations are also known as trajectory
C> recrossing calculations. Basically, the method relies on a dividing surface,
C> and running a molecular dynamics simulation it counts how often a trajectory
C> crosses that surface.
C>
C> \return Return .true. if successfull, and .false. otherwise.
C>
      logical function task_rfcf(rtdb)
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
c     Generic NWChem interface to perform rfcf (reactive flux correlation function) calculation
C     also known as trajectory  recrossing calulation
c     by sampling a dividing surface and running molecular dynamics simulations
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
c      logical  md_driver
c      external md_driver
      logical  mc_main
      external mc_main
      logical leps_gradient
      external leps_gradient
      logical scf_gradient
      external scf_gradient  
      logical dft_energy_gradient 
      external dft_energy_gradient
      logical task_freq
      external task_freq
      logical semi_dir_mp2,task_num_grad
      external semi_dir_mp2,task_num_grad
      logical task_gradient
      external task_gradient
c
      logical status
      double precision cpu,wall
      character*32 theory
c
      call ecce_print_module_entry('task rfcf')
c
      task_rfcf = .false.
c
      cpu  = util_cpusec()
      wall = util_wallsec()
c
c     Right now only have a QM component
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task:rfcf: theory not specified',0, RTDB_ERR)
c
c     Set output parameters as if failed
c
      status = rtdb_delete(rtdb, 'task:rfcf')
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_rfcf: failed to invalidate status',0,
     &       RTDB_ERR)
c
c     Actually do the deed
c
      if(theory.eq."dft") then
        status = task_freq(rtdb)
        status = mc_main(rtdb,task_gradient,theory)
ci        status = mc_main(rtdb,dft_energy_gradient,theory)
      else if(theory.eq.'leps') then
        status = task_freq(rtdb)
        status = mc_main(rtdb,leps_gradient,theory)
      else if (theory.eq.'scf') then
        status = task_freq(rtdb)
        status = mc_main(rtdb,scf_gradient,theory)
      else if (theory.eq.'mp2') then
        status = task_freq(rtdb)
        status = mc_main(rtdb,task_gradient,theory)
      else if (theory.eq.'direct_mp2') then
        status = task_freq(rtdb)
c        status = task_num_grad(rtdb)
        status = mc_main(rtdb,task_num_grad,theory)
      else
      call errquit('task_rfcf: invalid theory: only dft,scf and leps',0,
     &       RTDB_ERR)

      end if
c
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - wall
c
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, status))
     $   call errquit('task_rfcf: failed to set status',0, RTDB_ERR)
c
      if (status) then
c
         if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $        call errquit('task_rfcf: failed storing cputime',0,
     &       RTDB_ERR)
         if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $        call errquit('task_rfcf: failed storing walltime',0,
     &       RTDB_ERR)
c
      endif
c
      call ecce_print1('cpu time', mt_dbl, cpu, 1)
      call ecce_print1('wall time', mt_dbl, wall, 1)
      if (status) then
         call ecce_print_module_exit('task rfcf', 'ok')
      else
         call ecce_print_module_exit('task rfcf', 'failed')
      endif
c         
      task_rfcf = status
c
      end
C> @}

