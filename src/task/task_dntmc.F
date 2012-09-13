C> \ingroup task
C> @{
C>
C> \brief The driver for Dynamic Nucleation Theory Monte Carlo (DNTMC)
C>
C> \return Return .true. if the calculation was successfull and .false. 
C> otherwise.
C>
      logical function task_dntmc(rtdb)
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
c     Generic NWChem interface to perform dynamic nuclear theory
c     Monte Carlo
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
      logical dntmc_gibbs
      external dntmc_gibbs
c
      logical status
      double precision cpu,wall
      character*32 theory
c
c      call ecce_print_module_entry('task dynamics')
c
      task_dntmc = .false.
c
      cpu  = util_cpusec()
      wall = util_wallsec()
c
c     Right now only have a Shawn's potentials.
c     This will be important later.
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task:dntmc: theory not specified',0, RTDB_ERR)
c
c     Set output parameters as if failed
c
      status = rtdb_delete(rtdb, 'task:dntmc')
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_dntmc: failed to invalidate status',0,
     &       RTDB_ERR)
c
c     Actually do the deed
c
      status=dntmc_gibbs(rtdb)
c
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - wall
c
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, status))
     $   call errquit('task_dntmc: failed to set status',0, RTDB_ERR)
c
      if (status) then
c
         if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $        call errquit('task_dntmc: failed storing cputime',0,
     &       RTDB_ERR)
         if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $        call errquit('task_dntmc: failed storing walltime',0,
     &       RTDB_ERR)
c
      endif
c
c     call ecce_print1('cpu time', mt_dbl, cpu, 1)
c     call ecce_print1('wall time', mt_dbl, wall, 1)
c     if (status) then
c        call ecce_print_module_exit('task dynamics', 'ok')
c     else
c        call ecce_print_module_exit('task dynamics', 'failed')
c     endif
c         
      task_dntmc = status
c
      end
C> @}

