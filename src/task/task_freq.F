      logical function task_freq(rtdb)
c     
c     Central difference calculation of the hessian using
c     the generic energy/gradient interface.  Uses a routine inside
c     stepper to do the finite difference ... this needs to be
c     cleaned up to be independent of stepper.
c
c     Also will be hooked up to analytic methods as they are available.
c
c     Since this routine is directly invoked by application modules
c     no input is processed in this routine.
c
c     RTDB input parameters
c     ---------------------
c     task:theory
c
c     RTDB output parameters
c     ----------------------
c     task:hessian file name (string) - name of file containing hessian
c     task:status (logical)      - T/F on success/failure
c     task:cputime
c     task:walltime
c
      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "stdio.fh"
#include "global.fh"
#include "util.fh"
c     
      integer rtdb
c     
      double precision delta_pass
      double precision dbl_tmp
      integer nat, geom
      logical status
      character*(nw_max_path_len) filehess
c     
      double precision cpu, wall
      logical stpr_gen_hess_at
      external stpr_gen_hess_atc
c     
      cpu  = util_cpusec()
      wall = util_wallsec()
c
      status = rtdb_delete(rtdb, 'task:hessian file name')
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_freq: failed to invalidate status',0)
      call util_file_name('hess',  .false., .false.,filehess)
c     
c     check for change of finite difference delta
c     
      if (rtdb_get(rtdb,'stpr_gen:delta',MT_DBL,1,dbl_tmp)) then
         delta_pass = dbl_tmp
      else
         delta_pass = 0.01d00
      endif
c     
      if (ga_nodeid().eq.0 .and. 
     $     util_print('task_freq', print_low)) then
         write(6,*)
         write(6,*)
         call util_print_centered(6,
     $        'NWChem Finite-difference Hessian and Frequency Analysis',
     $        40,.true.)
         write(6,*)
         write(6,*)
      endif
c
      if (.not.stpr_gen_hess_at(rtdb,delta_pass))
     $     call errquit
     $     ('task_freq: stpr_gen_hess_at failed ',911)
c     
      if ((ga_nodeid()).eq.0) then
         status = rtdb_parallel(.false.)
c
c     create/load reference geometry
c
         if (.not.geom_create(geom,'geometry')) call errquit
     $        ('task_freq:geom_create failed?',1)
         if (.not.geom_rtdb_load(rtdb,geom,'geometry'))
     $        call errquit
     $        ('task_freq:geom_rtdb_load failed?',2)
         if (.not. geom_ncent(geom,nat)) call errquit
     $        ('task_freq:geom_ncent failed?',3)
         if (.not. geom_destroy(geom)) call errquit
     $        ('task_freq:geom_destroy failed?',911)
         call vib_vib(rtdb,filehess,.true.,
     $        0,.false.,0,.false.,nat)
         status = rtdb_parallel(.true.)
      endif
      call ga_sync()
c
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - cpu
c
      if (.not. rtdb_cput(rtdb, 'task:hessian file name', 1,
     $     filehess)) call errquit
     $     ('task_freq: failed storing hessian file name',0)
      if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $     call errquit('task_freq: failed storing cputime',0)
      if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $     call errquit('task_freq: failed storing walltime',0)
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .true.))
     $     call errquit('task_freq: failed to set status',0)
c
      task_freq = .true.
c
      end
