      logical function task_hessian(rtdb)
*
* $Id: task_hessian.F,v 1.11 2000-08-09 18:48:54 bjohnson Exp $
*
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "inp.fh"
#include "util.fh"
      integer rtdb
c
c     Generic NWChem interface to compute the analytic hessian.  
c
c     RTDB input parameters
c     ---------------------
c     task:theory (string) - name of (QM) level of theory to use
c     task:numerical (logical) - optional - if true use numerical 
c         differentiation. if 
c     task:analytic  (logical) - force analytic hessian
c     
c     RTDB output parameters no for analytic hessian at the moment.
c     ----------------------
c     task:hessian file name - that has a lower triangular 
C                              (double precision) array 
c                              derivative w.r.t. geometry cart. coords.
c     task:status (logical)  - T/F for success/failure
c     task:cputime (real)    - cpu time to execute the task
c     task:walltime (real)   - wall time to execute the task
c
c     Also returns status through the function value
c
c     If the method does not have analytic derivatives automatically call
c     the numerical derivative routine (not true in my case at the time).
c
      logical stpr_gen_hess_at
      external stpr_gen_hess_at
      logical  scf,scf_hessian
      external scf,scf_hessian
c
      integer ecce_old_print, ecce_junk_print
      logical status, ignore, numerical, analytic
      double precision cpu, wall, delta_pass
      double precision default_delta
      character*80 prefix
      character*32 theory
      character*(nw_max_path_len) filehess
c
      call ecce_print_module_entry('task hessian')
c
      task_hessian = .false.
c
      cpu  = util_cpusec()
      wall = util_wallsec()
c
c     Right now only have a QM component
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task:hessian: theory not specified',0)
      if (.not. rtdb_get(rtdb,'task:numerical',mt_log,1,numerical))
     &    numerical = .true.
      if (.not. rtdb_get(rtdb,'task:analytic',mt_log,1,analytic))
     &    analytic = .false.
c
c  for numerical hessian reset task:numerical so gradients are done analytically
c
*      if (numerical) then
*        if (.not.rtdb_put(rtdb,'task:numerical',mt_log,1,.false.))
*     &      call errquit('task_hessian: failed to reset numerical',911)
*      endif

      ignore = rtdb_delete(rtdb, 'task:numerical')
      ignore = rtdb_delete(rtdb, 'task:hessian file name')
      ignore = rtdb_delete(rtdb,'task:hessian')
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_hessian: failed to invalidate status',0)

      prefix = theory ! Most common

      if ((analytic.and.inp_compare(.false.,theory,'scf')).or.
     $    (analytic.and.inp_compare(.false.,theory,'dft'))) then
        numerical = .false.
      endif
c
c     Actually do the deed
c
      if (numerical) then
        if (ga_nodeid().eq.0 .and. 
     $      util_print('task_hessian', print_low)) then
          write(6,*)
          write(6,*)
          call util_print_centered(6,
     $        'NWChem Finite-difference Hessian ',
     $        40,.true.)
          write(6,*)
          write(6,*)
        endif
c
c define default delta value
c
        default_delta = 0.001d00  ! should it be 0.0025 
        if (inp_compare(.false.,theory,'dft')) default_delta = 0.01d00
c     
c     check for change of finite difference delta
c     
        if (.not.rtdb_get(rtdb,'stpr_gen:delta',MT_DBL,1,delta_pass))
     &      delta_pass = default_delta
c
        call ecce_print_control(0, ecce_old_print) ! disable ECCE printing
        status = stpr_gen_hess_at(rtdb,delta_pass)
        call ecce_print_control(ecce_old_print,ecce_junk_print) ! re-enable ECCE printing
      else
        
        if (ga_nodeid().eq.0 .and.
     $      util_print('task_hessian', print_low)) then
          write(6,*)
          write(6,*)
          call util_print_centered(6,
     $        'NWChem Analytic Hessian ',
     $        40,.true.)
          write(6,*)
        endif
c
        if (theory .eq. 'scf' .or. theory .eq. 'dft') then
          status = scf_hessian(rtdb)
        else
          call errquit('task_hessian: unknown analytic theory',0)
        endif
c
      endif
c
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - wall
c
c     if (status) then
c        key = ' '
c        write(key,'(a,'':hessian'')') prefix(1:inp_strlen(prefix))
c     endif
c
      task_hessian = status
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, task_hessian))
     $     call errquit('task_hessian: failed to store status',0)
c
      call util_file_name('hess',  .false., .false.,filehess)
      if (.not.rtdb_cput(rtdb,'task:hessian file name',1,filehess))
     &    call errquit('task_hessian: failed to store filename',911)
      if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $    call errquit('task_gradient: failed storing cputime',0)
      if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $    call errquit('task_gradient: failed storing walltime',0)
c
      call ecce_print1('cpu time', mt_dbl, cpu, 1)
      call ecce_print1('wall time', mt_dbl, wall, 1)
      if (task_hessian) then
         call ecce_print_module_exit('task hessian', 'ok')
      else
         call ecce_print_module_exit('task hessian', 'failed')
      endif
      end

