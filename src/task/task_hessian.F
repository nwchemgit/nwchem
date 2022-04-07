C> \ingroup task
C> @{
C>
C> \brief The main driver for analytical Hessian calculations
C>
C> This routine checks the level of theory requested and if no
C> analytical Hessians are available it will issue the evaluation of
C> numerical Hessians. The resulting Hessian is stored on a file,
C> rather than in the RTDB. The name of the Hessian file is stored
C> on the RTDB instead of the Hessian itself.
C>
C> \return Return .true. if the Hessian was calculated successfully,
C> and .false. otherwise.
C>
      logical function task_hessian(rtdb)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"
      integer rtdb !< [Input] The RTDB handle
c
c     Generic NWChem interface to compute the analytic hessian.  
c
c     RTDB input parameters
c     ---------------------
c     task:theory (string) - name of (QM) level of theory to use
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
      logical task_hessian_doit
      external task_hessian_doit
c
      logical task_bsse
      logical bsse_hessian
      external bsse_hessian
c
      double precision cpu, wall, delta_pass
      logical status
      logical lcgmin
      character*80 prefix
      character*(nw_max_path_len) filehess
c
      call ecce_print_module_entry('task hessian')
c
      task_hessian = .false.
c
      cpu  = util_cpusec()
      wall = util_wallsec()
c
      if (.not.rtdb_get(rtdb, 'dft:cgmin', mt_log, 1, lcgmin))
     &    lcgmin = .false.
c
c     Right now only have a QM component
c     but the calculation might be called with *doit as
c     the others task_*
c
      if (.not. rtdb_get(rtdb,'bsse',mt_log,1,task_bsse))
     $     task_bsse = .false.
c
      if(.not.task_bsse) then
        status= task_hessian_doit(rtdb)
        if(.not.status)
     $    call
     $  errquit('task_gradient: error task_hessian_doit',0,UNKNOWN_ERR)
      else
        status = bsse_hessian(rtdb)
        if(.not.status)
     $    call
     $  errquit('task_hessian: error call bsse hessian',911,UNKNOWN_ERR)
      endif
c
      if(.not.status)
     $    call errquit('task_hessian: error task_hessian_doit',0,0)
c
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - wall
c
      task_hessian = status
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, task_hessian))
     $     call errquit('task_hessian: failed to store status',0,
     &       RTDB_ERR)
c
      call util_file_name('hess',  .false., .false.,filehess)
      if (.not.rtdb_cput(rtdb,'task:hessian file name',1,filehess))
     &    call errquit('task_hessian: failed to store filename',911,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $    call errquit('task_gradient: failed storing cputime',0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $    call errquit('task_gradient: failed storing walltime',0,
     &       RTDB_ERR)
c
      call ecce_print1('cpu time', mt_dbl, cpu, 1)
      call ecce_print1('wall time', mt_dbl, wall, 1)
      if (task_hessian) then
         call ecce_print_module_exit('task hessian', 'ok')
      else
         call ecce_print_module_exit('task hessian', 'failed')
      endif
      end
C> @}
c
      logical function task_hessian_doit(rtdb)
c
      implicit none
#include "errquit.fh"
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
      logical  scf,hess_anal,hess_check
      external scf,hess_anal,hess_check
c
      logical task_bsse
c
      integer verb
      integer ecce_old_print, ecce_junk_print
      logical status, ignore, numerical, analytic
      double precision delta_pass
      double precision default_delta
      character*80 prefix
      character*32 theory
c
c     vdw contrib
      logical disp
      logical xc_chkdispauto
      external xc_chkdispauto
c
      call ecce_print_module_entry('task hessian')
c
c     Prevent BSSE calculations
      if (.not. rtdb_get(rtdb,'bsse',mt_log,1,task_bsse))
     $    task_bsse = .false.
c
      if (task_bsse) then
        if (.not. rtdb_put(rtdb,'bsse',mt_log,1,.false.))
     $    call errquit('task_hessian_doit:rtdb_put failed',911,RTDB_ERR)
      endif
c
c     Check for vdw (disp) calculations
      disp=.false.
      if (.not.rtdb_get(rtdb, 'dft:disp', mt_log, 1, disp))
     &      disp=.false.
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task:hessian: theory not specified',0,
     &       RTDB_ERR)
      if (.not. rtdb_get(rtdb,'task:numerical',mt_log,1,numerical))
     &    numerical = .false.
      if (.not. rtdb_get(rtdb,'task:analytic',mt_log,1,analytic))
     &    analytic = .true.
c
      ignore = rtdb_delete(rtdb, 'task:numerical')
      ignore = rtdb_delete(rtdb, 'task:hessian file name')
      ignore = rtdb_delete(rtdb,'task:hessian')
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_hessian: failed to invalidate status',0,
     &       RTDB_ERR)

      prefix = theory ! Most common
c
c  Make sure that we don't do analytic for anything other than
c  scf and dft.
c
      if ((.not.inp_compare(.false.,theory,'scf')).and.
     $    (.not.inp_compare(.false.,theory,'dft'))) then
        numerical = .true.
        analytic  = .false.
      end if
c
c     == make sure the theory supports analytic hessians ==
      if (.not.hess_check(rtdb)) then
          numerical = .true.
          analytic  = .false.
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
        if (inp_compare(.false.,theory,'dft'))
     &       default_delta = 0.01d00
        if (inp_compare(.false.,theory,'xtb')) then
          default_delta = 0.005d0   
          if (.not.rtdb_get(rtdb,'xtb:verb',mt_int,1,verb)) verb=1
          if (.not.rtdb_put(rtdb,'xtb:verb',mt_int,1,0))
     &       call errquit('task_hessian: verbosity put',0,RTDB_ERR)   
        endif
        if (inp_compare(.false.,theory,'pspw')   .or.
     &      inp_compare(.false.,theory,'band')   .or.
     &      inp_compare(.false.,theory,'paw'))
     &       default_delta = 0.01d00
c
c     check for change of finite difference delta
c
        if (.not.rtdb_get(rtdb,'stpr_gen:delta',MT_DBL,1,delta_pass))
     &      delta_pass = default_delta
c
        call ecce_print_control(0, ecce_old_print) ! disable ECCE printing
        status = stpr_gen_hess_at(rtdb,delta_pass)
        call ecce_print_control(ecce_old_print,ecce_junk_print) ! re-enable ECCE printing

        if (inp_compare(.false.,theory,'xtb')) then
          if (.not.rtdb_put(rtdb,'xtb:verb',mt_int,1,verb))
     &       call errquit('task_hessian: verbosity put',0,RTDB_ERR)   
        endif

      else  ! analytic hessian
c
c       == calculate analytic hessian ==      
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
          status = hess_anal(rtdb)
        else
          call errquit('task_hessian: unknown analytic theory',0,
     &       INPUT_ERR)
        endif
c
c       vdw contribution to analytic hessian 
c       for numerical hessians, this contribution is added via the gradients
        if(theory.eq.'dft') then
          if(disp.or.xc_chkdispauto()) call xc_vdw_to_hessian(rtdb)
        endif ! vdw
      endif  ! numerical or analytic
c
      if (task_bsse) then
        if (.not. rtdb_put(rtdb,'bsse',mt_log,1,task_bsse))
     $   call errquit('task_hessian_doit:rtdb_put failed',911,RTDB_ERR)
      endif
c
      task_hessian_doit = status
      end
