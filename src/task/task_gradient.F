      logical function task_gradient(rtdb)
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
      integer rtdb
c
c     Generic NWChem interface to compute the energy and gradient.  
c     Currently only the QM components are supported.
c
c     Since this routine is directly invoked by application modules
c     no input is processed in this routine.
c
c     RTDB input parameters
c     ---------------------
c     task:theory (string) - name of (QM) level of theory to use
c     task:numerical (logical) - optional - if true use numerical 
c         differentiation. If absent or false use default selection.
c     
c     RTDB output parameters
c     ----------------------
c     task:energy (real)   - total energy
c     task:gradient (real array) - derivative w.r.t. geometry cart. coords.
c     task:status (logical)- T/F for success/failure
c     task:cputime (real)  - cpu time to execute the task
c     task:walltime (real) - wall time to execute the task
c
c     Also returns status through the function value
c
c     If the method does not have analytic derivatives automatically call
c     the numerical derivative routine.
c
      logical  scf, nwdft, rimp2, oimp2, direct_mp2, ccsd, mcscf, selci
      external scf, nwdft, rimp2, oimp2, direct_mp2, ccsd, mcscf, selci
      logical semi_dir_mp2
      external semi_dir_mp2
c
      logical scf_gradient, grad_dft, task_num_grad
      external scf_gradient, grad_dft, task_num_grad
c
      integer k_grad, l_grad, type, nelem
      logical status, numerical
      double precision energy, cpu, wall
      character*80 theory, key
c
      task_gradient = .false.
c
      cpu  = util_cpusec()
      wall = util_wallsec()
c
c     Right now only have a QM component
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task:gradient: theory not specified',0)
      if (.not. rtdb_get(rtdb, 'task:numerical', mt_log, 1, 
     $     numerical)) numerical = .false.
c
c     Set output parameters as if failed
c
      status = rtdb_delete(rtdb, 'task:energy')
      status = rtdb_delete(rtdb, 'task:gradient')
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_gradient: failed to invalidate status',0)
c
c     Actually do the deed
c
      if (numerical) then
         status = task_num_grad(rtdb)
      else if (theory .eq. 'scf') then
         if (.not.rtdb_put(rtdb, 'scf:lagrangian',mt_log, 1, .true.))
     $        call errquit('task_gradient: rtdb is corrupt?',0)
         status = scf_gradient(rtdb)
      else if (theory .eq. 'dft') then
         status = nwdft(rtdb)
         if (status) status = grad_dft(rtdb)
      else if (theory .eq. 'semi_dir_mp2') then
         status = semi_dir_mp2(rtdb)
      else if (theory .eq. 'rimp2') then
         status = .false.       ! No numerical derivatives yet
      else if (theory .eq. 'oimp2') then
         status = .false.       ! No numerical derivatives yet
      else if (theory .eq. 'direct_mp2') then
         status = .false.       ! No numerical derivatives yet
      else if (theory .eq. 'ccsd') then
         status = .false.       ! No numerical derivatives yet
      else if (theory .eq. 'mcscf') then
         status = .false.       ! No numerical derivatives yet
      else if (theory .eq. 'selci') then
         status = .false.       ! No numerical derivatives yet
      else
         call errquit('task_gradient: unknown theory',0)
      endif
c
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - wall
c
      if (status) then
         key = ' '
         write(key,'(a,'':energy'')') theory(1:inp_strlen(theory))
         if (.not. rtdb_get(rtdb, key, mt_dbl, 1, energy))
     $        call errquit('task_gradient: no energy in rtdb?',0)
         write(key,'(a,'':gradient'')') theory(1:inp_strlen(theory))
         if (.not. rtdb_ma_get(rtdb, key, type, nelem, l_grad))
     $        call errquit('task_gradient: no gradient in rtdb?',0)
         if (.not. ma_get_index(l_grad, k_grad))
     $        call errquit('task_gradient: ma corrupt?',0)
c         
         if (.not. rtdb_put(rtdb, 'task:energy', mt_dbl, 1, energy))
     $        call errquit('task_gradient: failed storing energy',0)
         if (.not. rtdb_put(rtdb, 'task:gradient', mt_dbl, nelem,
     $        dbl_mb(k_grad)))
     $        call errquit('task_gradient: failed storing energy',0)
         if (.not. ma_free_heap(l_grad))
     $        call errquit('task_gradient: ma corrupt?',0)
      endif
      if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $     call errquit('task_gradient: failed storing cputime',0)
      if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $     call errquit('task_gradient: failed storing walltime',0)
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .true.))
     $     call errquit('task_gradient: failed to set status',0)
c
      task_gradient = status
c
      end

