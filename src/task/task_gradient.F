      logical function task_gradient(rtdb)
*
* $Id: task_gradient.F,v 1.29 2001-01-10 21:53:50 windus Exp $
*
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"
#include "global.fh"
#include "nwc_const.fh"
c
      integer rtdb
c
      logical status, ignore, numerical
      character*32 theory
      double precision cpu, wall, energy, dipole(3), 
     $     gradient(3,nw_max_atom), gmax, gnorm
      integer i, k, ncent, nelem
      logical task_gradient_doit, oniom_gradient, geom_rtdb_ncent
      external task_gradient_doit, oniom_gradient, geom_rtdb_ncent
c
      cpu  = util_cpusec()
      wall = util_wallsec()
c
c     Set output parameters as if failed
c
      ignore = rtdb_delete(rtdb, 'task:energy')
      ignore = rtdb_delete(rtdb, 'task:gradient')
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_gradient: failed to invalidate status',0)
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task:gradient: theory not specified',0)
c
      if (.not. geom_rtdb_ncent(rtdb, 'geometry', ncent))
     $     call errquit('task_gradient: failed getting #atoms',0)
      nelem = ncent*3
c
      if (.not. rtdb_get(rtdb, 'task:numerical', mt_log, 1, 
     $     numerical)) numerical = .false.
c
      if (theory .ne. 'oniom' .or. numerical) then
         status = task_gradient_doit(rtdb,theory,energy, gradient)
      else
         status = oniom_gradient(rtdb)
         if (status) then
            if (.not. rtdb_get(rtdb, 'oniom:energy', mt_dbl, 1, energy))
     $           call errquit('task_energy: no energy in rtdb?',0)
            if (.not. rtdb_get(rtdb, 'oniom:gradient', mt_dbl, nelem, 
     $           gradient))
     $           call errquit('task_gradient: no gradient in rtdb?',0)
            if (rtdb_get(rtdb, 'oniom:dipole', mt_dbl, 3, dipole)) then
               if (.not. rtdb_put(rtdb,'task:dipole',mt_dbl,3,dipole))
     $              call errquit('task_energy: failed storing dipole',0)
            else
               ignore = rtdb_delete(rtdb, 'task:dipole')
            endif
         endif
      endif
c
      if (status) then
         call geom_print_rtdb_ecce(rtdb)
         call ecce_print1('total energy', mt_dbl, energy, 1)
         call ecce_print2('total gradient', mt_dbl, gradient,
     $        3,3,nelem/3)
         gmax = 0.0d0
         gnorm = 0.0d0
         do i = 1, nelem/3
            do k = 1, 3
               gnorm = gnorm + gradient(k,i)**2
               gmax  = max(abs(gradient(k,i)),gmax)
            end do
         enddo
         gnorm = sqrt(gnorm)
         call ecce_print1('gradient norm', mt_dbl, gnorm, 1)
         call ecce_print1('gradient max', mt_dbl, gmax, 1)
      endif
c
      if (.not. rtdb_put(rtdb, 'task:energy', mt_dbl, 1, energy))
     $     call errquit('task_gradient: failed storing energy',0)
      if (.not. rtdb_put(rtdb, 'task:gradient', mt_dbl, nelem,
     $     gradient))
     $     call errquit('task_gradient: failed storing energy',0)
c
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - wall
c
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .true.))
     $     call errquit('task_gradient: failed to set status',0)
      if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $     call errquit('task_gradient: failed storing cputime',0)
      if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $     call errquit('task_gradient: failed storing walltime',0)
c
      call ecce_print1('cpu time', mt_dbl, cpu, 1)
      call ecce_print1('wall time', mt_dbl, wall, 1)
      if (status) then
         call ecce_print_module_exit('task gradient', 'ok')
      else
         call ecce_print_module_exit('task gradient', 'failed')
      endif
c
      task_gradient = status
c
      end
      logical function task_gradient_doit(rtdb,theory,energy,gradient)
*
* $Id: task_gradient.F,v 1.29 2001-01-10 21:53:50 windus Exp $
*
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"
#include "global.fh"
      integer rtdb
      character*32 theory
      double precision energy, gradient(3,*)
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
c     task:status (logical)- T/F for success/failure
c     if (status) then
c     .  task:energy (real)   - total energy
c     .  task:gradient (real array) - derivative w.r.t. geometry cart. coords.
c     .  task:dipole (real(3)) - total dipole if available
c     .  task:cputime (real)  - cpu time to execute the task
c     .  task:walltime (real) - wall time to execute the task
c
c     Also returns status through the function value
c
c     If the method does not have analytic derivatives automatically call
c     the numerical derivative routine.
c
      logical  scf, rimp2, oimp2, direct_mp2, ccsd, mcscf, selci
      external scf, rimp2, oimp2, direct_mp2, ccsd, mcscf, selci
      logical semi_dir_mp2, rimp2g, pspw_gradient
      external semi_dir_mp2, rimp2g, pspw_gradient
      logical dft_energy_gradient, sodft_energy_gradient
      external dft_energy_gradient, sodft_energy_gradient
c
      logical mcscf_gradient, scf_gradient, grad_dft, task_num_grad
      external mcscf_gradient, scf_gradient, grad_dft, task_num_grad
      logical geom_rtdb_ncent
c
      integer nelem, ncent
      logical status, numerical, ocosmo
      double precision dipole(3)
      character*80 key, prefix
c
      if (.not. geom_rtdb_ncent(rtdb, 'geometry', ncent))
     $     call errquit('task_gradient: failed getting #atoms',0)
      nelem = ncent*3
c
      call ecce_print_module_entry('task gradient')
c
      task_gradient_doit = .false.
c
c     Set the theory in the database so that all is consistent for
c     the lower level modules.  Should not really have to do this,
c     but this value is used (appropriately and inappropriately) and
c     things are not yet disentangled.
c
      if (.not. rtdb_cput(rtdb, 'task:theory', 1, theory))
     $     call errquit('task_gradient: setting theory?',0)
c
c     Right now only have a QM component
c
      if (.not. rtdb_get(rtdb, 'task:numerical', mt_log, 1, 
     $     numerical)) numerical = .false.
c
c Currently cosmo does not work with analytic gradients
c
      if (.not. rtdb_get(rtdb,'slv:cosmo',mt_log,1,ocosmo))
     $     ocosmo = .false.
      if (ocosmo) numerical = .true.
c
c     Actually do the deed
c
      prefix = theory           ! Most common
c
      if (numerical) then
         status = task_num_grad(rtdb)
      else if (theory .eq. 'scf') then
         status = scf_gradient(rtdb)
      else if (theory .eq. 'dft') then
         status = dft_energy_gradient(rtdb)
      else if (theory .eq. 'sodft') then
         status = sodft_energy_gradient(rtdb)
      else if (theory .eq. 'mp2' .or. theory .eq. 'semi_dir_mp2') then
         prefix = 'mp2'
         if(.not. rtdb_cput(rtdb,'sdmp2:task',1,'gradient'))
     $        call errquit('task_gradient: rtdb_put failed',0)
         status = semi_dir_mp2(rtdb)
      else if (theory .eq. 'mcscf') then
         status = mcscf_gradient(rtdb)
      else if (theory .eq. 'rimp2') then
         status = task_num_grad(rtdb) ! Only numerical
      else if (theory .eq. 'rimp2_grad') then
         status = task_num_grad(rtdb) ! Only numerical
      else if (theory .eq. 'direct_mp2') then
         status = task_num_grad(rtdb) ! Only numerical
      else if (theory .eq. 'selci') then
         status = task_num_grad(rtdb) ! Only numerical
      else if (theory .eq. 'oimp2') then
         status = .false.       ! No derivatives yet
      else if (theory .eq. 'ccsd') then
         status = task_num_grad(rtdb) ! Only numerical
      else if (theory .eq. 'ccsd(t)') then
         prefix = 'ccsd'
         status = task_num_grad(rtdb) ! Only numerical
      else if (theory .eq. 'ccsd+t(ccsd)') then
         prefix = 'ccsd'
         status = task_num_grad(rtdb) ! Only numerical
      else if (theory .eq. 'cg_pspw') then
         status = pspw_gradient(rtdb)
      else if (theory .eq. 'qmmm') then
         if (ga_nodeid().eq.0) then
           write(6,*)'task_gradient: qmmm should never be the theory '
           write(6,*)'task_gradient: you probably reversed the order ',
     &         'of the task input '
           write(6,*)'task [qmmm] <string theory> [<string operation ',
     &         '= energy>] [numerical] [ignore]'
         endif
         status = .false.
      else
         call errquit('task_gradient: unknown theory',0)
      endif
c
      if (status) then
         key = ' '
         write(key,'(a,'':energy'')') prefix(1:inp_strlen(prefix))
         if (.not. rtdb_get(rtdb, key, mt_dbl, 1, energy))
     $        call errquit('task_gradient: no energy in rtdb?',0)
c
         key = ' '
         write(key,'(a,'':dipole'')') prefix(1:inp_strlen(prefix))
         if (rtdb_get(rtdb, key, mt_dbl, 3, dipole)) then
           if (.not. rtdb_put(rtdb, 'task:dipole', mt_dbl, 3, dipole))
     $        call errquit('task_gradient: failed storing dipole',0)
	 endif
c
         write(key,'(a,'':gradient'')') prefix(1:inp_strlen(prefix))
         if (.not. rtdb_get(rtdb, key, mt_dbl, nelem, gradient))
     $        call errquit('task_gradient: no gradient in rtdb?',0)
      end if
c
      task_gradient_doit = status
c
      end

