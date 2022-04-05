C> \ingroup task
C> @{
C>
C> \brief The main driver for gradient evaluations
C>
C> This routine controls both analytical and numerical gradient
C> evaluations. Before starting any gradient evaluation it will check
C> whether the calculation was preceeded by a successfull energy
C> evaluation. If not the energy will be calculated first.
C>
C> Dependent on the theory setting a variety of gradient drivers
C> may be called. If no analytical gradients are available for a theory
C> the gradient will be evaluated using numerical differentiation.
C>
C> \return Return .true. if the gradient was calculated successfully,
C> and .false. otherwise.
C>
      logical function task_gradient(rtdb)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"
#include "global.fh"
#include "nwc_const.fh"
c
      integer rtdb !< [Input] The RTDB handle
c
      logical status, ignore, numerical
      logical ocons
      logical omm
      logical task_qmmm
      character*32 theory
      double precision cpu, wall, energy, dipole(3), 
     $     gradient(3,nw_max_atom), gmax, gnorm
      integer i, k, ncent, nelem
      logical task_gradient_doit, oniom_gradient, geom_rtdb_ncent
      external task_gradient_doit, oniom_gradient, geom_rtdb_ncent
      logical  task_qmmm_gradient
      external task_qmmm_gradient
      logical cons_active
      external cons_active
c
      logical task_bsse
      logical bsse_gradient
      external bsse_gradient
c     TP: added variable
      character*30 operation
c
c      logical task_embed
c      logical embed_gradient
c      external embed_gradient
c
      cpu  = util_cpusec()
      wall = util_wallsec()
c
c     Set output parameters as if failed
c
      ignore = rtdb_delete(rtdb, 'task:energy')
      ignore = rtdb_delete(rtdb, 'task:gradient')
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_gradient: failed to invalidate status',0,
     &       RTDB_ERR)
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task:gradient: theory not specified',0,
     &       RTDB_ERR)
c
      if (.not. geom_rtdb_ncent(rtdb, 'geometry', ncent))
     $     call errquit('task_gradient: failed getting #atoms',0,
     &       RTDB_ERR)
      nelem = ncent*3
c
      if (.not. rtdb_get(rtdb, 'task:numerical', mt_log, 1, 
     $     numerical)) numerical = .false.
c
      if( .not. rtdb_get(rtdb,'task:mm',mt_log,1,omm))
     &    omm = .false.
c
      if( .not. rtdb_get(rtdb,'task:QMMM',mt_log,1,task_qmmm))
     &    task_qmmm = .false.
c
c      if( .not. rtdb_get(rtdb,'task:embed',mt_log,1,task_embed))
c     &    task_embed = .false.
c
      if (.not. rtdb_get(rtdb,'bsse',mt_log,1,task_bsse))
     $     task_bsse = .false.
c
      if(task_qmmm) then
          status = task_qmmm_gradient(rtdb)
            if (.not. rtdb_get(rtdb, 'qmmm:energy', mt_dbl, 1, energy))
     $           call errquit('task_energy: no energy in rtdb?',0,
     &       RTDB_ERR)
            if (.not. rtdb_get(rtdb, 'qmmm:gradient', mt_dbl, nelem, 
     $           gradient))
     $           call errquit('task_gradient: no gradient in rtdb?',0,
     &       RTDB_ERR)
c
c      else if(task_embed) then
c          status = embed_gradient(rtdb)
c            if (.not. rtdb_get(rtdb, 'embed:energy', mt_dbl, 1, energy))
c     $           call errquit('task_gradient: no energy in rtdb?',0,
c     &       RTDB_ERR)
c            if (.not. rtdb_get(rtdb, 'embed:gradient', mt_dbl, nelem,
c     $           gradient))
c     $           call errquit('task_gradient: no gradient in rtdb?',0,
c     &       RTDB_ERR)
c
      else if(theory .ne. 'oniom' .or. numerical) then
c
        if(.not.task_bsse) then
c       TP: added the case where mm is present in task directive
          if(omm) then
            call mm_task_gradient(rtdb,theory,energy,gradient,status)
          else 
            status = task_gradient_doit(rtdb,theory,energy, gradient)
          end if ! if(omm)
        else
          status = bsse_gradient(rtdb,theory,energy,gradient)
          if(.not.status)
     $      call errquit('task_gradient: error call bsse gradient',911,
     &                   CALC_ERR)
        endif
c
      else
         status = oniom_gradient(rtdb)
         if (status) then
            if (.not. rtdb_get(rtdb, 'oniom:energy', mt_dbl, 1, energy))
     $           call errquit('task_energy: no energy in rtdb?',0,
     &       RTDB_ERR)
            if (.not. rtdb_get(rtdb, 'oniom:gradient', mt_dbl, nelem, 
     $           gradient))
     $           call errquit('task_gradient: no gradient in rtdb?',0,
     &       RTDB_ERR)
            if (rtdb_get(rtdb, 'oniom:dipole', mt_dbl, 3, dipole)) then
               if (.not. rtdb_put(rtdb,'task:dipole',mt_dbl,3,dipole))
     $              call errquit('task_energy: failed storing dipole',0,
     &       RTDB_ERR)
            else
               ignore = rtdb_delete(rtdb, 'task:dipole')
            endif
         endif
      endif
c
c     TP: moved omm above ocons
      if(omm.and.(.not.numerical)) 
     $    call mm_add_egrad(rtdb,energy,nelem/3,gradient)
       ocons = cons_active(rtdb)
c      if(.not.rtdb_get(rtdb, 'cons:active',mt_log, 1, ocons)) then
c        ocons = .false.
c      end if
      if(ocons)  call cons_add_egrad(rtdb,energy,nelem,gradient)
c
      call bq_force(rtdb)
c
c     if(omm.and.(.not.numerical)) 
c    $    call mm_add_egrad(rtdb,energy,nelem/3,gradient)
c
      if (status) then
#ifndef SKIP_ECCE
         call geom_print_rtdb_ecce(rtdb)
         call ecce_print1('total energy', mt_dbl, energy, 1)
         call ecce_print2('total gradient', mt_dbl, gradient,
     $        3,3,nelem/3)
#endif
         gmax = 0.0d0
         gnorm = 0.0d0
         do i = 1, nelem/3
            do k = 1, 3
               gnorm = gnorm + gradient(k,i)**2
               gmax  = max(abs(gradient(k,i)),gmax)
            end do
         enddo
         gnorm = sqrt(gnorm)
#ifndef SKIP_ECCE
         call ecce_print1('gradient norm', mt_dbl, gnorm, 1)
         call ecce_print1('gradient max', mt_dbl, gmax, 1)
#endif
      endif
c
      if (.not. rtdb_put(rtdb, 'task:energy', mt_dbl, 1, energy))
     $     call errquit('task_gradient: failed storing energy',0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:gradient', mt_dbl, nelem,
     $     gradient))
     $     call errquit('task_gradient: failed storing energy',0,
     &       RTDB_ERR)
c
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - wall
c
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .true.))
     $     call errquit('task_gradient: failed to set status',0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $     call errquit('task_gradient: failed storing cputime',0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $     call errquit('task_gradient: failed storing walltime',0,
     &       RTDB_ERR)
c
#ifndef SKIP_ECCE
      call ecce_print1('cpu time', mt_dbl, cpu, 1)
      call ecce_print1('wall time', mt_dbl, wall, 1)
      if (status) then
         call ecce_print_module_exit('task gradient', 'ok')
      else
         call ecce_print_module_exit('task gradient', 'failed')
      endif
#endif
c
      task_gradient = status
c
      end
C> @}
c
      logical function task_gradient_doit(rtdb,theory,energy,gradient)
*
* $Id$
*
      implicit none
#include "errquit.fh"
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
      logical semi_dir_mp2, rimp2g, pspw_gradient,band_gradient
      external semi_dir_mp2, rimp2g, pspw_gradient,band_gradient
      logical dft_energy_gradient, sodft_energy_gradient
      external dft_energy_gradient, sodft_energy_gradient
      logical  tddft_energy, tddft_gradient
      external tddft_energy, tddft_gradient
      logical tce_energy_fragment
      external tce_energy_fragment
      logical  paw_gradient
      external paw_gradient
      logical xtb_energy_gradient
      external xtb_energy_gradient
c
      logical mcscf_gradient, scf_gradient, grad_dft, task_num_grad
      external mcscf_gradient, scf_gradient, grad_dft, task_num_grad
      logical geom_rtdb_ncent
      logical leps_gradient 
      external leps_gradient
c
c      logical tce_gradient
c      external tce_gradient
c
      integer nelem, ncent
      logical status, numerical, ocosmo, xc_chktau
      external xc_chktau
      double precision dipole(3)
      character*80 key, prefix
c
      logical dftmp2
      logical rodft
c
      if (.not. geom_rtdb_ncent(rtdb, 'geometry', ncent))
     $     call errquit('task_gradient: failed getting #atoms',0,
     &       RTDB_ERR)
      nelem = ncent*3
c
#ifndef SKIP_ECCE
      call ecce_print_module_entry('task gradient')
#endif
c
      task_gradient_doit = .false.
c
c     Set the theory in the database so that all is consistent for
c     the lower level modules.  Should not really have to do this,
c     but this value is used (appropriately and inappropriately) and
c     things are not yet disentangled.
c
      if (.not. rtdb_cput(rtdb, 'task:theory', 1, theory))
     $     call errquit('task_gradient: setting theory?',0,
     &       RTDB_ERR)
c
c     Right now only have a QM component
c
      if (.not. rtdb_get(rtdb, 'task:numerical', mt_log, 1, 
     $     numerical)) numerical = .false.
c
c     Check if analytic gradients are supported within the theory
c
c     == double hybrid DFT ==
      dftmp2 = .false.
      if(theory.eq.'dft') then
        if (.not. rtdb_get(rtdb,'dft:dftmp2',mt_log,1,dftmp2))
     1   dftmp2 = .false.
        if (dftmp2) then
         numerical = .true.  ! do numerical gradients for fifth-rung (w/MP2) functionals
         if (ga_nodeid().eq.0) write(LuOut,*) 
     1    'DFT+MP2 based functional - numerical derivatives only'
        end if
      end if
c
c     == restricted open-shell DFT ==
      rodft = .false.
      if(theory.eq.'dft') then
        if (.not. rtdb_get(rtdb,'dft:rodft',mt_log,1,rodft))
     1   rodft = .false.
      end if
c
c     Currently cosmo does not work with analytic gradients
c
      if (.not. rtdb_get(rtdb,'slv:cosmo',mt_log,1,ocosmo))
     $     ocosmo = .false.
c
c     Actually do the deed
c
      prefix = theory           ! Most common
c
      if (theory .eq. 'scf') then
         if (numerical) then
           status = task_num_grad(rtdb)
         else
           status = scf_gradient(rtdb)
         endif
      else if (theory .eq. 'dft') then
         if (numerical) then
           status=task_num_grad(rtdb)
         else
           status = dft_energy_gradient(rtdb)
         endif
      else if (theory .eq. 'sodft') then
         if (numerical) then
           status=task_num_grad(rtdb)
         else
           status = sodft_energy_gradient(rtdb)
         endif
      else if (theory .eq. 'xtb') then
         if (numerical) then
           status=task_num_grad(rtdb)
         else
           status = xtb_energy_gradient(rtdb)
         endif
      else if (theory .eq. 'mp2' .or. theory .eq. 'semi_dir_mp2') then
         if (numerical) then
           status=task_num_grad(rtdb)
         else
           prefix = 'mp2'
           if(.not. rtdb_cput(rtdb,'sdmp2:task',1,'gradient'))
     $        call errquit('task_gradient: rtdb_put failed',0, RTDB_ERR)
           status = semi_dir_mp2(rtdb)
         endif
      else if (theory .eq. 'mcscf') then
         if (numerical) then
           status=task_num_grad(rtdb)
         else
           status = mcscf_gradient(rtdb)
         endif
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
      else if (theory .eq. 'tddft') then
         if (numerical) then
            status = task_num_grad(rtdb) ! Only numerical
         else
            status = tddft_gradient(rtdb)
         endif
      else if (theory .eq. 'tce') then
         status = task_num_grad(rtdb) ! Only numerical
c         if (numerical) then 
c            status = task_num_grad(rtdb) ! Only numerical
c         else
c            status = tce_gradient(rtdb)
c         endif
      else if (theory .eq. 'uccd') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'ccd'))
     1   call errquit('task_gradient: failed writing to rtdb',0,
     &       RTDB_ERR)
         status = task_num_grad(rtdb) ! Only numerical
      else if (theory .eq. 'ulccd') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'lccd'))
     1   call errquit('task_gradient: failed writing to rtdb',0,
     &       RTDB_ERR)
         status = task_num_grad(rtdb) ! Only numerical
      else if (theory .eq. 'uccsd') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'ccsd'))
     1   call errquit('task_gradient: failed writing to rtdb',0,
     &       RTDB_ERR)
         status = task_num_grad(rtdb) ! Only numerical
      else if (theory .eq. 'ulccsd') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'lccsd'))
     1   call errquit('task_gradient: failed writing to rtdb',0,
     &       RTDB_ERR)
         status = task_num_grad(rtdb) ! Only numerical
      else if (theory .eq. 'uqcisd') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'qcisd'))
     1   call errquit('task_gradient: failed writing to rtdb',0,
     &       RTDB_ERR)
         status = task_num_grad(rtdb) ! Only numerical
      else if (theory .eq. 'uccsdt') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'ccsdt'))
     1   call errquit('task_gradient: failed writing to rtdb',0,
     &       RTDB_ERR)
         status = task_num_grad(rtdb) ! Only numerical
      else if (theory .eq. 'uccsdtq') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'ccsdtq'))
     1   call errquit('task_gradient: failed writing to rtdb',0,
     &       RTDB_ERR)
         status = task_num_grad(rtdb) ! Only numerical
      else if (theory .eq. 'ucisd') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'cisd'))
     1   call errquit('task_gradient: failed writing to rtdb',0,
     &       RTDB_ERR)
         status = task_num_grad(rtdb) ! Only numerical
      else if (theory .eq. 'ucisdt') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'cisdt'))
     1   call errquit('task_gradient: failed writing to rtdb',0,
     &       RTDB_ERR)
         status = task_num_grad(rtdb) ! Only numerical
      else if (theory .eq. 'ucisdtq') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'cisdtq'))
     1   call errquit('task_gradient: failed writing to rtdb',0,
     &       RTDB_ERR)
         status = task_num_grad(rtdb) ! Only numerical
      else if (theory .eq. 'ump2') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'mbpt2'))
     1   call errquit('task_gradient: failed writing to rtdb',0,
     &       RTDB_ERR)
         status = task_num_grad(rtdb) ! Only numerical
      else if (theory .eq. 'ump3') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'mbpt3'))
     1   call errquit('task_gradient: failed writing to rtdb',0,
     &       RTDB_ERR)
         status = task_num_grad(rtdb) ! Only numerical
      else if (theory .eq. 'ump4') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'mbpt4'))
     1   call errquit('task_gradient: failed writing to rtdb',0,
     &       RTDB_ERR)
         status = task_num_grad(rtdb) ! Only numerical
      else if (theory .eq. 'umbpt2') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'mbpt2'))
     1   call errquit('task_gradient: failed writing to rtdb',0,
     &       RTDB_ERR)
         status = task_num_grad(rtdb) ! Only numerical
      else if (theory .eq. 'umbpt3') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'mbpt3'))
     1   call errquit('task_gradient: failed writing to rtdb',0,
     &       RTDB_ERR)
         status = task_num_grad(rtdb) ! Only numerical
      else if (theory .eq. 'umbpt4') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'mbpt4'))
     1   call errquit('task_gradient: failed writing to rtdb',0,
     &       RTDB_ERR)
         status = task_num_grad(rtdb) ! Only numerical
      else if (theory.eq.'pspw') then
         status = pspw_gradient(rtdb)
      else if (theory .eq. 'band') then
         status = band_gradient(rtdb)
      else if (theory .eq. 'paw') then
         status = paw_gradient(rtdb)
      else if (theory .eq. 'qmmm') then
         if (ga_nodeid().eq.0) then
           write(6,*)'task_gradient: qmmm should never be the theory '
           write(6,*)'task_gradient: you probably reversed the order ',
     &         'of the task input '
           write(6,*)'task [qmmm] <string theory> [<string operation ',
     &         '= energy>] [numerical] [ignore]'
         endif
         status = .false.
c      else if (theory .eq. 'embed') then
c         if (ga_nodeid().eq.0) then
c           write(6,*)'task [embed] <string theory> [<string operation ',
c     &         '= energy>] [numerical] [ignore]'
c         endif
c         status = .false.
      else if (theory .eq. 'oniom') then
         status = task_num_grad(rtdb) ! Doing Oniom in a numerical fashion
      else if (theory .eq. 'dangchang') then
         status = task_num_grad(rtdb) ! Dang Chang only supports numerical gradients
      else if (theory .eq. 'leps') then
         if (numerical) then
           status = task_num_grad(rtdb)
         else
           status = leps_gradient(rtdb)
         endif
      else
         call util_flush(luout)
         write(luout,*)'task_gradient: theory = ',
     &                 theory(1:inp_strlen(theory))
         call util_flush(luout)
         call errquit('task_gradient: unknown theory',0, INPUT_ERR)
      endif
c
      if (status) then
         key = ' '
         write(key,'(a,'':energy'')') prefix(1:inp_strlen(prefix))
         if (.not. rtdb_get(rtdb, key, mt_dbl, 1, energy))
     $        call errquit('task_gradient: no energy in rtdb?',0,
     &       RTDB_ERR)
c
         key = ' '
         write(key,'(a,'':dipole'')') prefix(1:inp_strlen(prefix))
         if (rtdb_get(rtdb, key, mt_dbl, 3, dipole)) then
           if (.not. rtdb_put(rtdb, 'task:dipole', mt_dbl, 3, dipole))
     $        call errquit('task_gradient: failed storing dipole',0,
     &       RTDB_ERR)
	 endif
c
         write(key,'(a,'':gradient'')') prefix(1:inp_strlen(prefix))
         if (.not. rtdb_get(rtdb, key, mt_dbl, nelem, gradient))
     $        call errquit('task_gradient: no gradient in rtdb?',0,
     &       RTDB_ERR)
      end if
c
      task_gradient_doit = status
c
      end

