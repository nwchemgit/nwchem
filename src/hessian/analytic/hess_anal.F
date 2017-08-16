C
C $Id$
C
      logical function hess_anal(rtdb)
C
C  This function is the molecular calculation module for calculating
C  analytic second derivatives.  Specifically, the second derivative
C  of the energy with respect to the nuclear coordinates will be
C  calculated.  Everything must be either read from or written to
C  the RTDB or files.
c
c  This routine assumes that task_hessian is the calling routine.  It
c  assumes that checks for unsupported wavefunctions has been done.
C
C  The notation here will follow that of Johnson and Frisch, J.Chem.Phys.
C  100, 7429-7442 (1994) unless otherwise annotated.
C
      implicit none
c
#include "errquit.fh"
#include "bas.fh"
#include "geom.fh"
#include "global.fh"
#include "hess_info.fh"
#include "mafdecls.fh"
#include "pstat.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "util.fh"
#include "prop.fh"
c
      integer rtdb

      integer restr, junk
      integer ps_total, ps_onedd, ps_nucdd, ps_twodd, ps_cphf
      logical oprintsome, oprint, ohessprof
      logical task_energy
      logical bgj_get_scf_method, hess_check
      external task_energy
      external bgj_get_scf_method, hess_check
      double precision t0
c
      logical ostart,ocontinue,orestart
      logical status
c
      hess_anal = .false.
c
c  Broadcast information to all of the nodes on read
c
      if (.not. rtdb_parallel(.true.)) call errquit
     *  ('hess_anal: unable to set rtdb_parallel to true',555, RTDB_ERR)
c
c  First check to see if we should fail the calculation because
c  the code can't handle it.
c
      if (.not. hess_check(rtdb)) then
        write(LuOut,*) 'The Analytic Hessian code cannot run ',
     *    'calculations open shell DFT, or charge fitted DFT ',
     *    'or relativistic all-electron methods,',
     *    'or metaGGA functionals '
        write(LuOut,*) 'Please use the numerical option'
        return
      endif
c
c  Set up the print mechanism (push down the stack and loads the print
c  information for hess)
c
      call util_print_push
      call util_print_rtdb_load(rtdb,'hess')
      oprintsome = util_print('hess-follow',print_low)
      oprintsome = oprintsome .and. (ga_nodeid().eq.0)
      oprint = util_print('hess_follow',print_high)
      oprint = oprint .and. (ga_nodeid().eq.0)
c
c  See if we are going to profile
c
      if (.not. rtdb_get(rtdb, 'hess:pstat',
     $     MT_LOG, 1, ohessprof)) ohessprof = .false.
      if (ohessprof) then
        if (.not.pstat_allocate('hess: total',pstat_qstat,0,junk,
     $    ps_total)) call errquit('hess_init: pstat_allocate failed',
     $    555, MEM_ERR)
        if (.not.pstat_allocate('hess: intdd-1e',pstat_qstat,0,junk,
     $    ps_onedd)) call errquit('hess_init: pstat_allocate failed',
     $    555, MEM_ERR)
        if (.not.pstat_allocate('hess: nuclear',pstat_qstat,0,junk,
     $    ps_nucdd)) call errquit('hess_init: pstat_allocate failed',
     $    555, MEM_ERR)
        if (.not.pstat_allocate('hess: intdd-2e',pstat_qstat,0,junk,
     $    ps_twodd)) call errquit('hess_init: pstat_allocate failed',
     $    555, MEM_ERR)
        if (.not.pstat_allocate('hess: cphf',pstat_qstat,0,junk,
     $    ps_cphf)) call errquit('hess_init: pstat_allocate failed',
     $    555, MEM_ERR)
        call pstat_on(ps_total)
      endif
c
c  Initial some variables for the energy calculation if they are
c  not sufficient for the Hessian calculation.
c
      call hess_energyinit(rtdb)
c
c  Find out if this is a restart
c
      restr = 0
      call util_get_rtdb_state(rtdb,ostart,ocontinue,orestart)
      if (orestart) restr = 1
      if (restr.gt.0.and. ga_nodeid().eq.0) 
     &    write(luout,*) 'Restarting calculation'
c
      call solver_setup(rtdb,restr)  ! restart parameters for solver (ga_lkain())
c
c  Get theory type
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('hess_init: theory not specified',555, RTDB_ERR)
      if (theory.eq.'hyb') theory = 'dft'
      if (theory.eq.'scf') theory = 'hf'
      use_theory = theory
c
c  Check to make sure that the energy and wavefunction are calculated
c
      if(.not.task_energy(rtdb)) call errquit
     *    ('hess_anal: energy failure', 555, CALC_ERR)
      if (oprint) write(LuOut,*) 'HESS: after the task_energy call'
c
c  2nd check to see if we should fail the calculation because
c  the code can't handle it.
c
      if (.not. hess_check(rtdb)) then
        write(LuOut,*) 'The Analytic Hessian code cannot run ',
     *    'calculations open shell DFT, or charge fitted DFT ',
     *    'or relativistic all-electron methods,',
     *    'or metaGGA functionals '
        write(LuOut,*) 'Please use the numerical option'
        return
      endif
c
c     initialize hessian variables
c
      if (oprint) write(LuOut,*) 'HESS: before the hess_init call'
      call hess_init(rtdb)
c
c  Get the one-electron second derivative contributions
c
      if (oprint) write(LuOut,*) 'HESS: before the onedd_cont call'
      if (ohessprof) call pstat_on(ps_onedd)
      t0=util_wallsec()
      call onedd_cont(rtdb)
      if (ohessprof) call pstat_off(ps_onedd)
c
c  Get the nuclear contribution to the hessian
c
      if (oprint) write(LuOut,*) 'HESS: before the nucdd_cont call'
      if (ohessprof) call pstat_on(ps_nucdd)
      call nucdd_cont(rtdb)
      if (ohessprof) call pstat_off(ps_nucdd)
      call hess_restart(rtdb,1)
        
      if (oprintsome) then
           write(LuOut,11) util_wallsec()-t0
           call util_flush(luout)
      endif
 11   format(/' HESSIAN: the one electron ',
     $    'contributions are done in ',f9.1,'s'/)
c
c  Get the two-electron second derivative contributions
c
      if (oprint) write(LuOut,*) 'HESS: before the twodd_cont call'
      if (ohessprof) call pstat_on(ps_twodd)
      t0=util_wallsec()
      call twodd_cont(rtdb)
      if (ohessprof) call pstat_off(ps_twodd)
      call hess_restart(rtdb,2)
      if (oprintsome) write(LuOut,12) util_wallsec()-t0
 12     format(/' HESSIAN: the two electron ',
     $    'contributions are done in ',f9.1,'s'/)
c
c  Copy fock matrix bits to g_rhs
c
      call hess_copymat(ndens,n3xyz,nbf,nbf,g_rhs,g_rhs_fx)
c
c  solve the CPSCF problem
c
      if (oprint) write(LuOut,*) 'HESS: before the hess_cphf call'
      if (ohessprof) call pstat_on(ps_cphf)
      call hess_cphf(rtdb,restr)
      if (ohessprof) call pstat_off(ps_cphf)
      if (oprintsome) write(LuOut,*) 'HESSIAN: the CPHF ',
     $  'contributions are done'
      if(oprint) call hess_print
c
c  dump hessian with the final result
c
      if(ga_nodeid().eq.0) then
         call hess_wrt()
         call stpr_print_ddipole(dbl_mb(k_dder),
     P        ' Derivative Dipole ',
     A       ncent,-1.0d00)
      endif
c     
c  clean up any hessian specific stuff
c
      if (oprint) write(LuOut,*) 'HESS: before the hess_tidy call'
      call hess_tidy()
c
c  clean up the changes we made because of the energy calculation
c
      if (oprint) write(LuOut,*) 
     $  'HESS: before the hess_energy_tidy call'
      call hess_energytidy(rtdb)
c
      if (ohessprof) then
        call pstat_off(ps_total)
        call pstat_print_all
        if (.not. pstat_free(ps_total))
     $    call errquit('hess_anal: pstat_free error',555, MEM_ERR)
        if (.not. pstat_free(ps_onedd))
     $    call errquit('hess_anal: pstat_free error',555, MEM_ERR)
        if (.not. pstat_free(ps_nucdd))
     $    call errquit('hess_anal: pstat_free error',555, MEM_ERR)
        if (.not. pstat_free(ps_twodd))
     $    call errquit('hess_anal: pstat_free error',555, MEM_ERR)
        if (.not. pstat_free(ps_cphf))
     $    call errquit('hess_anal: pstat_free error',555, MEM_ERR)
      endif
      call util_print_pop
c
      if (oprintsome) write(LuOut,*) 'HESSIAN: the Hessian ',
     $  'is done'
      hess_anal = .true.
c
      return 
      end
