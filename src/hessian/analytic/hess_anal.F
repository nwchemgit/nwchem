C
C $Id: hess_anal.F,v 1.10 2001-04-05 22:25:06 windus Exp $
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
#include "bas.fh"
#include "geom.fh"
#include "global.fh"
#include "hess_info.fh"
#include "mafdecls.fh"
#include "pstat.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "util.fh"
      integer rtdb

      integer restr, junk
      integer ps_total, ps_onedd, ps_nucdd, ps_twodd, ps_cphf
      logical oprint, ohessprof
      logical task_energy
      logical bgj_get_scf_method
      external task_energy
      external bgj_get_scf_method
c
      hess_anal = .false.
c
c     find out if this is a restart
c
      if (.not. rtdb_get(rtdb,'hessian:restart',mt_int,1,restr))
     *  restr = 0
c
c  Broadcast information to all of the nodes on read
c
      if (.not. rtdb_parallel(.true.)) call errquit
     *  ('hess_anal: unable to set rtdb_parallel to true',555)
c
c  Set up the print mechanism (push down the stack and loads the print
c  information for hess)
c
      call util_print_push
      call util_print_rtdb_load(rtdb,'hess')
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
     $    555)
        if (.not.pstat_allocate('hess: intdd-1e',pstat_qstat,0,junk,
     $    ps_onedd)) call errquit('hess_init: pstat_allocate failed',
     $    555)
        if (.not.pstat_allocate('hess: nuclear',pstat_qstat,0,junk,
     $    ps_nucdd)) call errquit('hess_init: pstat_allocate failed',
     $    555)
        if (.not.pstat_allocate('hess: intdd-2e',pstat_qstat,0,junk,
     $    ps_twodd)) call errquit('hess_init: pstat_allocate failed',
     $    555)
        if (.not.pstat_allocate('hess: cphf',pstat_qstat,0,junk,
     $    ps_cphf)) call errquit('hess_init: pstat_allocate failed',
     $    555)
        call pstat_on(ps_total)
      endif
c
c  Check to make sure that the energy and wavefunction are calculated
c
      if(.not.task_energy(rtdb)) call errquit 
     *  ('hess_anal: energy failure', 555)
      if (oprint) write(LuOut,*) 'HESS: after the task_energy call'
c
c     initialize hessian variables
c
      if (.not.bgj_get_scf_method(theory)) call errquit
     *  ('hess_anal: failure to get theory type', 555)
      if (theory.eq.'hyb') theory = 'dft'
      if (oprint) write(LuOut,*) 'HESS: before the hess_init call'
      call hess_init(rtdb)
c
c  Get the one-electron second derivative contributions
c
      if (oprint) write(LuOut,*) 'HESS: before the onedd_cont call'
      if (ohessprof) call pstat_on(ps_onedd)
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
c
c  Get the two-electron second derivative contributions
c
      if (oprint) write(LuOut,*) 'HESS: before the twodd_cont call'
      if (ohessprof) call pstat_on(ps_twodd)
      call twodd_cont(rtdb)
      if (ohessprof) call pstat_off(ps_twodd)
      call hess_restart(rtdb,2)
c
c  solve the CPSCF problem
c
      if (oprint) write(LuOut,*) 'HESS: before the hess_cphf call'
      if (ohessprof) call pstat_on(ps_cphf)
      call hess_cphf(rtdb)
      if (ohessprof) call pstat_off(ps_cphf)
c
c  dump hessian with the final result
c
      call hess_wrt()
c
c  clean up any hessian specific stuff
c
      if (oprint) write(LuOut,*) 'HESS: before the hess_tidy call'
      call hess_tidy()
c
      if (ohessprof) then
        call pstat_off(ps_total)
        call pstat_print_all
        if (.not. pstat_free(ps_total))
     $    call errquit('hess_anal: pstat_free error',555)
        if (.not. pstat_free(ps_onedd))
     $    call errquit('hess_anal: pstat_free error',555)
        if (.not. pstat_free(ps_nucdd))
     $    call errquit('hess_anal: pstat_free error',555)
        if (.not. pstat_free(ps_twodd))
     $    call errquit('hess_anal: pstat_free error',555)
        if (.not. pstat_free(ps_cphf))
     $    call errquit('hess_anal: pstat_free error',555)
      endif
      call util_print_pop
c
      if (oprint) write(LuOut,*) 'HESS: DONE'
      hess_anal = .true.
      return 
      end
