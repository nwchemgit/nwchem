C> \ingroup task
C> @{
C>
C> \brief The main driver for single point energy calculations
C>
C> Analyse the theory options on the `task <theory> [energy]` line.
C> Based on the theory found other drivers may be called. If the
C> energy calculation was successfull, i.e. in most computations the
C> calculation converged, then we return .true.
C>
C> Theories for which special energy drivers are called are: ONIOM,
C> BSSE, QMMM, and RISM.
C>
C> \return Return .true. if the energy was successfully calculated,
C> and .false. otherwise
C>
      logical function task_energy(rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
      integer rtdb !< [Input] The RTDB handle
c
c     Fortran does not like recursive calls, so for ease of coding
c     ONIOM do ONIOM out here so that it call call task_energy_doit().
c
      logical status, ignore
      logical ocons
      character*32 theory
      double precision cpu, wall, energy, dipole(3)
      logical task_energy_doit, oniom_energy
      external task_energy_doit, oniom_energy
      logical orism
      logical omm
      logical task_qmmm
      logical task_rism
      external task_rism
      logical task_fragment
      external task_qmmm_energy
      logical  task_qmmm_energy
      external task_fragment_energy
      logical  task_fragment_energy
      logical cons_active
      external cons_active
      logical  task_ccca_energy
      external task_ccca_energy
c     TP: added variable
      character*30 operation
c
c      logical task_embed
c      logical  embed_energy
c      external embed_energy
cc
      logical task_bsse
      logical bsse_energy
      external bsse_energy
      logical oecce
cc
      oecce = .true.
      cpu  = util_cpusec()
      wall = util_wallsec()
c
      ignore = rtdb_delete(rtdb, 'task:energy')
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_energy: failed to invalidate status',0,
     &       RTDB_ERR)
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task:energy: theory not specified',0, RTDB_ERR)
c
c      if( .not. rtdb_get(rtdb,'task:fragment',mt_log,1,task_fragment))
c     &    task_fragment = .false.
c
      if( .not. rtdb_get(rtdb,'task:rism',mt_log,1,orism))
     &    orism = .false.
c
      if( .not. rtdb_get(rtdb,'task:mm',mt_log,1,omm))
     &    omm = .false.
c
      if( .not. rtdb_get(rtdb,'task:QMMM',mt_log,1,task_qmmm))
     &    task_qmmm = .false.
c
c      if( .not. rtdb_get(rtdb,'task:embed',mt_log,1,task_embed))
c     &    task_embed = .false.
cc
      if (.not. rtdb_get(rtdb,'bsse',mt_log,1,task_bsse))
     $     task_bsse = .false.
cc 
      if(task_qmmm) then
          status = task_qmmm_energy(rtdb)
          if (.not. rtdb_get(rtdb, 'qmmm:energy', mt_dbl, 1, energy))
     $           call errquit('task_energy: no energy in rtdb?',0,
     &       RTDB_ERR)
       else if(orism) then
          status = task_rism(rtdb)
          energy=0.0d0
c      else if(task_fragment) then
c          write(*,*) "found task fragment"
c          status = task_fragment_energy(rtdb)
c      else if (task_embed) then
c          status = embed_energy(rtdb)
c          if (.not. rtdb_get(rtdb, 'embed:energy', mt_dbl, 1, energy))
c     $           call errquit('task_energy: no energy in rtdb?',0,
c     &       RTDB_ERR)
      else if (theory .ne. 'oniom') then
cc
        if(task_bsse) then
            status = bsse_energy(rtdb,theory,energy)
            if(.not.status)
     $        call errquit('task_energy: bsse_energy fail',0,CALC_ERR)
         else if (theory.eq.'ccca') then
            status = task_ccca_energy(rtdb,theory,energy)
         else
c          TP: added the case where mm is present in task directive
           if(omm) then
             call mm_task_energy(rtdb,theory,energy,status)
          else
            status = task_energy_doit(rtdb,theory,energy)
          endif ! if(omm)
        endif
cc
      else
         status = oniom_energy(rtdb)
         if (status) then
            if (.not. rtdb_get(rtdb, 'oniom:energy', mt_dbl, 1, energy))
     $           call errquit('task_energy: no energy in rtdb?',0,
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
       ocons = cons_active(rtdb)
c      if(.not.rtdb_get(rtdb, 'cons:active', mt_log, 1, ocons)) then
c        ocons=.false.
c      end if
c     TP: moved omm above ocons
      if(omm) call mm_add_energy(rtdb,energy)
      if(ocons) call cons_add_energy(rtdb,energy)
c     if(omm) call mm_add_energy(rtdb,energy)

c     write(*,*) "putting energy into rtdb",energy
      if (.not. rtdb_put(rtdb, 'task:energy', mt_dbl, 1, energy))
     $     call errquit('task_energy: failed storing energy',0,
     &       RTDB_ERR)
c
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - wall
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, status))
     $     call errquit('task_energy: failed to store status',0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $     call errquit('task_energy: failed storing cputime',0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $     call errquit('task_energy: failed storing walltime',0,
     &       RTDB_ERR)
c
#ifndef SKIP_ECCE
      call ecce_print1('cpu time', mt_dbl, cpu, 1)
      call ecce_print1('wall time', mt_dbl, wall, 1)
      if (theory .eq. 'md')  oecce=.false.
      if (theory .eq. 'smd') oecce=.false.
      if (orism) oecce=.false.
      if (status) then
         if (oecce) call geom_print_rtdb_ecce(rtdb)
         call ecce_print1('total energy', mt_dbl, energy, 1)
         call ecce_print_module_exit('task energy', 'ok')
      else
         call ecce_print_module_exit('task energy', 'failed')
      endif
#endif
c         
      task_energy = status
c
      end
C> @}
      logical function task_energy_doit(rtdb,theory,energy)
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
      integer rtdb
      character*(*) theory
      double precision energy
      
c
c     Generic NWChem interface to compute the energy.  Currently
c     only the QM components are supported.
c
c     RTDB input parameters
c     ---------------------
c     task:theory (string) - name of (QM) level of theory to use
c     
c     RTDB output parameters
c     ----------------------
c     task:status (logical)- T/F for success/failure
c     if (status) then
c     .  task:energy (real)   - total energy
c     .  task:dipole(real(3)) - total dipole moment if available
c     .  task:cputime (real)  - cpu time to execute the task
c     .  task:walltime (real) - wall time to execute the task
c
c     Also returns status through the function value
c
      logical  scf, rimp2, direct_mp2, ccsd, mcscf, selci
      external scf, rimp2, direct_mp2, ccsd, mcscf, selci
c     logical  oimp2
c     external oimp2
      logical semi_dir_mp2, rimp2g,nwmd,pspw_energy,band_energy
      external semi_dir_mp2, rimp2g,nwmd,pspw_energy,band_energy
      logical dft_energy, sodft_energy, tddft_energy
      external dft_energy, sodft_energy, tddft_energy
      logical tce_energy_fragment
      external tce_energy_fragment
      logical  paw_energy,ofpw_energy
      external paw_energy,ofpw_energy
      logical task_smd_energy
      external task_smd_energy
      logical dc_energy
      external dc_energy
      logical leps_energy
      external leps_energy
      logical  wfn1_task_energy
      external wfn1_task_energy
      logical xtb_energy
      external xtb_energy
c
      logical status
      double precision dipole(3)
      character*80 key, prefix
c
      logical dftmp2
      integer which_mp2
c
#ifndef SKIP_ECCE
      call ecce_print_module_entry('task energy')
#endif
c
      task_energy_doit = .false.
      status = .false.
c
c     Set the theory in the database so that all is consistent for
c     the lower level modules.  Should not really have to do this,
c     but this value is used (appropriately and inappropriately) and
c     things are not yet disentangled.
c
      if (.not. rtdb_cput(rtdb, 'task:theory', 1, theory))
     $     call errquit('task_energy: setting theory?',0,
     &       RTDB_ERR)
c
c     Actually do the deed
c
      prefix = theory           ! Most common scenario
      if (theory .eq. 'scf') then
         status = scf(rtdb)
      else if (theory .eq. 'smd') then
         status = task_smd_energy(rtdb)
      else if (theory .eq. 'dft') then
         if (.not. rtdb_get(rtdb,'dft:dftmp2',mt_log,1,dftmp2))
     1       dftmp2 = .false.
         if (dftmp2) then
           if (.not.rtdb_get(rtdb,'dft:which_mp2',mt_int,1,which_mp2))
     1       which_mp2 = 1
           if (which_mp2.eq.1) then        ! semidirect MP2
              prefix = 'mp2'
              if(.not. rtdb_cput(rtdb,'sdmp2:task',1,'energy')) then
                 call errquit('task_energy: rtdb_put failed',37,
     1                        RTDB_ERR)
              endif
              status = semi_dir_mp2(rtdb)
           elseif (which_mp2.eq.2) then    ! direct MP2
              status = direct_mp2(rtdb)
           elseif (which_mp2.eq.3) then    ! RI-MP2
              status = rimp2(rtdb)
           else
             write(6,*) 'task_energy: problem with DFT+MP2'
           endif
         else
            status = dft_energy(rtdb)
         endif
      else if (theory .eq. 'xtb')then
         status = xtb_energy(rtdb)
      else if (theory .eq. 'sodft')then 
         status = sodft_energy(rtdb)
      else if (theory .eq. 'tddft') then
         status = tddft_energy(rtdb)
      else if (theory .eq. 'tce') then
         status = tce_energy_fragment(rtdb)
      else if (theory .eq. 'uccd') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'ccd'))
     1   call errquit('task_energy: failed writing to rtdb',0,
     &       RTDB_ERR)
         status = tce_energy_fragment(rtdb)
      else if (theory .eq. 'ulccd') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'lccd'))
     1   call errquit('task_energy: failed writing to rtdb',0, RTDB_ERR)
         status = tce_energy_fragment(rtdb)
      else if (theory .eq. 'uccsd') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'ccsd'))
     1   call errquit('task_energy: failed writing to rtdb',0, RTDB_ERR)
         status = tce_energy_fragment(rtdb)
      else if (theory .eq. 'ulccsd') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'lccsd'))
     1   call errquit('task_energy: failed writing to rtdb',0, RTDB_ERR)
         status = tce_energy_fragment(rtdb)
      else if (theory .eq. 'uqcisd') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'qcisd'))
     1   call errquit('task_energy: failed writing to rtdb',0, RTDB_ERR)
         status = tce_energy_fragment(rtdb)
      else if (theory .eq. 'uccsdt') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'ccsdt'))
     1   call errquit('task_energy: failed writing to rtdb',0, RTDB_ERR)
         status = tce_energy_fragment(rtdb)
      else if (theory .eq. 'uccsdtq') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'ccsdtq'))
     1   call errquit('task_energy: failed writing to rtdb',0, RTDB_ERR)
         status = tce_energy_fragment(rtdb)
      else if (theory .eq. 'ucisd') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'cisd'))
     1   call errquit('task_energy: failed writing to rtdb',0, RTDB_ERR)
         status = tce_energy_fragment(rtdb)
      else if (theory .eq. 'ucisdt') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'cisdt'))
     1   call errquit('task_energy: failed writing to rtdb',0, RTDB_ERR)
         status = tce_energy_fragment(rtdb)
      else if (theory .eq. 'ucisdtq') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'cisdtq'))
     1   call errquit('task_energy: failed writing to rtdb',0, RTDB_ERR)
         status = tce_energy_fragment(rtdb)
      else if (theory .eq. 'ump2') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'mbpt2'))
     1   call errquit('task_energy: failed writing to rtdb',0, RTDB_ERR)
         status = tce_energy_fragment(rtdb)
      else if (theory .eq. 'ump3') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'mbpt3'))
     1   call errquit('task_energy: failed writing to rtdb',0, RTDB_ERR)
         status = tce_energy_fragment(rtdb)
      else if (theory .eq. 'ump4') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'mbpt4'))
     1   call errquit('task_energy: failed writing to rtdb',0, RTDB_ERR)
         status = tce_energy_fragment(rtdb)
      else if (theory .eq. 'umbpt2') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'mbpt2'))
     1   call errquit('task_energy: failed writing to rtdb',0, RTDB_ERR)
         status = tce_energy_fragment(rtdb)
      else if (theory .eq. 'umbpt3') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'mbpt3'))
     1   call errquit('task_energy: failed writing to rtdb',0, RTDB_ERR)
         status = tce_energy_fragment(rtdb)
      else if (theory .eq. 'umbpt4') then
         prefix = 'tce'
         if (.not.rtdb_cput(rtdb,'tce:model',1,'mbpt4'))
     1   call errquit('task_energy: failed writing to rtdb',0, RTDB_ERR)
         status = tce_energy_fragment(rtdb)
      else if (theory .eq. 'rimp2') then
         status = rimp2(rtdb)
      else if (theory .eq. 'rimp2_grad') then
         status = rimp2g(rtdb)
c     else if (theory .eq. 'oimp2') then
c        status = oimp2(rtdb)
      else if (theory .eq. 'direct_mp2') then
cedo         prefix = 'mp2'
         status = direct_mp2(rtdb)
      else if (theory .eq. 'ccsd') then
         status = ccsd(rtdb)
      else if (theory .eq. 'ccsd(t)') then
         prefix = 'ccsd'
         status = ccsd(rtdb)
      else if (theory .eq. 'ccsd+t(ccsd)') then
         prefix = 'ccsd'
         status = ccsd(rtdb)
      else if (theory .eq. 'mcscf') then
         status = mcscf(rtdb)
      else if (theory .eq. 'selci') then
         status = selci(rtdb)
      else if (theory .eq. 'mp2' .or. theory.eq.'semi_dir_mp2') then
         prefix = 'mp2'
         if(.not. rtdb_cput(rtdb,'sdmp2:task',1,'energy'))
     $        call errquit('task_energy: rtdb_put failed',0, RTDB_ERR)
         status = semi_dir_mp2(rtdb)
      else if (theory .eq. 'pspw') then
         status = pspw_energy(rtdb)
      else if (theory .eq. 'band') then
         status = band_energy(rtdb)
      else if (theory .eq. 'paw') then
         status = paw_energy(rtdb)
      else if (theory .eq. 'ofpw') then
         status = ofpw_energy(rtdb)
      else if (theory .eq. 'md') then
         status = nwmd(rtdb)
      else if (theory .eq. 'qmmm') then
         if(ga_nodeid().eq.0) then
           write(6,*)'task_energy: qmmm should never be the theory '
           write(6,*)'task_energy: you probably reversed the order ',
     &         'of the task input '
           write(6,*)'task [qmmm] <string theory> [<string operatio',
     &         'n = energy>] [numerical || analytic] [ignore]'
         endif
         status = .false.
c      else if (theory .eq. 'embed') then
c         if(ga_nodeid().eq.0) then
c           write(6,*)'task [embed] <string theory> [<string operatio',
c     &         'n = energy>] [numerical || analytic] [ignore]'
c         endif
c         status = .false.
      else if (theory .eq. 'leps') then
         status = leps_energy(rtdb)
      else if (theory .eq. 'dangchang') then
         status = dc_energy(rtdb)
      else if (theory .eq. 'wfn1') then
         status = wfn1_task_energy(rtdb)
      else
         call errquit('task_energy: unknown theory',0, INPUT_ERR)
      endif
c
      if (status) then
c
c     ENERGY is returned via the subroutine interface.
c     DIPOLE is currently quietly stored.
c
         key = ' '
         write(key,'(a,'':energy'')') prefix(1:inp_strlen(prefix))
         if (.not. rtdb_get(rtdb, key, mt_dbl, 1, energy))
     $        call errquit('task_energy: no energy in rtdb?',0,
     &       RTDB_ERR)
c
         key = ' '
         write(key,'(a,'':dipole'')') prefix(1:inp_strlen(prefix))
         if (rtdb_get(rtdb, key, mt_dbl, 3, dipole)) then
           if (.not. rtdb_put(rtdb, 'task:dipole', mt_dbl, 3, dipole))
     $          call errquit('task_energy: failed storing dipole',0,
     &       RTDB_ERR)
	 endif
      endif
c
      task_energy_doit = status
c
      end
      subroutine rtdb_clean(rtdb, prefix)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "inp.fh"
      integer rtdb
      character*(*) prefix
c     
c     Delete all entries in the database beginning with prefix
c
      integer lp
      logical status
      character*255 name
c
      lp = inp_strlen(prefix)
c
      status = rtdb_first(rtdb, name)
 10   if (status) then
         if (inp_compare(.true.,prefix(1:lp),name(1:lp))) then
            if (.not. rtdb_delete(rtdb,name)) call errquit
     $           ('rtdb_clean: failed deleting existing entry',0,
     &       RTDB_ERR)
         endif
         status = rtdb_next(rtdb, name)
         goto 10
      endif
c     
      end

