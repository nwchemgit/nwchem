      logical function task_energy(rtdb)
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
      integer rtdb
c
c     Fortran does not like recursive calls, so for ease of coding
c     ONIOM do ONIOM out here so that it call call task_energy_doit().
c
      logical status, ignore
      character*32 theory
      double precision cpu, wall, energy, dipole(3)
      logical task_energy_doit, oniom_energy
      external task_energy_doit, oniom_energy
c
      cpu  = util_cpusec()
      wall = util_wallsec()
c
      ignore = rtdb_delete(rtdb, 'task:energy')
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_energy: failed to invalidate status',0)
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task:energy: theory not specified',0)
c
      if (theory .ne. 'oniom') then
         status = task_energy_doit(rtdb,theory,energy)
      else
         status = oniom_energy(rtdb)
         if (status) then
            if (.not. rtdb_get(rtdb, 'oniom:energy', mt_dbl, 1, energy))
     $           call errquit('task_energy: no energy in rtdb?',0)
            if (rtdb_get(rtdb, 'oniom:dipole', mt_dbl, 3, dipole)) then
               if (.not. rtdb_put(rtdb,'task:dipole',mt_dbl,3,dipole))
     $              call errquit('task_energy: failed storing dipole',0)
            else
               ignore = rtdb_delete(rtdb, 'task:dipole')
            endif
         endif
      endif
      if (.not. rtdb_put(rtdb, 'task:energy', mt_dbl, 1, energy))
     $     call errquit('task_energy: failed storing energy',0)
c
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - wall
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, status))
     $     call errquit('task_energy: failed to store status',0)
      if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $     call errquit('task_energy: failed storing cputime',0)
      if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $     call errquit('task_energy: failed storing walltime',0)
c
      call ecce_print1('cpu time', mt_dbl, cpu, 1)
      call ecce_print1('wall time', mt_dbl, wall, 1)
      if (status) then
         if (theory .ne. 'md') call geom_print_rtdb_ecce(rtdb)
         call ecce_print1('total energy', mt_dbl, energy, 1)
         call ecce_print_module_exit('task energy', 'ok')
      else
         call ecce_print_module_exit('task energy', 'failed')
      endif
c         
      task_energy = status
c
      end
      logical function task_energy_doit(rtdb,theory,energy)
*
* $Id: task_energy.F,v 1.25 2001-02-20 23:31:29 edo Exp $
*
      implicit none
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
      logical  scf, rimp2, oimp2, direct_mp2, ccsd, mcscf, selci
      external scf, rimp2, oimp2, direct_mp2, ccsd, mcscf, selci
      logical semi_dir_mp2, nwargos, rimp2g,nwmd,pspw_energy
      external semi_dir_mp2, nwargos, rimp2g,nwmd,pspw_energy
      logical dft_energy, sodft_energy
      external dft_energy, sodft_energy
c
      logical status
      double precision dipole(3)
      character*80 key, prefix
c
      call ecce_print_module_entry('task energy')
c
      task_energy_doit = .false.
c
c     Set the theory in the database so that all is consistent for
c     the lower level modules.  Should not really have to do this,
c     but this value is used (appropriately and inappropriately) and
c     things are not yet disentangled.
c
      if (.not. rtdb_cput(rtdb, 'task:theory', 1, theory))
     $     call errquit('task_energy: setting theory?',0)
c
c     Actually do the deed
c
      prefix = theory           ! Most common scenario
      if (theory .eq. 'scf') then
         status = scf(rtdb)
      else if (theory .eq. 'dft') then
         status = dft_energy(rtdb)
      else if (theory .eq. 'sodft')then 
         status = sodft_energy(rtdb)
      else if (theory .eq. 'rimp2') then
         status = rimp2(rtdb)
      else if (theory .eq. 'rimp2_grad') then
         status = rimp2g(rtdb)
      else if (theory .eq. 'oimp2') then
         status = oimp2(rtdb)
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
     $        call errquit('task_energy: rtdb_put failed',0)
         status = semi_dir_mp2(rtdb)
      else if (theory .eq. 'cg_pspw') then
         status = pspw_energy(rtdb)
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
      else if (theory .eq. 'nwargos') then
         status = nwargos(rtdb)
      else
         call errquit('task_energy: unknown theory',0)
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
     $        call errquit('task_energy: no energy in rtdb?',0)
c
         key = ' '
         write(key,'(a,'':dipole'')') prefix(1:inp_strlen(prefix))
         if (rtdb_get(rtdb, key, mt_dbl, 3, dipole)) then
           if (.not. rtdb_put(rtdb, 'task:dipole', mt_dbl, 3, dipole))
     $          call errquit('task_energy: failed storing dipole',0)
	 endif
      endif
c
      task_energy_doit = status
c
      end
      subroutine rtdb_clean(rtdb, prefix)
      implicit none
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
     $           ('rtdb_clean: failed deleting existing entry',0)
         endif
         status = rtdb_next(rtdb, name)
         goto 10
      endif
c     
      end

