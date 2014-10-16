      logical function rohf(rtdb, energy)
C$Id$
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"
#include "tcgmsg.fh"
#include "crohf.fh"
#include "cscf.fh"
#include "geom.fh"
#include "bas.fh"
c
c     Compute ROHF wavefunction.
c     Initial vectors etc. are assumed to be loaded in common/cscf/
c     by the driver routines
c
      integer rtdb
      double precision energy !< [Output] Total SCF energy
      double precision ecosmo !< The COSMO solvation energy
c
c     local variables
c
      logical oprint, oprint_sym
      logical converged
      logical odbug
      integer nprint
      double precision eone     ! One electron energy
      double precision etwo     ! Two electron energy
      double precision enrep    ! Effective nuclear repulsion energy
      double precision rohf_time
c
      integer ilo, ihi          ! For printing movecs analysis
      double precision eval_pr_tol_lo, eval_pr_tol_hi
      parameter (eval_pr_tol_lo = -1.5d0, eval_pr_tol_hi=0.5)
C
*      Integer NAtoms, NShells
c
      logical rohf_nr_solve
      logical rohf_diis_solve
      external ga_create_atom_blocked, rohf_nr_solve
      external rohf_diis_solve
c
      oprint = util_print('information',print_low)
      odbug  = util_print('rohf debug',print_high)
c
c     Initialize the ROHF data structures and allocate memory
c
      call rohf_init(rtdb)
c
      call ga_sync()
c
      if(odbug) then
         write(luout,*) 'ready to solve ...',ga_nodeid()
      endif
c
c     Start timer
c
      rohf_time = -util_cpusec()
c
c     Newton-Raphson based SCF
c
      ecosmo = 0.0d0
      if (maxiter .le. 0) then
         converged = .false.
         energy = 0.0d0
         eone = 0.0d0
         etwo = 0.0d0
      else if (odiis) then
         converged = rohf_diis_solve(rtdb, energy, eone, etwo, ecosmo,
     $                               enrep)
      else
         converged = rohf_nr_solve(rtdb, energy, eone, etwo, ecosmo,
     $                             enrep)
      endif
c
      call ga_sync()
      rohf_time = rohf_time + util_cpusec()
c
c     Print out results
c
      if (ga_nodeid().eq.0 .and. (oprint .or. .not.converged)) then
         if (.not. converged) then
            write(luout,*)
            call util_print_centered(luout,
     $           'Calculation failed to converge', 20, .true.)
            write(luout,*)
         endif
         write(luout,2) scftype, energy, eone, etwo, enrep
         if (abs(ecosmo).gt.0.0d0) then
           write(luout,3) ecosmo
         endif
         write(luout,4) rohf_time
 2       format(//
     $        '       Final ',a4,' results '/
     $        '       ------------------ '//
     $        '         Total SCF energy =', f20.12/
     $        '      One-electron energy =', f20.12/
     $        '      Two-electron energy =', f20.12/
     $        ' Nuclear repulsion energy =', f20.12)
 3       format(
     $        '             COSMO energy =', f20.12)
 4       format(/
     $        '        Time for solution =', f9.1,'s'//)
         call util_flush(6)
      endif
c
      call ecce_print1('total energy', mt_dbl, energy, 1)
      call ecce_print1('one-electron energy', mt_dbl, eone, 1)
      call ecce_print1('two-electron energy', mt_dbl, etwo, 1)
      call ecce_print1('nuclear repulsion energy', mt_dbl, enrep, 1)
c
c     These statistics are currently junk and need to be reworked
c
c$$$      if (util_print('schwarz',print_high)) then
c$$$         if ( .NOT. geom_ncent(geom, natoms) ) Call ErrQuit(
c$$$     $        'rohf: problem with call to geom_ncent', geom)
c$$$         if ( .NOT. bas_numcont(basis, nshells) ) Call ErrQuit(
c$$$     $        'rohf: problem with call to bas_numcont', basis)
c$$$         call schwarz_print(natoms, nshells)
c$$$      endif
C
      if (oadapt) then
         nprint = min(nclosed+nopen+10,nmo)
         if (util_print('all vector symmetries', print_high))
     $        nprint = nmo
         oprint_sym = util_print('final vector symmetries',
     $        print_default)
         call scf_movecs_sym_adapt(basis, g_movecs, oprint_sym,
     $        nprint, '- final', .true., int_mb(k_irs))
         call ecce_print1('orbital symmetries', mt_int,
     $        int_mb(k_irs), nmo)
      endif
c
      if (ga_nodeid() .eq. 0) then
         if (util_print('final evals', print_default)) then
            nprint = min(nclosed+nopen+10,nmo)
            if (util_print('final evals', print_high)) nprint = nmo
            call util_print_centered(luout,'Final eigenvalues',20,
     $                               .true.)
            call output(dbl_mb(k_eval),
     $           1, nprint,
     $           1, 1, nmo, 1, 1)
*            call ecce_print1('orbital energies', mt_dbl, 
*     $           dbl_mb(k_eval), nmo)
            call util_flush(6)
         endif
      endif
c
      call movecs_fix_phase(g_movecs)
      if (util_print('final vectors analysis', print_default)) then
         do ilo = 1,max(1,nclosed-10)
            if (dbl_mb(k_eval+ilo-1) .ge. eval_pr_tol_lo) goto 961
         enddo
 961     do ihi = min(nclosed+nopen+10,nmo),nmo
            if (dbl_mb(k_eval+ihi-1) .ge. eval_pr_tol_hi) goto 9611
         enddo
         ihi = max(ihi-1,1)
 9611    continue
         if (util_print('final vectors analysis', print_high)) then
            ilo = 1
            ihi = nmo
         endif
         call movecs_print_anal(basis, ilo, ihi, 0.15d0, g_movecs, 
     $        'ROHF Final Molecular Orbital Analysis', 
     $        .true., dbl_mb(k_eval), oadapt, int_mb(k_irs),
     $        .true., dbl_mb(k_occ))
      endif
c
      if (ga_nodeid() .eq. 0) then
         if (util_print('final vectors', print_debug)) then
            write(luout,*)
            call util_print_centered(luout,'Final MO vectors',40,
     $                               .true.)
            write(luout,*)
            call util_flush(6)
         endif
      endif
c
      if (util_print('final vectors', print_debug)) 
     $     call ga_print(g_movecs)
c
      call rohf_analyze(rtdb)
c
c     Clean up common blocks and allocated data and return
c
      call rohf_tidy
c
      rohf = converged
c
      end
      subroutine rohf_analyze(rtdb)
      implicit none
#include "errquit.fh"
#include "cscf.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "bas.fh"
      integer rtdb
c
c     Analyze the ROHF wavefunction.  Optionally print the mulliken
c     analysis and multipoles, always store the dipole in the RTDB
c     
      integer g_cdens, g_odens, g_over
c     
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      logical oprintmulliken, oprintmultipole,oprintinertia
c
      integer len_moments, lmax
      parameter (lmax = 2)
      parameter (len_moments = (lmax+1)*(lmax+2)*(lmax+3)/6)
      double precision totalmoments(len_moments) 
      double precision openmoments (len_moments)
      double precision nuclmoments (len_moments)
      double precision center(3)
      integer k, l, m, ltotal, ind
      data center/3*0.0d0/
c
      oprintmulliken  = util_print('mulliken', print_default)
      oprintmultipole = util_print('multipole', print_default)
      oprintinertia = util_print('inertia', print_default).and.
     . ga_nodeid().eq.0                                    
c
c     moments of inertia
c
      if(oprintinertia) call geom_momint(geom)
c     
c     First analyze the total density
c     
      g_cdens = ga_create_atom_blocked(geom, basis, 'rohf_mull:cdens')
      g_over  = ga_create_atom_blocked(geom, basis, 'rohf_mull:over')
      call ga_zero(g_cdens)
      call ga_zero(g_over)
      call int_1e_ga(basis, basis, g_over, 'overlap', .false.)
c     
      if (nopen .gt. 0) then
         g_odens = ga_create_atom_blocked(geom, basis,'rohf_mull:odens')
         call ga_zero(g_odens)
      else
         g_odens = 0            ! Invalid handle
      endif
      call rohf_dens(g_movecs, nbf, nclosed, nopen, g_cdens, g_odens)
      if (nopen .gt. 0)
     $     call ga_dadd(1.0d0, g_cdens, 1.0d0, g_odens, g_cdens)
c     
      if (oprintmulliken) then
         if (ga_nodeid() .eq. 0) then
            write(luout,*)
            call util_print_centered(luout,
     $           'Mulliken analysis of the total density', 20,.true.)
         endif
         call mull_pop(geom, basis, g_cdens, g_over, 'total')
c     
         if (nopen .gt. 0) then
            if (ga_nodeid() .eq. 0) then
               write(luout,*)
               call util_print_centered(luout,
     $              'Mulliken analysis of the open-shell density', 
     $              20,.true.)
            endif
            call mull_pop(geom, basis, g_odens, g_over, 'open shell')
         endif
      endif
c     
      if (.not. bas_cando_mpoles(basis)) goto 100
c     
      call dfill(len_moments, 0.0d0, totalmoments, 1)
      call dfill(len_moments, 0.0d0, nuclmoments, 1)
      call dfill(len_moments, 0.0d0, openmoments, 1)
      call geom_nuc_mpole(geom, center, lmax, nuclmoments, len_moments)
      call multipole_density(basis, center, lmax, g_cdens, 
     $     totalmoments, len_moments)
      if (nopen .gt. 0) 
     $     call multipole_density(basis, center, lmax, g_odens, 
     $     openmoments, len_moments)
      do k = 1, len_moments
         totalmoments(k) = totalmoments(k) + nuclmoments(k)
      enddo
c     
      if (oprintmultipole .and. ga_nodeid().eq.0) then
         write(luout,*)
         call util_print_centered(luout,
     $     'Multipole analysis of the density wrt the origin', 
     $      30, .true.)
         write(luout,*)
         write(luout,*) '    L   x y z        total         open',
     $                 '         nuclear'
         write(luout,*) '    -   - - -        -----         ----',
     $        '         -------'
         ind = 0
         do ltotal = 0, lmax
            do k = ltotal, 0, -1
               do l = ltotal-k, 0, -1
                  m = ltotal - k - l
                  ind = ind + 1
                  write(luout,12) ltotal, k, l, m, totalmoments(ind),
     $                 openmoments(ind), nuclmoments(ind)
 12               format(4x,i2,2x,3i2,3f14.6)
               enddo
            enddo
            write(luout,*)
         enddo
      endif
      call ecce_print1('total dipole', mt_dbl, totalmoments(2), 3)
      call ecce_print1('open shell dipole', mt_dbl, openmoments(2), 3)
      call ecce_print1('nuclear dipole', mt_dbl, nuclmoments(2), 3)
      call ecce_print1('total quadrupole', mt_dbl, totalmoments(5), 6)
      call ecce_print1('open shell quadrupole',mt_dbl,openmoments(5),6)
      call ecce_print1('nuclear quadrupole', mt_dbl, nuclmoments(5), 6)
c
c     Store the dipole moment
c
      if (.not. rtdb_put(rtdb, 'scf:dipole', mt_dbl, 
     $     3, totalmoments(2))) call errquit('multipole: rtdb?',0,
     &       RTDB_ERR)
c     
 100  if (nopen .gt. 0) then
         if (.not. ga_destroy(g_odens)) 
     $        call errquit('rohf_mull: ga?',0, GA_ERR)
      endif
      if (.not. ga_destroy(g_over)) call errquit('rohf_mull: ga?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_cdens)) call errquit('rohf_mull: ga?',0,
     &       GA_ERR)
c     
      end
