      logical function rohf(rtdb, energy)
C$Id: rohf.F,v 1.9 1996-05-14 18:00:05 d3g681 Exp $
      implicit none
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
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
      double precision energy   ! [output] Total SCF energy
c
c     local variables
c
      logical oprint, oprint_sym
      logical converged
      integer nprint
      double precision eone     ! One electron energy
      double precision etwo     ! Two electron energy
      double precision enrep    ! Effective nuclear repulsion energy
      double precision rohf_time
C
      Integer NAtoms, NShells
c
      logical rohf_nr_solve
      logical rohf_diis_solve
      external ga_create_atom_blocked, rohf_nr_solve
      external rohf_diis_solve
c
      oprint = util_print('information',print_low)
c
c     Initialize the ROHF data structures and allocate memory
c
      call rohf_init(rtdb)
c
c     Start timer
c
      call ga_sync()
      rohf_time = -util_cpusec()
c
c     Newton-Raphson based SCF
c
      if (odiis) then
        converged = rohf_diis_solve( rtdb, energy, eone, etwo, enrep )
      else
        converged = rohf_nr_solve(rtdb, energy,  eone, etwo, enrep)
      endif
c
      call ga_sync()
      rohf_time = rohf_time + util_cpusec()
c
c     Print out results
c
      if (ga_nodeid().eq.0 .and. (oprint .or. .not.converged)) then
         if (.not. converged) then
            write(6,*)
            call util_print_centered(6,
     $           'Calculation failed to converge', 20, .true.)
            write(6,*)
         endif
         write(6,2) scftype, energy, eone, etwo, enrep, rohf_time
 2       format(//
     $        '       Final ',a4,' results '/
     $        '       ------------------ '//
     $        '         Total SCF energy =', f20.12/
     $        '      One electron energy =', f20.12/
     $        '      Two electron energy =', f20.12/
     $        ' Nuclear repulsion energy =', f20.12//
     $        '        Time for solution =', f9.1,'s'//)
         call util_flush(6)
      endif
C
      if (util_print('schwarz',print_high)) then
         if ( .NOT. geom_ncent(geom, natoms) ) Call ErrQuit(
     $        'rohf: problem with call to geom_ncent', geom)
         if ( .NOT. bas_numcont(basis, nshells) ) Call ErrQuit(
     $        'rohf: problem with call to bas_numcont', basis)
         call schwarz_print(natoms, nshells)
      endif
C
      if (ga_nodeid() .eq. 0) then
         if (util_print('final evals', print_default)) then
            call util_print_centered(6,'Final eigenvalues',20,.true.)
            call output(dbl_mb(k_eval),
     $           1, min(nclosed+nopen+10,nmo),
     $           1, 1, nmo, 1, 1)
            call util_flush(6)
         endif
C
         if (util_print('final vectors', print_debug)) then
            write(6,*)
            call util_print_centered(6,'Final MO vectors',40,.true.)
            write(6,*)
            call util_flush(6)
         endif
      endif
c
      if (oadapt) then
         nprint = min(nclosed+nopen+10,nmo)
         if (util_print('all vector symmetries', print_high))
     $        nprint = nmo
         oprint_sym = util_print('final vector symmetries',
     $        print_default)
         call scf_movecs_sym_adapt(basis, g_movecs, oprint_sym,
     $        nprint, '- final')
      endif
c
      if (util_print('final vectors', print_debug)) then
         call ga_print(g_movecs)
      endif
c
      if (util_print('mulliken', print_high)) then
         call rohf_mulliken()
      endif
c
c     Clean up common blocks and allocated data and return
c
      call rohf_tidy
c
      rohf = converged
c
      end
      subroutine rohf_mulliken
      implicit none
#include "cscf.fh"
#include "global.fh"
c
      integer g_cdens, g_odens, g_over
c
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
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
      write(6,*)
      call util_print_centered(6,
     $     'Mulliken analysis of the total density', 20,.true.)
      write(6,*)
      call mull_pop(geom, basis, g_cdens, g_over)
c      
c     Do the open shell(s)
c
      if (nopen .gt. 0) then
         write(6,*)
         call util_print_centered(6,
     $        'Mulliken analysis of the open-shell density', 20,.true.)
         write(6,*)
         call mull_pop(geom, basis, g_odens, g_over)
         if (.not. ga_destroy(g_odens)) call errquit('rohf_mull: ga?',0)
      endif
c
      if (.not. ga_destroy(g_over)) call errquit('rohf_mull: ga?',0)
      if (.not. ga_destroy(g_cdens)) call errquit('rohf_mull: ga?',0)
c
      end
