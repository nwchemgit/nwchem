       subroutine rhf_dens_to_mo(rtdb, geom, basis, nelec, nbf, tol2e,
     $     g_dens, g_mocf, oskel)
C$Id: rhf_dens_mo.F,v 1.1 1995-12-13 01:25:36 d3g681 Exp $
       implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "tcgmsg.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "cscfps.fh"
#include "util.fh"
       integer rtdb, geom, basis, nbf, nelec
       double precision tol2e   ! [input] selection threshold for 2e integs
       integer g_dens           ! [input] AO density in GA
       integer g_mocf           ! [output] MO vectors in GA
c
c     Given a density in g_dens return an orthonormal set of
c     mo vectors in g_mocf by doing an approximate Fock build
c     and a diagonalization
c
       logical oprint, osplit, oskel
       integer g_fock, g_over, g_hcore
       integer nocc
       integer l_ev, k_ev
       integer i,j
       double precision escf, eone, etwo, enrep, dens_norm
       integer ga_create_atom_blocked
       external ga_create_atom_blocked
c
      oprint = util_print('mo guess', print_default)
c
c     The guess density must be symmetrized if sym is used
c
      if (oskel) then
         if (oscfps) call pstat_on(ps_sym_sym)
         call sym_symmetrize(geom, basis, .true., g_dens)
         if (oscfps) call pstat_off(ps_sym_sym)
      endif
c
       if (.not. geom_nuc_rep_energy(geom, enrep))
     $      call errquit('rhf_dens_to_mo: no enrep?', 0)
       if (.not. ma_push_get(MT_DBL, nbf, 'Evals', l_ev, k_ev))
     $      call errquit('rhf_dens_to_mo: push of evals failed', nbf)
c
       g_fock  = ga_create_atom_blocked(geom, basis,'Temp Fock')
       g_hcore = ga_create_atom_blocked(geom, basis,'Temp HCore')
       g_over  = ga_create_atom_blocked(geom, basis,'Temp Over')
c
c     If the norm of the density matrix is non-zero then adjust
c     it so that it specifies the correct no. of electrons
c
       call ga_zero(g_over)
       call int_1e_ga(basis, basis, g_over, 'overlap', .false.)
       if (util_print('ao overlap',print_debug))
     $      call ga_print(g_over)
       dens_norm = ga_ddot(g_over, g_dens)
       if (abs(dens_norm) .gt. 1.0d-2 .and.
     $      abs(dens_norm-dble(nelec)).gt.1.0d-6) then
          if (ga_nodeid() .eq. 0 .and. oprint) then
             write(6,17) dens_norm, nelec
 17          format(/' Renormalizing density from ',f8.2,' to ',i4)
             call util_flush(6)
          endif
          dens_norm = dble(nelec)/dens_norm
          call ga_dscal(g_dens, dens_norm)
       endif
c
c     Diagonalize the overlap to examine eigenvalues for linear dependence
c
       if (oprint) then
#ifdef PARALLEL_DIAG
          call ga_diag_std(g_over, g_mocf, dbl_mb(k_ev))
#else
          call ga_diag_std_seq(g_over, g_mocf, dbl_mb(k_ev))
#endif
          if (ga_nodeid() .eq. 0) then
             do i = 0, nbf
                if (dbl_mb(k_ev+i).gt.1d-4) goto 333
             enddo
 333         if (i .gt. 0) then
                write(6,334) i
 334            format(/' The overlap matrix has ',i3,' eigenvalues ',
     $               'less than 10^-4'/)
                write(6,33) (dbl_mb(k_ev+j),j=0,i-1)
 33             format(1p,8d9.2)
                write(6,*)
                call util_flush(6)
             endif
          endif
       endif
c
c     Build the fock_matrix and overlap matrix
c
       call ga_zero(g_fock)
       call ga_zero(g_hcore)
       call ga_zero(g_over)
       call int_1e_ga(basis, basis, g_over, 'overlap', .false.)
       call int_1e_ga(basis, basis, g_hcore,'kinetic', oskel)
       call int_1e_ga(basis, basis, g_hcore,'potential', oskel)
c
       if (abs(dens_norm) .gt. 1d-2)
     $      call rhf_fock_2e(geom, basis, g_dens, g_fock, tol2e,
     &        .true., .true., oskel)
c
c     Compute contributions to the energy and diagonalize ... note
c     that we can do this with the skelton matrices
c
       call ga_screen(g_hcore, 1.0d-13) ! Force zeroes due to symmetries
       call ga_screen(g_fock, 1.0d-13) ! Force zeroes due to symmetries
       eone = ga_ddot(g_hcore,g_dens)
       etwo = 0.5d0*ga_ddot(g_fock,g_dens)
       escf = eone + etwo + enrep
       call ga_dadd(1.d0,g_hcore,1.d0,g_fock,g_fock)
c
       if (oskel) then
         if (oscfps) call pstat_on(ps_sym_sym)
          call sym_symmetrize(geom, basis, .false., g_fock)
         if (oscfps) call pstat_off(ps_sym_sym)
       endif
c
c$$$c
c$$$c     Debug for george
c$$$c
c$$$       do i = 1, nbf
c$$$          do j = 1, i
c$$$             call ga_get(g_fock, i, i, j, j, ff, 1)
c$$$             call ga_get(g_over, i, i, j, j, oo, 1)
c$$$             write(6,111) i, j, ff, oo
c$$$ 111         format(2i5, 2d22.14)
c$$$          enddo
c$$$       enddo
c$$$       call util_flush(6)
c
       if (oscfps) call pstat_on(ps_diag)
#if defined(PARALLEL_DIAG)
       call ga_diag(g_fock,g_over,g_mocf,dbl_mb(k_ev))
#else
       call ga_diag_seq(g_fock,g_over,g_mocf,dbl_mb(k_ev))
#endif
       if (oscfps) call pstat_off(ps_diag)
       call ga_sync()
c
       nocc = nelec/2           !
**       osplit = abs(dbl_mb(k_ev+nocc-1)-dbl_mb(k_ev+nocc)).lt.1.0e-3
       osplit = .false.
       if (ga_nodeid().eq.0 .and. (oprint .or. osplit)) then
          write(6,*)
          call util_print_centered(6, 'Non-variational initial energy',
     $         20, .true.)
          write(6,901) escf, eone, etwo,
     $         dbl_mb(k_ev+nocc-1), dbl_mb(k_ev+nocc)
 901     format(/
     $         ' Total energy = ',f14.6/
     $         ' 1-e energy   = ',f14.6/
     $         ' 2-e energy   = ',f14.6/
     $         ' HOMO         = ',f14.6/
     $         ' LUMO         = ',f14.6/)
c
         if (osplit) then
            write(6,902)
 902        format(' !! WARNING: occupation is splitting degenerate',
     $           ' orbitals'/)
            if (oskel) write(6,903)
 903        format(' !! Cannot use symmetry with ',
     $           'non-symmetric occupation ... swap orbitals or '/
     $           ' !! must disable skeleton Fock build'/)
            write(6,*) ' Eigenvalues '
            write(6,*)
            call output(dbl_mb(k_ev), 1, nbf, 1, 1, nbf, 1, 1)
         endif
         call util_flush(6)
       endif
c
       if (.not. ga_destroy(g_fock)) call errquit
     $      ('rhf_dens_to_mo: destroy of fock?', 0)
       if (.not. ga_destroy(g_over)) call errquit
     $      ('rhf_dens_to_mo: destroy of over?', 0)
       if (.not. ga_destroy(g_hcore)) call errquit
     $      ('rhf_dens_to_mo: destroy of hcore?', 0)
       if (.not. ma_pop_stack(l_ev)) call errquit
     $      ('rhf_dens_to_mo: destroy of ev?', 0)
c
       end
