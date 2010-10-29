      subroutine rhf_dens_to_mo(rtdb, geom, basis, nelec, nbf, nmo,
     $     tol2e, g_dens, g_mocf, oskel, evals )
C     $Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "tcgmsg.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "cscfps.fh"
#include "util.fh"
#include "stdio.fh"
      integer rtdb, geom, basis, nbf, nmo, nelec
      double precision tol2e    ! [input] selection threshold for 2e integs
      integer g_dens            ! [input] AO density in GA
      integer g_mocf            ! [input/output] Indep vectors/MO vectors in GA
      double precision evals(nbf)  ! [output] Orbital energies
c     
c     Given a density in g_dens and the nmo independent eigenvectors of S 
c     in g_mocf, return an orthonormal set of
c     mo vectors in g_mocf by doing an approximate Fock build
c     and a diagonalization
c     
      logical oprint, oskel
      integer g_fock, g_over, g_hcore, g_tmp, g_mofock, g_u
      integer nocc
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
     $     call errquit('rhf_dens_to_mo: no enrep?', 0,
     &       GEOM_ERR)
c     
      g_fock  = ga_create_atom_blocked(geom, basis,'Temp Fock')
      g_hcore = ga_create_atom_blocked(geom, basis,'Temp HCore')
      g_over  = g_hcore
c     
c     If the norm of the density matrix is non-zero then adjust
c     it so that it specifies the correct no. of electrons
c     
      if (util_print('ao kinetic',print_debug)) then
        call ga_zero(g_over)
        call int_1e_ga(basis, basis, g_over, 'kinetic', .false.)
        call ga_print(g_over)
      endif
c
      call ga_zero(g_over)
      call int_1e_ga(basis, basis, g_over, 'overlap', .false.)
      if (util_print('ao overlap',print_debug))
     $     call ga_print(g_over)
      dens_norm = ga_ddot(g_over, g_dens)
      if (abs(dens_norm) .gt. 1.0d-2 .and.
     $     abs(dens_norm-dble(nelec)).gt.1.0d-6) then
         if (ga_nodeid() .eq. 0 .and. oprint) then
            write(luout,17) dens_norm, nelec
 17         format(/' Renormalizing density from ',f10.2,' to ',i6)
            call util_flush(luout)
         endif
         dens_norm = dble(nelec)/dens_norm
         call ga_dscal(g_dens, dens_norm)
      endif
c     
c     (Ortho)normalize the eigenvectors of S while we still have S
c     
      call ga_orthog(g_mocf, g_over, .true.)
c     
c     Build the fock_matrix
c     
      call ga_zero(g_fock)
      call ga_zero(g_hcore)
      call int_1e_ga(basis, basis, g_hcore,'kinetic', oskel)
      call int_1e_ga(basis, basis, g_hcore,'potential', oskel)
c     
      if (nelec.gt.1 .and. abs(dens_norm) .gt. 1d-2)
     $     call rhf_fock_2e(geom, basis, g_dens, g_fock, tol2e,
     &     .true., .true., oskel)
c     
c     Compute contributions to the energy and symmetrize
c     
*      call ga_screen(g_hcore, 1.0d-12) ! Force zeroes
*      call ga_screen(g_fock, 1.0d-12) ! Force zeroes
      eone = ga_ddot(g_hcore,g_dens)
      etwo = 0.5d0*ga_ddot(g_fock,g_dens)
      escf = eone + etwo + enrep
      call ga_dadd(1.d0,g_hcore,1.d0,g_fock,g_fock)
      if (.not. ga_destroy(g_hcore)) call errquit
     $     ('rhf_dens_to_mo: destroy of hcore?', 0, GA_ERR)
c     
      if (oskel) then
         if (oscfps) call pstat_on(ps_sym_sym)
         call sym_symmetrize(geom, basis, .false., g_fock)
         if (oscfps) call pstat_off(ps_sym_sym)
      endif
c     
c     Transform the fock matrix to the orthonormal subspace and diagonalize
c     
      if (.not. ga_create(MT_DBL, nbf, nmo, 'rhf_d_m: tmp',
     $     32, 32, g_tmp)) call errquit('rhf_d_m: tmp',0, GA_ERR)
      if (.not. ga_create(MT_DBL, nmo, nmo, 'rhf_d_m: mof',
     $     32, 32, g_mofock)) call errquit('rhf_d_m: mof',0, GA_ERR)
      call two_index_transf(g_fock, g_mocf, g_mocf, g_tmp, g_mofock )
      if (.not. ga_destroy(g_fock)) call errquit('rhf_d_m: ga?',0,
     &       GA_ERR)
c     
      call ga_copy(g_mocf, g_tmp) ! Save the orthonormal basis
      if (.not. ga_create(MT_DBL, nmo, nmo, 'rhf_d_m: u',
     $     32, 32, g_u)) call errquit('rhf_d_m: u',0, GA_ERR)
*      write(luout,*) ' g_mofock ', g_mofock
*      write(luout,*) ' g_u ', g_u
*      call ga_summarize(0)
c     
      if (oscfps) call pstat_on(ps_diag)
#ifdef PARALLEL_DIAG
#ifdef SCALAPACK
      call ga_pdsyev(g_mofock, g_u, evals, 0)
#else
      call ga_diag_std(g_mofock, g_u, evals )
#endif
#else
      call ga_diag_std_seq(g_mofock, g_u, evals )
#endif
      if (oscfps) call pstat_off(ps_diag)
      call ga_sync()
c     
      call ga_dgemm('n', 'n', nbf, nmo, nmo, 1.0d0, g_tmp, g_u,
     $     0.0d0, g_mocf)
      if (.not. ga_destroy(g_u)) call errquit('rhf_d_m: ga?',0, GA_ERR)
      if (.not. ga_destroy(g_mofock)) call errquit('rhf_d_m: ga?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_tmp)) call errquit('rhf_d_m: ga?',0,
     &       GA_ERR)
c     
      nocc = (nelec+1)/2
      if (ga_nodeid().eq.0 .and. oprint) then
         write(luout,*)
         call util_print_centered(luout,
     $        'Non-variational initial energy',
     $        20, .true.)
         write(luout,901) escf, eone, etwo, evals(nocc)
         if (nmo.ne.nocc)  write(luout,9011) evals(nocc+1)
         write(luout,*)
 901     format(/
     $        ' Total energy = ',f14.6/
     $        ' 1-e energy   = ',f14.6/
     $        ' 2-e energy   = ',f14.6/
     $        ' HOMO         = ',f14.6)
 9011    format(' LUMO         = ',f14.6)
c
         if ( oskel .and. nmo.ne.nocc .and. 
     $        (abs(evals(nocc)-evals(nocc+1)) .lt. 1e-3)) then
            write(luout,902) 
 902        format(/,1x,70('!'),/,
     $           ' Warning - the HOMO and LUMO are degenerate and you',
     $           ' are using symmetry.',/,
     $           ' This can lead to non-variational energies and poor',
     $           ' convergence.',/,
     $           ' Modify the initial guess, or use an open-shell',
     $           ' wavefunction, or turn',/,
     $           ' off symmetry.',/,1x,70('!'),/)
         endif
c     
         call util_flush(luout)
      endif
c     
      end
