      subroutine dft_cg_guess(rtdb)
      implicit none
c
c     The purpose of this routine is to get the guess vectors into
c     g_movecs so that we can start the quasi-Newton-Raphson solver.
c
#include "global.fh"
#include "cscf.fh"
c
      integer rtdb
      integer g_dens(2)
      logical status
c
      status = ga_duplicate(g_movecs(1),g_dens(1),'tmp dens')
      if (scftype.eq.'UHF') then
        status = ga_duplicate(g_movecs(2),g_dens(2),'tmp dens')
      endif
      call scf_vectors_guess(rtdb, tol2e, geom, basis,
     &                       trans, movecs_in, movecs_out,
     &                       movecs_guess, scftype, nclosed, nopen,
     &                       nbf, nmo, nalpha, nbeta,  k_eval, k_occ,
     &                       k_irs, g_movecs, g_dens, 'vectors',
     &                       'dft', title, oskel, oadapt,
     &                       .true.)
      status = ga_destroy(g_dens(1))
      if (scftype.eq.'UHF') then
        status = ga_destroy(g_dens(2))
      endif
c
      end
c $Id$
