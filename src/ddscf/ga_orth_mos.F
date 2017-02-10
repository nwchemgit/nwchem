      subroutine ga_orthog_mos(basis, g_vecs)
C$Id$
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "cscfps.fh"
      integer basis
      integer g_vecs
c
c     Convenience routine that allocates temporary global memory,
c     makes the overlap matrix and orthogonalizes the MOs in
c     g_vecs
c
      integer g_over, nbf, geom
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c
      if (.not. bas_numbf(basis, nbf)) call errquit
     $  ('orthog mos: bas_numbf', 0, BASIS_ERR)
      if (.not. bas_geom(basis, geom)) call errquit
     $  ('orthog mos: bas_geom?', 0, BASIS_ERR)
c
*      if (.not. ga_create(MT_DBL, nbf, nbf, 'orthog_mos:overlap',
*     $    nbf, 0, g_over)) call errquit
*     $    ('orthog_mos: failed to create overlap matrix', 0)
      g_over = ga_create_atom_blocked(geom, basis, 
     $     'orthog_mos: overlap')
c
      call ga_zero(g_over)
      call int_1e_ga(basis, basis, g_over, 'overlap', .false.)
      if (oscfps) call pstat_on(ps_orthog)
      call ga_orthog(g_vecs, g_over, .true.)
      if (oscfps) call pstat_off(ps_orthog)
c
      if (.not. ga_destroy(g_over))
     $  call errquit('orthog mos: ga destroy ?', 0, GA_ERR)
c
      end
