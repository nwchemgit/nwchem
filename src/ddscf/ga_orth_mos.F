      subroutine ga_orthog_mos(basis, g_vecs)
C$Id: ga_orth_mos.F,v 1.1 1995-12-13 01:25:19 d3g681 Exp $
      implicit none
#include "bas.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "cscfps.fh"
      integer basis
      integer g_vecs
c
c     Convenience routine that allocates temporary global memory,
c     makes the overlap matrix and orthogonalizes the MOs in
c     g_ves
c
      integer g_over, nbf
c
      if (.not. bas_numbf(basis, nbf)) call errquit
     $  ('orthog mos: bas_numbf ', 0)
c
      if (.not. ga_create(MT_DBL, nbf, nbf, 'orthog_mos:overlap',
     $    nbf, 1, g_over)) call errquit
     $    ('orthog_mos: failed to create overlap matrix', 0)
c
      call ga_zero(g_over)
      call int_1e_ga(basis, basis, g_over, 'overlap', .false.)
      if (oscfps) call pstat_on(ps_orthog)
      call ga_orthog(g_vecs, g_over, .true.)
      if (oscfps) call pstat_off(ps_orthog)
c
      if (.not. ga_destroy(g_over))
     $  call errquit('orthog mos: ga destroy ?', 0)
c
      end
