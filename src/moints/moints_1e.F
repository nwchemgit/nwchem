      subroutine moints_1e( nbf, basis, g_movecs, g_onel )
*
* $Id$
*
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer nbf
      integer basis
      integer g_movecs
      integer g_onel
      logical oskel, status
      integer g_tmp
      data oskel/.false./
c
c     1-e component
c     
      status = ga_duplicate( g_onel, g_tmp, 'crap')
      call ga_zero(g_onel)
      call int_1e_ga(basis, basis, g_onel, 'kinetic', oskel)
      call int_1e_ga(basis, basis, g_onel, 'potential', oskel)
      call ga_dgemm('n', 'n', nbf, nbf, nbf, 1.d0, g_onel, g_movecs,
     $               0.d0, g_tmp )
      call ga_dgemm('t', 'n', nbf, nbf, nbf, 1.d0, g_movecs, g_tmp,
     $               0.d0, g_onel )

      status = ga_destroy(g_tmp)
c
c
c
      return
      end

      
      
