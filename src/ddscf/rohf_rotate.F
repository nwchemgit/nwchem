      subroutine rohf_rotate_movecs(step, g_search)
*
* $Id: rohf_rotate.F,v 1.4 1997-11-07 11:35:17 d3e129 Exp $
*
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "cscf.fh"
c
      double precision step
      integer g_search
c
c     Rotate the movecs (ga handle in common) by exp(step*search)
c
      integer g_tmp
      double precision scale
c
      if (step .ne. 0.0d0) then
         if (.not. ga_create(MT_DBL, nbf, nmo, 
     $        'rohf_rotate_movecs: tmp', 32, 32, g_tmp)) call errquit
     $        ('rohf_rotate_movecs: ga create failed', nbf)
         call ga_dscal(g_search, step)
         call rohf_k2cf(basis, nbf, nmo, nclosed,
     $        nopen, g_search, g_movecs, g_tmp)
         scale = 1.0d0/step
         call ga_dscal(g_search, scale)
         call ga_copy(g_tmp, g_movecs)
         if (.not. ga_destroy(g_tmp)) call errquit
     $        ('rohf_nr_solve: tmp?', 0)
      endif
c     
      end
