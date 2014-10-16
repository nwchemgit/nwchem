      subroutine rohf_rotate_movecs(step, g_search)
*
* $Id$
*
      implicit none
#include "errquit.fh"
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
     $        ('rohf_rotate_movecs: ga create failed', nbf, GA_ERR)
         call ga_dscal(g_search, step)
         call rohf_k2cf(basis, nbf, nmo, nclosed,
     $        nopen, g_search, g_movecs, g_tmp)
         scale = 1.0d0/step
         call ga_dscal(g_search, scale)
         call ga_copy(g_tmp, g_movecs)
         if (.not. ga_destroy(g_tmp)) call errquit
     $        ('rohf_rotate_movecs: tmp?', 0, GA_ERR)
      endif
c     
      end
