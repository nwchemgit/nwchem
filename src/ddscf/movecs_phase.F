      subroutine movecs_fix_phase(g_vecs)
*
* $Id: movecs_phase.F,v 1.3 2003-10-17 22:54:33 carlfahl Exp $
*
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer g_vecs
c
      integer l_x, k_x, n, m, i, j, type
      integer idamax
      external idamax
c
c     g_vecs is a matrix (n, m) regarded as m vectors of length n
c     ... modify the phases so that the largest component of
c     each vector is positive
c
      call ga_sync
c
      call ga_inquire(g_vecs, type, n, m)
      if (.not. ma_push_get(MT_DBL, n, 'movecs_fp', l_x, k_x))
     $     call errquit('movecs_fp: failed to allocate ', n, MA_ERR)
c
      do i = ga_nodeid()+1,m, ga_nnodes()
         call ga_get(g_vecs, 1, n, i, i, dbl_mb(k_x), 1)
         j = idamax(n, dbl_mb(k_x), 1)
         if (dbl_mb(k_x+j-1) .lt. 0.0d0) then
            call dscal(n, -1.0d0, dbl_mb(k_x), 1)
            call ga_put(g_vecs, 1, n, i, i, dbl_mb(k_x), 1)
         endif
      enddo
c
      if (.not. ma_pop_stack(l_x)) call errquit
     $     ('movecs_fp: failed to pop stack', 0, MA_ERR)
c
      call ga_sync
c
      end

      
