      subroutine movecs_fix_phase(g_vecs)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer g_vecs
c
      integer l_x, k_x, n, m, i, j, k, type
      double precision v
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
         call ga_get(g_vecs, 1, n, i, i, dbl_mb(k_x), n)
         v = 0.0d0
         k = 0
         do j = 0, n-1
c
c          The next if-statement deals with coefficients that should
c          be equally large, e.g. due to symmetry, but in practice
c          are not exactly equally large. Now a coefficient has to be
c          at least 10% larger than the previously largest one to be 
c          considered larger.
c
           if (abs(dbl_mb(k_x+j)).gt.(v*1.1d0)) then
             v = abs(dbl_mb(k_x+j))
             k = j
           endif
         enddo
         if (dbl_mb(k_x+k) .lt. 0.0d0) then
            call dscal(n, -1.0d0, dbl_mb(k_x), 1)
            call ga_put(g_vecs, 1, n, i, i, dbl_mb(k_x), n)
         endif
      enddo
c
      if (.not. ma_pop_stack(l_x)) call errquit
     $     ('movecs_fp: failed to pop stack', 0, MA_ERR)
c
      call ga_sync
c
      end

      
