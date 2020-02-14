      subroutine civecs_fix_phase(g_x,g_y,ipol,tda)
*
* $Id: movecs_phase.F 22758 2012-08-22 21:39:19Z d3y133 $
*
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer ipol
      integer g_x(2)
      integer g_y(2)
      logical tda
c
      integer l_x, k_x, n, m, i, j, k, type
      integer l_y, k_y
      integer l_x2, k_x2, l_y2, k_y2
      integer n2, m2
      double precision v
      logical scale_flag
      integer idamax
      external idamax
c
c     g_x is a matrix (n, m) regarded as m vectors of length n
c     ... modify the phases so that the largest component of
c     each vector is positive, g_y will be scaled by the same factor
c
      call ga_sync
c
      call ga_inquire(g_x(1),type,n,m)
      if (.not. ma_push_get(MT_DBL,n,'movecs_fp',l_x,k_x))
     $     call errquit('civecs_fp: failed to allocate ',n,MA_ERR)
      if (.not.tda) then
        if (.not. ma_push_get(MT_DBL,n,'movecs_fp',l_y,k_y))
     $       call errquit('civecs_fp: failed to allocate ',n,MA_ERR)
      end if
      if (ipol.gt.1) then
        call ga_inquire(g_x(2),type,n2,m2)
        if (.not. ma_push_get(MT_DBL,n2,'movecs_fp',l_x2,k_x2))
     $       call errquit('civecs_fp: failed to allocate ',n2,MA_ERR)
        if (.not.tda) then
          if (.not. ma_push_get(MT_DBL,n2,'movecs_fp',l_y2,k_y2))
     $         call errquit('civecs_fp: failed to allocate ',n2,MA_ERR)
        end if
      end if
c
      do i = ga_nodeid()+1,m, ga_nnodes()
         call ga_get(g_x(1),1,n,i,i,dbl_mb(k_x),n)
         if (.not.tda) call ga_get(g_y(1),1,n,i,i,dbl_mb(k_y),n)
         if (ipol.gt.1) then
           call ga_get(g_x(2),1,n2,i,i,dbl_mb(k_x2),n2)
           if (.not.tda) call ga_get(g_y(2),1,n2,i,i,dbl_mb(k_y2),n2)
         end if
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
         if (dbl_mb(k_x+k).lt.0.0d0) then
            call dscal(n,-1.0d0,dbl_mb(k_x),1)
            call ga_put(g_x(1),1,n,i,i,dbl_mb(k_x),n)
            if (.not.tda) then
              call dscal(n,-1.0d0,dbl_mb(k_y),1)
              call ga_put(g_y(1),1,n,i,i,dbl_mb(k_y),n)
            end if
            if (ipol.gt.1) then
              call dscal(n2,-1.0d0,dbl_mb(k_x2),1)
              call ga_put(g_x(2),1,n2,i,i,dbl_mb(k_x2),n2)
              if (.not.tda) then
                call dscal(n2,-1.0d0,dbl_mb(k_y2),1)
                call ga_put(g_y(2),1,n2,i,i,dbl_mb(k_y2),n2)
              end if
            end if
         endif
      enddo
c
      if (.not. ma_chop_stack(l_x)) call errquit
     $     ('civecs_fp: failed to chop stack', 0, MA_ERR)
c
      call ga_sync
c
      end

      
