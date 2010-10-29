      subroutine ga_mix(g_a, n, nvec, b, ld)
C$Id$
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
      integer g_a
      integer n, nvec, ld
      double precision b(ld,nvec)
c
c     This routine is set up to optimize the rotation of
c     a (small) set of vectors amoung themselves.  
c
c     A must be distributed by columns so that an entire
c     row is present on a processor.
c
c     A(i,j) <- sum(l=1..k) A(i,l)*B(l,j)
c
      integer ilo, ihi, jlo, jhi
      integer i, j, idim, jdim
      integer l_c, k_c, l_a, k_a
c
c     Determine the distribution
c
      call ga_sync()
      call ga_distribution(g_a, ga_nodeid(), ilo, ihi, jlo, jhi)
      if (ihi.gt.0 .and. ihi.ge.ilo) then
c
c     This process holds some data ... check have an entire row
c     and then access
c
         idim = ihi - ilo + 1
         if (jhi .gt. nvec) jhi = nvec
         jdim = jhi - jlo + 1
         if (jlo.ne.1 .or. jhi.ne.nvec) call errquit
     $        ('ga_mix: must distribute by columns', 0, GA_ERR)
c
         if (.not.ma_push_get(MT_DBL, nvec, 'ga_mix:c', l_c, k_c))
     $        call errquit('ga_mix: ma_push_get? c', 0, MA_ERR)
         if (.not.ma_push_get(MT_DBL, nvec, 'ga_mix:a', l_a, k_a))
     $        call errquit('ga_mix: ma_push_get? a', 0, MA_ERR)
c
         do i = ilo, ihi
            call ga_get(g_a, i, i, 1, nvec, dbl_mb(k_a), 1)
            do j = 1, nvec
               dbl_mb(k_c+j-1) =
     $              ddot(nvec, dbl_mb(k_a), 1, b(1,j), 1)
            enddo
            call ga_put(g_a, i, i, 1, nvec, dbl_mb(k_c), 1)
         enddo
c
         if (.not. ma_pop_stack(l_a)) call errquit('ga_mix: pop?',0,
     &       MA_ERR)
         if (.not. ma_pop_stack(l_c)) call errquit('ga_mix: pop?',0,
     &       MA_ERR)
       endif
c
       call ga_sync()
c
       end

