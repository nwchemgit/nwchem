      subroutine ga_local_mdot(n, nva, nvb, c, ldc, g_a, g_b)
*
* $Id: ga_local_mdot.F,v 1.4 2000-10-26 19:21:44 edo Exp $
*
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer n, nva, nvb
      integer ldc
      double precision c(ldc,*)
      integer g_a, g_b
c
c     c(i,j) = sum(k=1,n) a(k,i) b(k,j) or C = At . B
c
c     A and B identically distributed down the columns.  All 
c     of the second index local for both
c
      integer ailo, aihi, ajlo, ajhi, aind, ald, aoff
      integer bilo, bihi, bjlo, bjhi, bind, bld, boff
      integer i, j, k, ij, kdim
      double precision sum, work(2500)
c
      call ga_distribution(g_a, ga_nodeid(), ailo, aihi, ajlo, ajhi)
      call ga_distribution(g_b, ga_nodeid(), bilo, bihi, bjlo, bjhi)
c
      if (ailo.ne.bilo .or. aihi.ne.bihi) call errquit
     $     ('ga_l_mdot: distribution not the same?',0)
c
      call ga_sync
c
      if (aihi.gt.0 .and. aihi.ge.ailo) then
         if (ajlo.ne.1 .or. ajhi.lt.nva) call errquit('ga_l_mdot:a2?',0)
         if (bjlo.ne.1 .or. bjhi.lt.nvb) call errquit('ga_l_mdot:b2?',0)
c
         call ga_access(g_a, ailo, aihi, ajlo, ajhi, aind, ald)
         call ga_access(g_b, bilo, bihi, bjlo, bjhi, bind, bld)
         call ga_lmdotb(work,dbl_mb(aind),dbl_mb(bind),
     ,        aihi,ailo,ald,nva,nvb,bld)
         call ga_release(g_a, ailo, aihi, ajlo, ajhi)
         call ga_release(g_b, bilo, bihi, bjlo, bjhi)
      else
         call dfill(nva*nvb, 0.0d0, work, 1)
      endif
c
      call ga_sync
c
      call ga_dgop(1, work, nva*nvb, '+')
c
      do j = 1, nvb
         do i = 1, nva
            ij   = (j-1)*nva + i
            c(i,j) = work(ij)
         enddo
      enddo
c
c
      end
      subroutine ga_lmdotb(work,a,b,aihi,ailo,ald,nva,nvb,bld)
      implicit none
#include "mafdecls.fh"
      double precision a(*),b(*),work(*)
      integer ailo,aihi,ald,nva,nvb,bld
c     
      integer i,j,k,aoff,boff,ij,kdim
      double precision sum
c
      kdim = aihi - ailo + 1
      do j = 1, nvb
         do i = 1, nva
            aoff = (i-1)*ald 
            boff = (j-1)*bld
            ij   = (j-1)*nva + i
            if (ij.gt.2500) call errquit('Robert is lazy',0)
            sum = 0.0d0
            do k = 1, kdim
               sum = sum + a(aoff+k)*b(boff+k)
            enddo
            work(ij) = sum
         enddo
      enddo
      return
      end
