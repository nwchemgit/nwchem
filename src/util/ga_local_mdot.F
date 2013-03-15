      subroutine ga_local_mdot(n, nva, nvb, c, ldc, g_a, g_b)
*
* $Id$
*
      implicit none
#include "errquit.fh"
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
      integer ailo, aihi, ajlo, ajhi, aind, ald
      integer bilo, bihi, bjlo, bjhi, bind, bld
      integer i, j, ij
      integer k_work,l_work
      integer dim1,dim2,type
#ifdef  BAD_GACCESS 
      integer g_c
      if(.not.ga_duplicate(g_b, g_c, ' g_b clone')) call
     $     errquit('galocalmdot: gaduplicate  failed ',0)
      call ga_inquire(g_a, type, dim1, dim2)
      call ga_inquire(g_b, type, dim1, dim2)
      call ga_matmul_patch('t','n',1d0,0d0,
     1     g_a,1,nva,1,dim1,
     2     g_b,1,dim1,1,nvb,
     3     g_c,1,nva ,1,nvb )
      call ga_get(g_c,1,nva,1,nvb,c,ldc)
      if(.not.ga_destroy(g_c)) call
     $     errquit('galocalmdot: gadestroy  failed ',0)
#else
c
      call ga_distribution(g_a, ga_nodeid(), ailo, aihi, ajlo, ajhi)
      call ga_distribution(g_b, ga_nodeid(), bilo, bihi, bjlo, bjhi)
c
      if (ailo.ne.bilo .or. aihi.ne.bihi) call errquit
     $     ('ga_l_mdot: distribution not the same?',0, GA_ERR)
c
      call ga_sync
c
      if (.not.ma_push_get(mt_dbl,nva*nvb,'work',l_work,k_work))
     1         call errquit('ga_local_mdot: push work',nva*nvb,MA_ERR)
c
      if (aihi.gt.0 .and. aihi.ge.ailo) then
         if (ajlo.ne.1 .or. ajhi.lt.nva) call errquit('ga_l_mdot:a2?',0,
     &       GA_ERR)
         if (bjlo.ne.1 .or. bjhi.lt.nvb) call errquit('ga_l_mdot:b2?',0,
     &       GA_ERR)
c
         call ga_access(g_a, ailo, aihi, ajlo, ajhi, aind, ald)
         call ga_access(g_b, bilo, bihi, bjlo, bjhi, bind, bld)
         call ga_lmdotb(dbl_mb(k_work),dbl_mb(aind),dbl_mb(bind),
     ,        aihi,ailo,ald,nva,nvb,bld)
         call ga_release(g_a, ailo, aihi, ajlo, ajhi)
         call ga_release(g_b, bilo, bihi, bjlo, bjhi)
      else
         call dfill(nva*nvb, 0.0d0, dbl_mb(k_work), 1)
      endif
c
      call ga_sync
c
      call ga_dgop(1, dbl_mb(k_work), nva*nvb, '+')
c
      do j = 1, nvb
         do i = 1, nva
            ij   = (j-1)*nva + i
            c(i,j) = dbl_mb(k_work+ij-1)
         enddo
      enddo
c
      if (.not.ma_pop_stack(l_work))
     1  call errquit('ga_local_mdot: pop work',1,MA_ERR)
c
c
#endif
      end
      subroutine ga_local_zmdot(n, nva, nvb, c, ldc, g_a, g_b)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer n, nva, nvb
      integer ldc
      double complex c(ldc,*)
      double complex one, zero
      integer g_a, g_b
c
c     c(i,j) = sum(k=1,n) a(k,i) b(k,j) or C = At . B
c
c     A and B identically distributed down the columns.  All 
c     of the second index local for both
c
      integer ailo, aihi, ajlo, ajhi, aind, ald
      integer bilo, bihi, bjlo, bjhi, bind, bld
      integer i, j, ij
      integer k_work,l_work
      integer dim1,dim2,type
#ifdef  BAD_GACCESS 
      integer g_c
      one  = complex(1.0d0,0.0d0)
      zero = complex(0.0d0,0.0d0)
      if(.not.ga_duplicate(g_b, g_c, ' g_b clone')) call
     $     errquit('galocalzmdot: gaduplicate  failed ',0)
      call ga_inquire(g_a, type, dim1, dim2)
      call ga_inquire(g_b, type, dim1, dim2)
      call ga_matmul_patch('c','n',one,zero,
     1     g_a,1,nva,1,dim1,
     2     g_b,1,dim1,1,nvb,
     3     g_c,1,nva ,1,nvb )
      call ga_get(g_c,1,nva,1,nvb,c,ldc)
      if(.not.ga_destroy(g_c)) call
     $     errquit('galocalzmdot: gadestroy  failed ',0)
#else
c
      call ga_distribution(g_a, ga_nodeid(), ailo, aihi, ajlo, ajhi)
      call ga_distribution(g_b, ga_nodeid(), bilo, bihi, bjlo, bjhi)
c
      if (ailo.ne.bilo .or. aihi.ne.bihi) call errquit
     $     ('ga_l_mdot: distribution not the same?',0, GA_ERR)
c
      call ga_sync
c
      if (.not.ma_push_get(mt_dcpl,nva*nvb,'work',l_work,k_work))
     1         call errquit('ga_local_mdot: push work',nva*nvb,MA_ERR)
c
      if (aihi.gt.0 .and. aihi.ge.ailo) then
         if (ajlo.ne.1.or.ajhi.lt.nva) call errquit('ga_l_zmdot:a2?',0,
     &       GA_ERR)
         if (bjlo.ne.1.or.bjhi.lt.nvb) call errquit('ga_l_zmdot:b2?',0,
     &       GA_ERR)
c
         call ga_access(g_a, ailo, aihi, ajlo, ajhi, aind, ald)
         call ga_access(g_b, bilo, bihi, bjlo, bjhi, bind, bld)
         call ga_lmzdotb(dcpl_mb(k_work),dcpl_mb(aind),dcpl_mb(bind),
     ,        aihi,ailo,ald,nva,nvb,bld)
         call ga_release(g_a, ailo, aihi, ajlo, ajhi)
         call ga_release(g_b, bilo, bihi, bjlo, bjhi)
      else
         call dfill(nva*nvb, 0.0d0, dcpl_mb(k_work), 1)
      endif
c
      call ga_sync
c
      call ga_zgop(1, dcpl_mb(k_work), nva*nvb, '+')
c
      do j = 1, nvb
         do i = 1, nva
            ij   = (j-1)*nva + i
            c(i,j) = dcpl_mb(k_work+ij-1)
         enddo
      enddo
c
      if (.not.ma_pop_stack(l_work))
     1  call errquit('ga_local_mdot: pop work',1,MA_ERR)
c
c
#endif
      end
      subroutine ga_lmdotb(work,a,b,aihi,ailo,ald,nva,nvb,bld)
      implicit none
#include "errquit.fh"
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
            sum = 0.0d0
            do k = 1, kdim
               sum = sum + a(aoff+k)*b(boff+k)
            enddo
            work(ij) = sum
         enddo
      enddo
      return
      end
      subroutine ga_lmzdotb(work,a,b,aihi,ailo,ald,nva,nvb,bld)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
      double complex a(*),b(*),work(*)
      integer ailo,aihi,ald,nva,nvb,bld
c     
      integer i,j,k,aoff,boff,ij,kdim
      double complex sum
c
      kdim = aihi - ailo + 1
      do j = 1, nvb
         do i = 1, nva
            aoff = (i-1)*ald 
            boff = (j-1)*bld
            ij   = (j-1)*nva + i
            sum = cmplx(0.0d0,0.0d0)
            do k = 1, kdim
               sum = sum + conjg(a(aoff+k))*b(boff+k)
            enddo
            work(ij) = sum
         enddo
      enddo
      return
      end
