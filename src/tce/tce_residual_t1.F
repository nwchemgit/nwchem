      subroutine tce_residual_t1(d_r1,k_r1_offset,residual)
c
c $Id$
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer d_r1
      integer p1b
      integer h2b
      integer k_r1_offset
      integer size
      integer l_r1,k_r1
      integer g_residual
      double precision residual
      integer nprocs
      integer count
      integer next
cc      integer nxtask
cc      external nxtask
      INTEGER NXTASK
      EXTERNAL NXTASK
      logical nodezero
c
c     new
c
      if (fileisga(d_r1)) then
        residual = ga_ddot(d_r1,d_r1)
        return
      endif
c
c     =====================
c     Zero scratch residual
c     =====================
c
      nodezero = (ga_nodeid().eq.0)
      if (.not.ga_create(mt_dbl,1,1,'residual',1,1,g_residual))
     1  call errquit('tce_residual_t1: GA problem',0,GA_ERR)
      residual = 0.0d0
      call ga_put(g_residual,1,1,1,1,residual,1)
      call ga_sync()
c
c     ================
c     Loop over blocks
c     ================
c
      nprocs = ga_nnodes()
      count = 0
cc      next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
      do p1b = noab+1,noab+nvab
        do h2b = 1,noab
          if (next.eq.count) then
            if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h2b-1)) then
            if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
     1        +int_mb(k_spin+h2b-1).ne.4)) then
            if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h2b-1)).eq.0) then
              size = int_mb(k_range+p1b-1) * int_mb(k_range+h2b-1)
              if (.not.ma_push_get(mt_dbl,size,'r1',l_r1,k_r1))
     1          call errquit('tce_residual_t1: MA problem',0,MA_ERR)
              call get_hash_block(d_r1,dbl_mb(k_r1),size,
     1          int_mb(k_r1_offset),((p1b-noab-1)*noab+h2b-1))
              if (nodezero.and.util_print('residual',
     1          print_debug)) call ma_print_compact
     2          (dbl_mb(k_r1),size,1,'t1 residual')
              residual = ddot(size,dbl_mb(k_r1),1,dbl_mb(k_r1),1)
              call ga_acc(g_residual,1,1,1,1,residual,1,1.0d0)
              if (.not.ma_pop_stack(l_r1))
     1          call errquit('tce_residual_t1: MA problem',1,MA_ERR)
            endif
            endif
            endif
cc            next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
          endif
          count = count + 1
        enddo
      enddo
cc      next = nxtask(-nprocs,1)
      next = NXTASK(-nprocs, 1)
      call ga_sync()
      call ga_get(g_residual,1,1,1,1,residual,1)
      residual = 0.5d0*dsqrt(residual)
      if (.not.ga_destroy(g_residual))
     1  call errquit('tce_residual_t1: GA problem',1,GA_ERR)
      return
      end










      subroutine tce_residual_tr1(d_r1,k_r1_offset,residual)
c
c $Id$
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer d_r1
      integer p1b
      integer h2b
      integer k_r1_offset
      integer size
      integer l_r1,k_r1
      integer g_residual
      double precision residual
      integer nprocs
      integer count
      integer next
cc      integer nxtask
cc      external nxtask
      INTEGER NXTASK
      EXTERNAL NXTASK
      logical nodezero
c
c     new
c
      if (fileisga(d_r1)) then
        residual = ga_ddot(d_r1,d_r1)
        return
      endif
c
c     =====================
c     Zero scratch residual
c     =====================
c
      nodezero = (ga_nodeid().eq.0)
      if (.not.ga_create(mt_dbl,1,1,'residual',1,1,g_residual))
     1  call errquit('tce_residual_tr1: GA problem',0,GA_ERR)
      residual = 0.0d0
      call ga_put(g_residual,1,1,1,1,residual,1)
      call ga_sync()
c
c     ================
c     Loop over blocks
c     ================
c
      nprocs = ga_nnodes()
      count = 0
cc      next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
      do p1b = noab+1,noab+nvab
        do h2b = 1,noab
          if (next.eq.count) then
            if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h2b-1)) then
            if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
     1        +int_mb(k_spin+h2b-1).ne.4)) then
            if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h2b-1)).eq.
     1        irrep_x) then
              size = int_mb(k_range+p1b-1) * int_mb(k_range+h2b-1)
              if (.not.ma_push_get(mt_dbl,size,'rr1',l_r1,k_r1))
     1          call errquit('tce_residual_t1: MA problem',0,MA_ERR)
              call get_hash_block(d_r1,dbl_mb(k_r1),size,
     1          int_mb(k_r1_offset),((p1b-noab-1)*noab+h2b-1))
              if (nodezero.and.util_print('residual',
     1          print_debug)) call ma_print_compact
     2          (dbl_mb(k_r1),size,1,'tr1 residual')
              residual = ddot(size,dbl_mb(k_r1),1,dbl_mb(k_r1),1)
              call ga_acc(g_residual,1,1,1,1,residual,1,1.0d0)
              if (.not.ma_pop_stack(l_r1))
     1          call errquit('tce_residual_tr1: MA problem',1,MA_ERR)
            endif
            endif
            endif
cc            next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
          endif
          count = count + 1
        enddo
      enddo
cc      next = nxtask(-nprocs,1)
      next = NXTASK(-nprocs, 1)
      call ga_sync()
      call ga_get(g_residual,1,1,1,1,residual,1)
      residual = 0.5d0*dsqrt(residual)
      if (.not.ga_destroy(g_residual))
     1  call errquit('tce_residual_tr1: GA problem',1,GA_ERR)
      return
      end
