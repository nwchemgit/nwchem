      subroutine tce_residual_t2(d_r2,k_r2_offset,residual)
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
      integer d_r2
      integer p1b
      integer p2b
      integer h3b
      integer h4b
      integer k_r2_offset
      integer size
      integer l_r2,k_r2
      integer g_residual
      double precision factor
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
      if (fileisga(d_r2)) then
        residual = ga_ddot(d_r2,d_r2)
        return
      endif
c
c     =====================
c     Zero scratch residual
c     =====================
c
      nodezero = (ga_nodeid().eq.0)
      if (.not.ga_create(mt_dbl,1,1,'residual',1,1,g_residual))
     1  call errquit('tce_residual_t2: GA problem',0,GA_ERR)
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
        do p2b = p1b,noab+nvab
          do h3b = 1,noab
            do h4b = h3b,noab
              if (next.eq.count) then
                if (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) 
     1            .eq. int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) then
                if ((.not.restricted).or.
     1            (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+
     2            int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) then
                if (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),
     1            ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1))))
     2            .eq. 0) then
                  factor = 1.0d0
                  if (p2b .gt. p1b) factor = factor * 2.0d0
                  if (h4b .gt. h3b) factor = factor * 2.0d0
                  size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
     1                 * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
                  if (.not.ma_push_get(mt_dbl,size,'r2',l_r2,k_r2))
     1              call errquit('tce_residual_t2: MA problem',0,MA_ERR)
                  call get_hash_block(d_r2,dbl_mb(k_r2),size,
     1              int_mb(k_r2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
                  if (nodezero.and.util_print('residual',
     1              print_debug)) call ma_print_compact
     2              (dbl_mb(k_r2),size,1,'t2 residual')
                  residual = ddot(size,dbl_mb(k_r2),1,dbl_mb(k_r2),1)
                  call ga_acc(g_residual,1,1,1,1,residual,1,factor)
                  if (.not.ma_pop_stack(l_r2))
     1              call errquit('tce_residual_t2: MA problem',1,MA_ERR)
                endif
                endif
                endif
cc                next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
              endif
              count = count + 1
            enddo
          enddo
        enddo
      enddo
cc      next = nxtask(-nprocs,1)
      next = NXTASK(-nprocs, 1)
      call ga_sync()
      call ga_get(g_residual,1,1,1,1,residual,1)
      residual = 0.5d0*dsqrt(residual)
      if (.not.ga_destroy(g_residual))
     1  call errquit('tce_residual_t2: GA problem',1,GA_ERR)
      return
      end
      subroutine tce_residual_tr2(d_r2,k_r2_offset,residual)
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
      integer d_r2
      integer p1b
      integer p2b
      integer h3b
      integer h4b
      integer k_r2_offset
      integer size
      integer l_r2,k_r2
      integer g_residual
      double precision factor
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
      if (fileisga(d_r2)) then
        residual = ga_ddot(d_r2,d_r2)
        return
      endif
c
c     =====================
c     Zero scratch residual
c     =====================
c
      nodezero = (ga_nodeid().eq.0)
      if (.not.ga_create(mt_dbl,1,1,'residual',1,1,g_residual))
     1  call errquit('tce_residual_tr2: GA problem',0,GA_ERR)
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
        do p2b = p1b,noab+nvab
          do h3b = 1,noab
            do h4b = h3b,noab
              if (next.eq.count) then
                if (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) 
     1            .eq. int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) then
                if ((.not.restricted).or.
     1            (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+
     2            int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) then
                if (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),
     1            ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1))))
     2            .eq. irrep_x) then
                  factor = 1.0d0
                  if (p2b .gt. p1b) factor = factor * 2.0d0
                  if (h4b .gt. h3b) factor = factor * 2.0d0
                  size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
     1                 * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
                  if (.not.ma_push_get(mt_dbl,size,'rr2',l_r2,k_r2))
     1             call errquit('tce_residual_tr2: MA problem',0,MA_ERR)
                  call get_hash_block(d_r2,dbl_mb(k_r2),size,
     1              int_mb(k_r2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
                  if (nodezero.and.util_print('residual',
     1              print_debug)) call ma_print_compact
     2              (dbl_mb(k_r2),size,1,'tr2 residual')
                  residual = ddot(size,dbl_mb(k_r2),1,dbl_mb(k_r2),1)
                  call ga_acc(g_residual,1,1,1,1,residual,1,factor)
                  if (.not.ma_pop_stack(l_r2))
     1             call errquit('tce_residual_tr2: MA problem',1,MA_ERR)
                endif
                endif
                endif
cc                next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
              endif
              count = count + 1
            enddo
          enddo
        enddo
      enddo
cc      next = nxtask(-nprocs,1)
      next = NXTASK(-nprocs, 1)
      call ga_sync()
      call ga_get(g_residual,1,1,1,1,residual,1)
      residual = 0.5d0*dsqrt(residual)
      if (.not.ga_destroy(g_residual))
     1  call errquit('tce_residual_tr2: GA problem',1,GA_ERR)
      return
      end
c
c
c
c
      subroutine tce_residual_t2a(d_r2,k_r2_offset,residual)
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
      integer d_r2
      integer p1b
      integer p2b
      integer h3b
      integer h4b
      integer k_r2_offset
      integer size
      integer l_r2,k_r2
      integer g_residual
      double precision factor
      double precision residual
      LOGICAL is_active_1,is_active_2,is_active_3,is_active_4
      integer nprocs
      integer count
      integer next
cc      integer nxtask
cc      external nxtask
      INTEGER NXTASK
      EXTERNAL NXTASK
      logical nodezero
c
c     =====================
c     Zero scratch residual
c     =====================
c
      nodezero = (ga_nodeid().eq.0)
      if (.not.ga_create(mt_dbl,1,1,'residual',1,1,g_residual))
     1  call errquit('tce_residual_t2: GA problem',0,GA_ERR)
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
        do p2b = p1b,noab+nvab
          do h3b = 1,noab
            do h4b = h3b,noab
              if (next.eq.count) then
              IF(is_active_4(p1b,p2b,h3b,h4b)) THEN
                if (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) 
     1            .eq. int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) then
                if ((.not.restricted).or.
     1            (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+
     2            int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) then
                if (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),
     1            ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1))))
     2            .eq. 0) then
                  factor = 1.0d0
                  if (p2b .gt. p1b) factor = factor * 2.0d0
                  if (h4b .gt. h3b) factor = factor * 2.0d0
                  size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
     1                 * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
                  if (.not.ma_push_get(mt_dbl,size,'r2',l_r2,k_r2))
     1              call errquit('tce_residual_t2: MA problem',0,MA_ERR)
                  call get_hash_block(d_r2,dbl_mb(k_r2),size,
     1              int_mb(k_r2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
                  if (nodezero.and.util_print('residual',
     1              print_debug)) call ma_print_compact
     2              (dbl_mb(k_r2),size,1,'t2 residual')
                  residual = ddot(size,dbl_mb(k_r2),1,dbl_mb(k_r2),1)
                  call ga_acc(g_residual,1,1,1,1,residual,1,factor)
                  if (.not.ma_pop_stack(l_r2))
     1              call errquit('tce_residual_t2: MA problem',1,MA_ERR)
                endif
                endif
                endif
                endif
cc                next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
              endif
              count = count + 1
            enddo
          enddo
        enddo
      enddo
cc      next = nxtask(-nprocs,1)
      next = NXTASK(-nprocs, 1)
      call ga_sync()
      call ga_get(g_residual,1,1,1,1,residual,1)
      residual = 0.5d0*dsqrt(residual)
      if (.not.ga_destroy(g_residual))
     1  call errquit('tce_residual_t2: GA problem',1,GA_ERR)
      return
      end

