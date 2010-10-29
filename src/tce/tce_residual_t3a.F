      subroutine tce_residual_t3a(d_r3,k_r3_offset,residual)
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
      integer d_r3
      integer p1b
      integer p2b
      integer p3b
      integer h4b
      integer h5b
      integer h6b
      integer k_r3_offset
      integer size
      integer l_r3,k_r3
      integer g_residual
      double precision factor
      double precision residual
      integer nprocs
      integer count
      integer next
      integer nxtask
      external nxtask
      logical nodezero
      logical acolo
c
c     new
c
      if (fileisga(d_r3)) then
        residual = ga_ddot(d_r3,d_r3)
        return
      endif
c
c     =====================
c     Zero scratch residual
c     =====================
c
      nodezero = (ga_nodeid().eq.0)
      if (.not.ga_create(mt_dbl,1,1,'residual',1,1,g_residual))
     1  call errquit('tce_residual_t3: GA problem',0,GA_ERR)
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
      next = nxtask(nprocs,1)
      do p1b = noab+1,noab+nvab
        do p2b = p1b,noab+nvab
          do p3b = p2b,noab+nvab
            do h4b = 1,noab
              do h5b = h4b,noab
                do h6b = h5b,noab
                  if (next.eq.count) then
                    if(acolo(p1b,p2b,p3b,h4b,h5b,h6b)) then
                    if (int_mb(k_spin+p1b-1)
     1                 +int_mb(k_spin+p2b-1)
     2                 +int_mb(k_spin+p3b-1)
     3              .eq.int_mb(k_spin+h4b-1)
     4                 +int_mb(k_spin+h5b-1)
     5                 +int_mb(k_spin+h6b-1)) then
                    if ((.not.restricted).or.
     1                 (int_mb(k_spin+p1b-1)
     2                 +int_mb(k_spin+p2b-1)
     3                 +int_mb(k_spin+p3b-1)
     4                 +int_mb(k_spin+h4b-1)
     5                 +int_mb(k_spin+h5b-1)
     6                 +int_mb(k_spin+h6b-1).ne.12)) then
                    if (ieor(int_mb(k_sym+p1b-1),
     1                  ieor(int_mb(k_sym+p2b-1),
     2                  ieor(int_mb(k_sym+p3b-1),
     3                  ieor(int_mb(k_sym+h4b-1),
     4                  ieor(int_mb(k_sym+h5b-1),
     5                       int_mb(k_sym+h6b-1)))))) .eq. 0) then
                      factor = 1.0d0
                      if ((p3b .gt. p2b).and.(p2b .gt. p1b)) then
                        factor = factor * 6.0d0
                      else if (p3b .gt. p2b) then
                        factor = factor * 2.0d0
                      else if (p2b .gt. p1b) then
                        factor = factor * 2.0d0
                      endif
                      if ((h6b .gt. h5b).and.(h5b .gt. h4b)) then
                        factor = factor * 6.0d0
                      else if (h6b .gt. h5b) then
                        factor = factor * 2.0d0
                      else if (h5b .gt. h4b) then
                        factor = factor * 2.0d0
                      endif
                      size = int_mb(k_range+p1b-1)
     1                     * int_mb(k_range+p2b-1)
     2                     * int_mb(k_range+p3b-1)
     3                     * int_mb(k_range+h4b-1)
     4                     * int_mb(k_range+h5b-1)
     5                     * int_mb(k_range+h6b-1)
                      if (.not.ma_push_get(mt_dbl,size,'r3',
     1                  l_r3,k_r3))
     2                  call errquit('tce_residual_t3: MA problem',0,
     3                  MA_ERR)
                      call get_hash_block(d_r3,dbl_mb(k_r3),size,
     1                  int_mb(k_r3_offset),((((((p1b-noab-1)*nvab
     2                  +p2b-noab-1)*nvab+p3b-noab-1)*noab
     3                  +h4b-1)*noab+h5b-1)*noab+h6b-1))
                      if (nodezero.and.util_print('residual',
     1                  print_debug)) call ma_print_compact
     2                  (dbl_mb(k_r3),size,1,'t3 residual')
                      residual = ddot(size,dbl_mb(k_r3),1,
     1                                     dbl_mb(k_r3),1)
                      call ga_acc(g_residual,1,1,1,1,
     1                              residual,1,factor)
                      if (.not.ma_pop_stack(l_r3))
     1                  call errquit('tce_residual_t3: MA problem',1,
     2                  MA_ERR)
                    endif
                    endif
                    endif
                    endif         !active
                    next = nxtask(nprocs,1)
                  endif
                  count = count + 1
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      next = nxtask(-nprocs,1)
      call ga_sync()
      call ga_get(g_residual,1,1,1,1,residual,1)
      residual = 0.5d0*dsqrt(residual)
      if (.not.ga_destroy(g_residual))
     1  call errquit('tce_residual_t3: GA problem',1,GA_ERR)
      return
      end
