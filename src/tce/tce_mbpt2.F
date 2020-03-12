      subroutine tce_mbpt2(d_mo2e,k_2e_offset,
     1                     d_t2,k_t2_offset,mbpt2)
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
      integer d_mo2e
      integer d_t2
      integer p1b
      integer p2b
      integer h3b
      integer h4b
      integer k_2e_offset
      integer k_t2_offset
      integer size
      integer l_t2,k_t2
      integer l_v2,k_v2
      integer g_mbpt2
      double precision factor
      double precision mbpt2
      integer nprocs
      integer count
      integer next
cc      integer nxtask
cc      external nxtask
      INTEGER NXTASK
      EXTERNAL NXTASK
      double precision cpu
      double precision wall
      logical nodezero
c
      nodezero=(ga_nodeid().eq.0)
      cpu=-util_cpusec()
      wall=-util_wallsec()
c
c     =====================
c     Zero scratch residual
c     =====================
c
      if (.not.ga_create(mt_dbl,1,1,'mbpt2',1,1,g_mbpt2))
     1  call errquit('tce_mbpt2: GA problem',0,GA_ERR)
      mbpt2 = 0.0d0
      call ga_put(g_mbpt2,1,1,1,1,mbpt2,1)
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
                if ((.not.restricted).or.
     1             (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)
     2             +int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8))
     3             then
                if (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) 
     1            .eq. int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) then
                if (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),
     1            ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1))))
     2            .eq. 0) then
                  factor = 0.25d0
                  if (p2b .gt. p1b) factor = factor * 2.0d0
                  if (h4b .gt. h3b) factor = factor * 2.0d0
                if ((restricted).and.
     1             (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)
     2             +int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).eq.4))
     3             factor = factor * 2.0d0
                  size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
     1                 * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
                  if (.not.ma_push_get(mt_dbl,size,'v2',l_v2,k_v2))
     1              call errquit('tce_mbpt2: MA problem',0,MA_ERR)
                  call get_hash_block(d_mo2e,dbl_mb(k_v2),size,
     1              int_mb(k_2e_offset),((((p1b-1)*(noab+nvab)+p2b-1)
     2              *(noab+nvab)+h3b-1)*(noab+nvab)+h4b-1))
                  if (.not.ma_push_get(mt_dbl,size,'t2',l_t2,k_t2))
     1              call errquit('tce_mbpt2: MA problem',1,MA_ERR)
                  call get_hash_block(d_t2,dbl_mb(k_t2),size,
     1              int_mb(k_t2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
                  mbpt2 = ddot(size,dbl_mb(k_t2),1,dbl_mb(k_v2),1)
                  if (util_print('mbpt2',print_debug))
     1              write(LuOut,9200) ga_nodeid(),p1b,p2b,h3b,h4b,
     2              mbpt2 * factor
                  call ga_acc(g_mbpt2,1,1,1,1,mbpt2,1,factor)
                  if (.not.ma_pop_stack(l_t2))
     1              call errquit('tce_mbpt2: MA problem',2,MA_ERR)
                  if (.not.ma_pop_stack(l_v2))
     1              call errquit('tce_mbpt2: MA problem',3,MA_ERR)
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
      call ga_get(g_mbpt2,1,1,1,1,mbpt2,1)
      if (.not.ga_destroy(g_mbpt2))
     1  call errquit('tce_mbpt2: GA problem',1,GA_ERR)
      cpu=cpu+util_cpusec()
      wall=wall+util_wallsec()
      if (nodezero) write(LuOut,9020) cpu, wall
 9020 format(1x,'Cpu & wall time / sec',2f15.1)
 9200 format(/,1x,'node',i3,' blocks:',4i3,' MBPT(2) contribution',
     1       f20.10)
      return
      end
