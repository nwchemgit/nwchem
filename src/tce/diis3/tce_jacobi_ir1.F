      subroutine tce_jacobi_ir1(d_r1,d_t1,k_t1_offset,omega,
     1                          shift,damping)
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
#include "tce_diis.fh"
      integer d_r1
      integer d_t1
      integer p1b
      integer h2b
      integer p1
      integer h2
      integer k_t1_offset
      integer size
      integer l_r1,k_r1
      integer i
      integer nprocs
      integer count
      integer next
      integer nxtask
      external nxtask
      logical noloadbalance
      logical nodezero         ! True if node 0
      double precision shift,omega,damping
      double precision denom
c
c     ================
c     Loop over blocks
c     ================
c
      nodezero = (ga_nodeid().eq.0)
      noloadbalance = ((ioalg.eq.4).or.
     1                ((ioalg.eq.6).and.(.not.fileisga(d_r1))))
      nprocs = ga_nnodes()
      count = 0
      next = nxtask(nprocs,1)
      do p1b = noab+1,noab+nvab
        do h2b = 1,noab
          if (noloadbalance.or.(next.eq.count)) then
            if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h2b-1)) then 
            if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
     1        +int_mb(k_spin+h2b-1).ne.4)) then
            if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h2b-1)).eq.
     &          irrep_x) then
              size = int_mb(k_range+p1b-1) * int_mb(k_range+h2b-1)
              if (.not.ma_push_get(mt_dbl,size,'rr1',l_r1,k_r1))
     1          call errquit('tce_jacobi_ir1: MA problem',0,MA_ERR)
              call get_hash_block(d_r1,dbl_mb(k_r1),size,
     1          int_mb(k_t1_offset),((p1b-noab-1)*noab+h2b-1))
              i = 0
              do p1 = 1,int_mb(k_range+p1b-1)
               do h2 = 1,int_mb(k_range+h2b-1)
                i = i + 1
                denom = ( 
     1             -dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1)
     2             +dbl_mb(k_evl_sorted+int_mb(k_offset+h2b-1)+h2-1) 
     3                   )
                dbl_mb(k_r1+i-1) = damping*dbl_mb(k_r1+i-1)
     1                           / ( shift + denom - omega )
               enddo
              enddo
              call add_hash_block(d_t1,dbl_mb(k_r1),size,
     1          int_mb(k_t1_offset),((p1b-noab-1)*noab+h2b-1))
              if (nodezero.and.util_print('tr1',print_debug)) then
                call get_hash_block(d_t1,dbl_mb(k_r1),size,
     1          int_mb(k_t1_offset),((p1b-noab-1)*noab+h2b-1))
                call ma_print_compact(dbl_mb(k_r1),size,1,'tr1')
              endif
              if (.not.ma_pop_stack(l_r1))
     1          call errquit('tce_jacobi_ir1: MA problem',1,MA_ERR)
            endif
            endif
            endif
            next = nxtask(nprocs,1)
          endif
          count = count + 1
        enddo
      enddo
      next = nxtask(-nprocs,1)
      call ga_sync()
      return
      end
