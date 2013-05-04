      subroutine tce_jacobi2_ir2(d_r2,d_t2,k_t2_offset,omega,
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
      integer d_r2
      integer d_t2
      integer p1b
      integer p2b
      integer h3b
      integer h4b
      integer p1
      integer p2
      integer h3
      integer h4
      integer k_t2_offset
      integer size
      integer l_r2,k_r2
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
     1                ((ioalg.eq.6).and.(.not.fileisga(d_r2))))
      nprocs = ga_nnodes()
      count = 0
      next = nxtask(nprocs,1)
      do p1b = noab+1,noab+nvab
        do p2b = p1b,noab+nvab
          do h3b = 1,noab
            do h4b = h3b,noab
              if (noloadbalance.or.(next.eq.count)) then
                if (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) 
     1            .eq. int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) then
                if ((.not.restricted).or.
     1            (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+
     2            int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) then
                if (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),
     1            ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1))))
     2            .eq. irrep_x) then
                  size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
     1                 * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
                  if (.not.ma_push_get(mt_dbl,size,'rr2',l_r2,k_r2))
     1              call errquit('tce_jacobi2_ir2: MA problem',0,MA_ERR)
                  call get_hash_block(d_r2,dbl_mb(k_r2),size,
     1              int_mb(k_t2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
                  i = 0
                  do p1 = 1,int_mb(k_range+p1b-1)
                   do p2 = 1,int_mb(k_range+p2b-1)
                    do h3 = 1,int_mb(k_range+h3b-1)
                     do h4 = 1,int_mb(k_range+h4b-1)
                       i = i + 1
                denom = ( 
     1             -dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1)
     2             -dbl_mb(k_evl_sorted+int_mb(k_offset+p2b-1)+p2-1)
     3             +dbl_mb(k_evl_sorted+int_mb(k_offset+h3b-1)+h3-1) 
     4             +dbl_mb(k_evl_sorted+int_mb(k_offset+h4b-1)+h4-1) 
     5                   )
                dbl_mb(k_r2+i-1) = damping*dbl_mb(k_r2+i-1)
     1                           / ( shift - denom*denom - omega*omega )
                     enddo
                    enddo
                   enddo
                  enddo
                  call add_hash_block(d_t2,dbl_mb(k_r2),size,
     1              int_mb(k_t2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
                  if (nodezero.and.util_print('tr2',print_debug)) 
     1              then
                    call get_hash_block(d_t2,dbl_mb(k_r2),size,
     1              int_mb(k_t2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
                    call ma_print_compact
     2              (dbl_mb(k_r2),size,1,'tr2')
                  endif
                  if (.not.ma_pop_stack(l_r2))
     1              call errquit('tce_jacobi2_ir2: MA problem',1,MA_ERR)
                endif
                endif
                endif
                next = nxtask(nprocs,1)
              endif
              count = count + 1
            enddo
          enddo
        enddo
      enddo
      next = nxtask(-nprocs,1)
      call ga_sync()
      return
      end
