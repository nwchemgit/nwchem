      subroutine tce_jacobi_l2(d_r2,d_l2,k_l2_offset,shift)
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
      integer d_l2
      integer p1b
      integer p2b
      integer h3b
      integer h4b
      integer p1
      integer p2
      integer h3
      integer h4
      integer k_l2_offset
      integer size
      integer l_r2,k_r2
      integer i
      integer nprocs
      integer count
      integer next
      integer nxtask
      external nxtask
      logical nodezero
      logical noloadbalance
      double precision shift
c
      nodezero = (ga_nodeid().eq.0)
      noloadbalance = ((ioalg.eq.4).or.
     1                ((ioalg.eq.6).and.(.not.fileisga(d_l2)))) 
      nprocs = ga_nnodes()
      count = 0
      next = nxtask(nprocs,1)
      do h3b = 1,noab
        do h4b = h3b,noab
          do p1b = noab+1,noab+nvab
            do p2b = p1b,noab+nvab
              if (noloadbalance.or.(next.eq.count)) then
                if (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) 
     1            .eq. int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) then
                if ((.not.restricted).or.
     1            (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+
     2            int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) then
                if (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),
     1            ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1))))
     2            .eq. irrep_y) then
                  size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
     1                 * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
                  if (.not.ma_push_get(mt_dbl,size,'r2',l_r2,k_r2))
     1              call errquit('tce_jacobi_l2: MA problem',0,MA_ERR)
                  call get_hash_block(d_r2,dbl_mb(k_r2),size,
     1              int_mb(k_l2_offset),((((h3b-1)*noab+h4b-1)
     2              *nvab+p1b-noab-1)*nvab+p2b-noab-1))
                  i = 0
                  do h3 = 1,int_mb(k_range+h3b-1)
                    do h4 = 1,int_mb(k_range+h4b-1)
                      do p1 = 1,int_mb(k_range+p1b-1)
                        do p2 = 1,int_mb(k_range+p2b-1)
                          i = i + 1
                          dbl_mb(k_r2+i-1) = dbl_mb(k_r2+i-1)
     1           / (-dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1)
     2              -dbl_mb(k_evl_sorted+int_mb(k_offset+p2b-1)+p2-1)
     3              +dbl_mb(k_evl_sorted+int_mb(k_offset+h3b-1)+h3-1)
     4              +dbl_mb(k_evl_sorted+int_mb(k_offset+h4b-1)+h4-1)
     5             +shift)
                        enddo
                      enddo
                    enddo
                  enddo
                  call add_hash_block(d_l2,dbl_mb(k_r2),size,
     1              int_mb(k_l2_offset),((((h3b-1)*noab+h4b-1)
     2              *nvab+p1b-noab-1)*nvab+p2b-noab-1))
                  if (.not.ma_pop_stack(l_r2))
     1              call errquit('tce_jacobi_l2: MA problem',1,MA_ERR)
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
