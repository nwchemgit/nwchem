      subroutine tce_jacobi_l3(d_r3,d_l3,k_l3_offset,shift)
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
      integer d_l3
      integer p1b
      integer p2b
      integer p3b
      integer h4b
      integer h5b
      integer h6b
      integer p1
      integer p2
      integer p3
      integer h4
      integer h5
      integer h6
      integer k_l3_offset
      integer size
      integer l_r3,k_r3
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
     1                ((ioalg.eq.6).and.(.not.fileisga(d_l3))))
      nprocs = ga_nnodes()
      count = 0
      next = nxtask(nprocs,1)
      do h4b = 1,noab
        do h5b = h4b,noab
          do h6b = h5b,noab
            do p1b = noab+1,noab+nvab
              do p2b = p1b,noab+nvab
                do p3b = p2b,noab+nvab
                  if (noloadbalance.or.(next.eq.count)) then
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
     5                  int_mb(k_sym+h6b-1)))))) .eq. irrep_y) then
                      size = int_mb(k_range+p1b-1)
     1                     * int_mb(k_range+p2b-1)
     2                     * int_mb(k_range+p3b-1)
     3                     * int_mb(k_range+h4b-1)
     4                     * int_mb(k_range+h5b-1)
     5                     * int_mb(k_range+h6b-1)
                      if (.not.ma_push_get(mt_dbl,size,'r3',
     1                  l_r3,k_r3))
     2                  call errquit('tce_jacobi_l3: MA problem',0,
     3                  MA_ERR)
                      call get_hash_block(d_r3,dbl_mb(k_r3),size,
     1                  int_mb(k_l3_offset),((((((h4b-1)*noab
     2                  +h5b-1)*noab+h6b-1)*nvab+p1b-noab-1)*nvab
     3                  +p2b-noab-1)*nvab+p3b-noab-1))
                      i = 0
                      do h4 = 1,int_mb(k_range+h4b-1)
                        do h5 = 1,int_mb(k_range+h5b-1)
                          do h6 = 1,int_mb(k_range+h6b-1)
                            do p1 = 1,int_mb(k_range+p1b-1)
                              do p2 = 1,int_mb(k_range+p2b-1)
                                do p3 = 1,int_mb(k_range+p3b-1)
                                  i = i + 1
                     dbl_mb(k_r3+i-1) = dbl_mb(k_r3+i-1)
     1           / (-dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1)
     2              -dbl_mb(k_evl_sorted+int_mb(k_offset+p2b-1)+p2-1)
     3              -dbl_mb(k_evl_sorted+int_mb(k_offset+p3b-1)+p3-1)
     4              +dbl_mb(k_evl_sorted+int_mb(k_offset+h4b-1)+h4-1)
     5              +dbl_mb(k_evl_sorted+int_mb(k_offset+h5b-1)+h5-1)
     6              +dbl_mb(k_evl_sorted+int_mb(k_offset+h6b-1)+h6-1)
     7             +shift)
                                enddo
                              enddo
                            enddo
                          enddo
                        enddo
                      enddo
                      call add_hash_block(d_l3,dbl_mb(k_r3),size,
     1                  int_mb(k_l3_offset),((((((h4b-1)*noab
     2                  +h5b-1)*noab+h6b-1)*nvab+p1b-noab-1)*nvab
     3                  +p2b-noab-1)*nvab+p3b-noab-1))
                      if (.not.ma_pop_stack(l_r3))
     1                  call errquit('tce_jacobi_l3: MA problem',1,
     2                  MA_ERR)
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
        enddo
      enddo
      next = nxtask(-nprocs,1)
      call ga_sync()
      return
      end
