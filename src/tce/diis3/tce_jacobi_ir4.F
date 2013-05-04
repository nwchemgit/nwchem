      subroutine tce_jacobi_ir4(d_r4,d_t4,k_t4_offset,omega,
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
      integer d_r4
      integer d_t4
      integer p1b
      integer p2b
      integer p3b
      integer p4b
      integer h5b
      integer h6b
      integer h7b
      integer h8b
      integer p1
      integer p2
      integer p3
      integer p4
      integer h5
      integer h6
      integer h7
      integer h8
      integer k_t4_offset
      integer size
      integer l_r4,k_r4
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
     1                ((ioalg.eq.6).and.(.not.fileisga(d_r4))))
      nprocs = ga_nnodes()
      count = 0
      next = nxtask(nprocs,1)
      do p1b = noab+1,noab+nvab
        do p2b = p1b,noab+nvab
          do p3b = p2b,noab+nvab
            do p4b = p3b,noab+nvab
              do h5b = 1,noab
                do h6b = h5b,noab
                  do h7b = h6b,noab
                    do h8b = h7b,noab
                      if (noloadbalance.or.(next.eq.count)) then
                      if (int_mb(k_spin+p1b-1)
     1                   +int_mb(k_spin+p2b-1)
     2                   +int_mb(k_spin+p3b-1)
     3                   +int_mb(k_spin+p4b-1)
     4                .eq.int_mb(k_spin+h5b-1)
     5                   +int_mb(k_spin+h6b-1)
     6                   +int_mb(k_spin+h7b-1)
     7                   +int_mb(k_spin+h8b-1)) then
                      if ((.not.restricted).or.
     1                   (int_mb(k_spin+p1b-1)
     2                   +int_mb(k_spin+p2b-1)
     3                   +int_mb(k_spin+p3b-1)
     4                   +int_mb(k_spin+p4b-1)
     5                   +int_mb(k_spin+h5b-1)
     6                   +int_mb(k_spin+h6b-1)
     7                   +int_mb(k_spin+h7b-1)
     8                   +int_mb(k_spin+h8b-1).ne.16)) then
                      if (ieor(int_mb(k_sym+p1b-1),
     1                    ieor(int_mb(k_sym+p2b-1),
     2                    ieor(int_mb(k_sym+p3b-1),
     3                    ieor(int_mb(k_sym+p4b-1),
     4                    ieor(int_mb(k_sym+h5b-1),
     5                    ieor(int_mb(k_sym+h6b-1),
     6                    ieor(int_mb(k_sym+h7b-1),
     7                    int_mb(k_sym+h8b-1)))))))) .eq. irrep_x) then
                        size = int_mb(k_range+p1b-1)
     1                       * int_mb(k_range+p2b-1)
     2                       * int_mb(k_range+p3b-1)
     3                       * int_mb(k_range+p4b-1)
     4                       * int_mb(k_range+h5b-1)
     5                       * int_mb(k_range+h6b-1)
     6                       * int_mb(k_range+h7b-1)
     7                       * int_mb(k_range+h8b-1)
                        if (.not.ma_push_get(mt_dbl,size,'rr4',
     1                    l_r4,k_r4))
     2                    call errquit('tce_jacobi_ir4: MA problem',0,
     3                    MA_ERR)
                        call get_hash_block(d_r4,dbl_mb(k_r4),size,
     1                    int_mb(k_t4_offset),((((((((p1b-noab-1)*nvab
     2                    +p2b-noab-1)*nvab+p3b-noab-1)*nvab
     3                    +p4b-noab-1)*noab+h5b-1)*noab
     4                    +h6b-1)*noab+h7b-1)*noab+h8b-1))
                        i = 0
                        do p1 = 1,int_mb(k_range+p1b-1)
                          do p2 = 1,int_mb(k_range+p2b-1)
                            do p3 = 1,int_mb(k_range+p3b-1)
                              do p4 = 1,int_mb(k_range+p4b-1)
                                do h5 = 1,int_mb(k_range+h5b-1)
                                  do h6 = 1,int_mb(k_range+h6b-1)
                                    do h7 = 1,int_mb(k_range+h7b-1)
                                      do h8 = 1,int_mb(k_range+h8b-1)
                                    i = i + 1
                denom = ( 
     1             -dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1)
     2             -dbl_mb(k_evl_sorted+int_mb(k_offset+p2b-1)+p2-1)
     3             -dbl_mb(k_evl_sorted+int_mb(k_offset+p3b-1)+p3-1)
     4             -dbl_mb(k_evl_sorted+int_mb(k_offset+p4b-1)+p4-1)
     5             +dbl_mb(k_evl_sorted+int_mb(k_offset+h5b-1)+h5-1) 
     6             +dbl_mb(k_evl_sorted+int_mb(k_offset+h6b-1)+h6-1) 
     7             +dbl_mb(k_evl_sorted+int_mb(k_offset+h7b-1)+h7-1) 
     8             +dbl_mb(k_evl_sorted+int_mb(k_offset+h8b-1)+h8-1) 
     9                   )
                dbl_mb(k_r4+i-1) = damping*dbl_mb(k_r4+i-1)
     1                           / ( shift + denom - omega )
                                      enddo
                                    enddo
                                  enddo
                                enddo
                              enddo
                            enddo
                          enddo
                        enddo
                        call add_hash_block(d_t4,dbl_mb(k_r4),size,
     1                    int_mb(k_t4_offset),((((((((p1b-noab-1)*nvab
     2                    +p2b-noab-1)*nvab+p3b-noab-1)*nvab
     3                    +p4b-noab-1)*noab+h5b-1)*noab
     4                    +h6b-1)*noab+h7b-1)*noab+h8b-1))
                        if (nodezero.and.util_print('tr4',print_debug))
     1                    then
                          call get_hash_block(d_t4,dbl_mb(k_r4),size,
     1                    int_mb(k_t4_offset),((((((((p1b-noab-1)*nvab
     2                    +p2b-noab-1)*nvab+p3b-noab-1)*nvab
     3                    +p4b-noab-1)*noab+h5b-1)*noab
     4                    +h6b-1)*noab+h7b-1)*noab+h8b-1))
                          call ma_print_compact
     1                    (dbl_mb(k_r4),size,1,'tr4')
                        endif
                        if (.not.ma_pop_stack(l_r4))
     1                    call errquit('tce_jacobi_ir4: MA problem',1,
     2                    MA_ERR)
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
        enddo
      enddo
      next = nxtask(-nprocs,1)
      call ga_sync()
      return
      end
