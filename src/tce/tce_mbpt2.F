      subroutine tce_mbpt2(d_mo2e,k_2e_offset,
     1                     d_t2,k_t2_offset,mbpt2)
c
c $Id: tce_mbpt2.F,v 1.1 2002-10-18 23:42:08 sohirata Exp $
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
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
      integer nxtval
      external nxtval
      double precision ddot
      external ddot
c
c     ================
c     Loop over blocks
c     ================
c
      nprocs = ga_nnodes()
      count = 0
      next = nxtval(nprocs)
      if (.not.ga_create(mt_dbl,1,1,'mbpt2',1,1,g_mbpt2))
     1  call errquit('tce_mbpt2: GA problem',0)
      mbpt2 = 0.0d0
      call ga_put(g_mbpt2,1,1,1,1,mbpt2,1)
      do p1b = noab+1,noab+nvab
        do p2b = p1b,noab+nvab
          do h3b = 1,noab
            do h4b = h3b,noab
              if (next.eq.count) then
                if (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) 
     1            .eq. int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) then
                if (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),
     1            ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1))))
     2            .eq. 0) then
                    factor = 0.25d0
                    if (p2b .gt. p1b) factor = factor * 2.0d0
                    if (h4b .gt. h3b) factor = factor * 2.0d0
                    size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
     1                   * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
                    if (.not.ma_push_get(mt_dbl,size,'v2',l_v2,k_v2))
     1                call errquit('tce_mbpt2: MA problem',0)
                    call get_block(d_mo2e,dbl_mb(k_v2),size,
     1                int_mb(k_2e_offset+(((p1b-1)*(noab+nvab)+p2b-1)
     2                *(noab+nvab)+h3b-1)*(noab+nvab)+h4b-1))
                    if (.not.ma_push_get(mt_dbl,size,'t2',l_t2,k_t2))
     1                call errquit('tce_mbpt2: MA problem',1)
                    call get_block(d_t2,dbl_mb(k_t2),size,
     1                int_mb(k_t2_offset+(((p1b-noab-1)*nvab+p2b-noab-1)
     2                *noab+h3b-1)*noab+h4b-1))
                    mbpt2 = ddot(size,dbl_mb(k_t2),1,dbl_mb(k_v2),1)
                    call ga_acc(g_mbpt2,1,1,1,1,mbpt2,1,factor)
                    if (.not.ma_pop_stack(l_t2))
     1                call errquit('tce_mbpt2: MA problem',2)
                    if (.not.ma_pop_stack(l_v2))
     1                call errquit('tce_mbpt2: MA problem',3)
                  endif
                endif
                next = nxtval(nprocs)
              endif
              count = count + 1
            enddo
          enddo
        enddo
      enddo
      call ga_get(g_mbpt2,1,1,1,1,mbpt2,1)
      if (.not.ga_destroy(g_mbpt2))
     1  call errquit('tce_mbpt2: GA problem',1)
      next = nxtval(-nprocs)
      call ga_sync()
      return
      end
