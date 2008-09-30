      subroutine tce_make_denom2(d_denom2,k_denom2_offset,irrep,
     1                           denom_power,omega,shift)
c
c $Id: tce_make_denom2.F,v 1.2 2008-09-30 18:43:17 jhammond Exp $
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
      integer d_denom2
      integer p1b
      integer p2b
      integer h3b
      integer h4b
      integer p1
      integer p2
      integer h3
      integer h4
      integer k_2e_offset
      integer k_denom2_offset
      integer size
      integer l_denom2,k_denom2
      integer i,irrep
      integer nprocs
      integer count
      integer next
      INTEGER NXTASK
      EXTERNAL NXTASK
      integer denom_power
      double precision shift,omega
c
c     ================
c     Loop over blocks
c     ================
c
      nprocs = ga_nnodes()
      count = 0
      next = NXTASK(nprocs, 1)
      do p1b = noab+1,noab+nvab
        do p2b = p1b,noab+nvab
          do h3b = 1,noab
            do h4b = h3b,noab
              if (next.eq.count) then
                if (util_print('denom2',print_debug))
     1            write(LuOut,9200) p1b,p2b,h3b,h4b,ga_nodeid()
                if (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) 
     1            .eq. int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) then
                if ((.not.restricted).or.
     1            (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+
     2            int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) then
                if (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),
     1            ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1))))
     2            .eq. irrep) then
                 size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
     1                * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
                 if (.not.ma_push_get(mt_dbl,size,'denom2',l_denom2,
     1                                                     k_denom2))
     2             call errquit('tce_make_denom2: MA problem',0,MA_ERR)
                 i = 0
                 do p1 = 1,int_mb(k_range+p1b-1)
                   do p2 = 1,int_mb(k_range+p2b-1)
                     do h3 = 1,int_mb(k_range+h3b-1)
                       do h4 = 1,int_mb(k_range+h4b-1)
                         i = i + 1
c                         dbl_mb(k_denom2+i-1) = dbl_mb(k_denom2+i-1)
c     1          / (-dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1)
c     2             -dbl_mb(k_evl_sorted+int_mb(k_offset+p2b-1)+p2-1)
c     3             +dbl_mb(k_evl_sorted+int_mb(k_offset+h3b-1)+h3-1)
c     4             +dbl_mb(k_evl_sorted+int_mb(k_offset+h4b-1)+h4-1))
                       dbl_mb(k_denom2+i-1) =
     1         (( -dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1)
     2            -dbl_mb(k_evl_sorted+int_mb(k_offset+p2b-1)+p2-1)
     3            +dbl_mb(k_evl_sorted+int_mb(k_offset+h3b-1)+h3-1)
     4            +dbl_mb(k_evl_sorted+int_mb(k_offset+h4b-1)+h4-1) )
     5           **denom_power-omega**denom_power+shift)
                       enddo
                     enddo
                   enddo
                 enddo
                 if (util_print('denom2',print_debug))
     1             call ma_print_compact(dbl_mb(k_denom2),size,1,
     2             'denom2')
                 call put_hash_block(d_denom2,dbl_mb(k_denom2),size,
     1             int_mb(k_denom2_offset),((((p1b-noab-1)*nvab+
     2             p2b-noab-1)*noab+h3b-1)*noab+h4b-1))
                 if (.not.ma_pop_stack(l_denom2))
     1             call errquit('tce_make_denom2: MA problem',1,MA_ERR)
                endif
                endif
                endif
                next = NXTASK(nprocs, 1)
              endif
              count = count + 1
            enddo
          enddo
        enddo
      enddo
      next = NXTASK(-nprocs, 1)
      call ga_sync()
 9200 format(/,4i3,' has been assigned to node ',i3)
      return
      end
