      subroutine tce_guess_t2(d_mo2e,k_2e_offset,
     1                        d_t2,k_t2_offset)
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
      integer d_mo2e
      integer d_t2
      integer p1b
      integer p2b
      integer h3b
      integer h4b
      integer p1
      integer p2
      integer h3
      integer h4
      integer k_2e_offset
      integer k_t2_offset
      integer size
      integer l_t2,k_t2
      integer i
      integer nprocs
      integer count
      integer next
cc      integer nxtask
cc      external nxtask
      INTEGER NXTASK
      EXTERNAL NXTASK
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
                if (util_print('t2',print_debug))
     1            write(LuOut,9200) p1b,p2b,h3b,h4b,ga_nodeid()
                if (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1) 
     1            .eq. int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) then
                if ((.not.restricted).or.
     1            (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+
     2            int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) then
                if (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),
     1            ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1))))
     2            .eq. 0) then
                 size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
     1                * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
                 if (.not.ma_push_get(mt_dbl,size,'t2',l_t2,k_t2))
     1             call errquit('tce_guess_t2: MA problem',0,MA_ERR)
c
        if(.not.intorb) then
                 call get_hash_block(d_mo2e,dbl_mb(k_t2),size,
     1             int_mb(k_2e_offset),((((p1b-1)*(noab+nvab)+p2b-1)
     2             *(noab+nvab)+h3b-1)*(noab+nvab)+h4b-1))
        else
                 call get_hash_block_i(d_mo2e,dbl_mb(k_t2),size,
     1             int_mb(k_2e_offset),((((p1b-1)*(noab+nvab)+p2b-1)
     2             *(noab+nvab)+h3b-1)*(noab+nvab)+h4b-1),
     3             h4b,h3b,p2b,p1b)
        end if
                 i = 0
                 do p1 = 1,int_mb(k_range+p1b-1)
                   do p2 = 1,int_mb(k_range+p2b-1)
                     do h3 = 1,int_mb(k_range+h3b-1)
                       do h4 = 1,int_mb(k_range+h4b-1)
                         i = i + 1
                         dbl_mb(k_t2+i-1) = dbl_mb(k_t2+i-1)
     1          / (-dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1)
     2             -dbl_mb(k_evl_sorted+int_mb(k_offset+p2b-1)+p2-1)
     3             +dbl_mb(k_evl_sorted+int_mb(k_offset+h3b-1)+h3-1)
     4             +dbl_mb(k_evl_sorted+int_mb(k_offset+h4b-1)+h4-1))
                       enddo
                     enddo
                   enddo
                 enddo
                 if (util_print('t2',print_debug))
     1             call ma_print_compact(dbl_mb(k_t2),size,1,'t2')
                 call put_hash_block(d_t2,dbl_mb(k_t2),size,
     1             int_mb(k_t2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2             *noab+h3b-1)*noab+h4b-1))
                 if (.not.ma_pop_stack(l_t2))
     1             call errquit('tce_guess_t2: MA problem',1,MA_ERR)
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
 9200 format(/,4i3,' has been assigned to node ',i3)
      return
      end


      subroutine tce_guess_t2_zero(d_t2,k_t2_offset)
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
      integer l_t2,k_t2
      integer i
      integer nprocs
      integer count
      integer next
      integer nxtask
      external nxtask

      nprocs = ga_nnodes()
      count = 0
      next = NXTASK(nprocs,1)
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
                 size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
     1                * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
                 if (.not.ma_push_get(mt_dbl,size,'t2',l_t2,k_t2))
     1             call errquit('tce_guess_t2: MA problem',0,MA_ERR)

                 i = 0
                 do p1 = 1,int_mb(k_range+p1b-1)
                   do p2 = 1,int_mb(k_range+p2b-1)
                     do h3 = 1,int_mb(k_range+h3b-1)
                       do h4 = 1,int_mb(k_range+h4b-1)
                         i = i + 1
                         dbl_mb(k_t2+i-1) = 0.0d0
                       enddo
                     enddo
                   enddo
                 enddo

                 call put_hash_block(d_t2,dbl_mb(k_t2),size,
     1             int_mb(k_t2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2             *noab+h3b-1)*noab+h4b-1))
                 if (.not.ma_pop_stack(l_t2))
     1             call errquit('tce_guess_t2: MA problem',1,MA_ERR)
                endif
                endif
                endif
                next = NXTASK(nprocs,1)
              endif
              count = count + 1
            enddo
          enddo
        enddo
      enddo
      next = NXTASK(-nprocs,1)
      call ga_sync()
      return
      end

