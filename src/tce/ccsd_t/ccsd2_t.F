      SUBROUTINE ccsd2_t(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                   d_f1,k_f1_offset,d_v2,k_v2_offset,
     2                   d_y1,k_y1_offset,d_y2,k_y2_offset,
     3                   d_e,k_e_offset,energy,size_t1)
C
C     $Id$
C
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer d_t1
      integer k_t1_offset
      integer d_t2
      integer k_t2_offset
      integer d_f1
      integer k_f1_offset
      integer d_v2
      integer k_v2_offset
      integer d_y1
      integer k_y1_offset
      integer d_y2
      integer k_y2_offset
      integer d_e
      integer k_e_offset
      integer t_h1b, t_h1
      integer t_h2b, t_h2
      integer t_h3b, t_h3
      integer t_p4b, t_p4
      integer t_p5b, t_p5
      integer t_p6b, t_p6
      integer k_right,l_right
      integer k_left,l_left
      integer k_left_sorted,l_left_sorted
      integer size,i
      integer g_energy
      integer nxtask
      integer next
      integer nprocs
      integer count
      integer d_i1_1,d_i1_2,d_i1_3,d_i1_4,d_i1_5
      integer k_i1_offset_1,l_i1_offset_1
      integer k_i1_offset_2,l_i1_offset_2
      integer k_i1_offset_3,l_i1_offset_3
      integer k_i1_offset_4,l_i1_offset_4
      integer k_i1_offset_5,l_i1_offset_5
c - T1/X1 LOCALIZATION -------------------
      integer l_t1_local,k_t1_local
      integer size_t1
c ---------------------------------------
      double precision energy
      double precision factor
      external nxtask
c
c - T1/X1 LOCALIZATION ----------
c    opening l_t1_local and l_x1_local
        if (.not.MA_PUSH_GET(mt_dbl,size_t1,'t1_local',
     1      l_t1_local,k_t1_local))
     1      call errquit('t1_local',1,MA_ERR)
        call ma_zero(dbl_mb(k_t1_local),size_t1)
c    copy d_t1 ==> l_t1_local
        call ga_get(d_t1,1,size_t1,1,1,dbl_mb(k_t1_local),1)
c -------------------------------
c
c
c     Caution! k_right & k_left are not even allocated yet
c     but they won't be used.
c
      call cr_ccsd_t_N(dbl_mb(k_right),d_f1,d_i1_1,d_i1_2,
     1  k_t1_local,d_t2,d_v2,k_f1_offset,k_i1_offset_1,k_i1_offset_2,
     2  k_t1_offset,k_t2_offset,k_v2_offset,l_i1_offset_1,
     3  l_i1_offset_2,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,1)
      call ccsd2_t_left(dbl_mb(k_left),d_f1,d_i1_3,d_i1_4,d_i1_5,
     1  d_t1,d_v2,d_y1,d_y2,k_f1_offset,k_i1_offset_3,
     2  k_i1_offset_4,k_i1_offset_5,k_t1_offset,k_v2_offset,
     3  k_y1_offset,k_y2_offset,l_i1_offset_3,l_i1_offset_4,
     4  l_i1_offset_5,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,1)
c
c     Get the numerator
c
      if (.not.ga_create(mt_dbl,1,1,'perturbative',1,1,g_energy))
     1  call errquit('ccsd_t: GA problem',0,GA_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = nxtask(nprocs,1)
      energy=0.0d0
      do t_p4b = noab+1,noab+nvab
       do t_p5b = t_p4b,noab+nvab
        do t_p6b = t_p5b,noab+nvab
         do t_h1b = 1,noab
          do t_h2b = t_h1b,noab
           do t_h3b = t_h2b,noab
            if (int_mb(k_spin+t_p4b-1)
     1         +int_mb(k_spin+t_p5b-1)
     2         +int_mb(k_spin+t_p6b-1)
     3      .eq.int_mb(k_spin+t_h1b-1)
     4         +int_mb(k_spin+t_h2b-1)
     5         +int_mb(k_spin+t_h3b-1)) then
            if ((.not.restricted).or.
     1         (int_mb(k_spin+t_p4b-1)
     1         +int_mb(k_spin+t_p5b-1)
     2         +int_mb(k_spin+t_p6b-1)
     3         +int_mb(k_spin+t_h1b-1)
     4         +int_mb(k_spin+t_h2b-1)
     5         +int_mb(k_spin+t_h3b-1).le.8)) then
            if (ieor(int_mb(k_sym+t_p4b-1),
     1          ieor(int_mb(k_sym+t_p5b-1),
     2          ieor(int_mb(k_sym+t_p6b-1),
     3          ieor(int_mb(k_sym+t_h1b-1),
     4          ieor(int_mb(k_sym+t_h2b-1),
     5               int_mb(k_sym+t_h3b-1)))))).eq.0) then
c
            if (next.eq.count) then
c
            size = int_mb(k_range+t_p4b-1)
     1           * int_mb(k_range+t_p5b-1)
     2           * int_mb(k_range+t_p6b-1)
     3           * int_mb(k_range+t_h1b-1)
     4           * int_mb(k_range+t_h2b-1)
     5           * int_mb(k_range+t_h3b-1)
            if (.not.MA_PUSH_GET(mt_dbl,size,'right moment 2,3',
     1        l_right,k_right)) call errquit('ccsd_t',3,MA_ERR)
            if (.not.MA_PUSH_GET(mt_dbl,size,'left moment 2,3',
     1        l_left,k_left)) call errquit('ccsd_t',3,MA_ERR)
            if (.not.MA_PUSH_GET(mt_dbl,size,'left moment 2,3',
     1        l_left_sorted,k_left_sorted)) call errquit('ccsd_t',3,
     2        MA_ERR)
            do i = 1, size
             dbl_mb(k_right+i-1) = 0.0d0
             dbl_mb(k_left+i-1) = 0.0d0
            enddo
            call cr_ccsd_t_N(dbl_mb(k_right),d_f1,d_i1_1,d_i1_2,
     1        k_t1_local,d_t2,
     &        d_v2,k_f1_offset,k_i1_offset_1,k_i1_offset_2,
     2        k_t1_offset,k_t2_offset,k_v2_offset,l_i1_offset_1,
     3        l_i1_offset_2,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
            call ccsd2_t_left(dbl_mb(k_left),d_f1,d_i1_3,d_i1_4,d_i1_5,
     1        d_t1,d_v2,d_y1,d_y2,k_f1_offset,k_i1_offset_3,
     2        k_i1_offset_4,k_i1_offset_5,k_t1_offset,k_v2_offset,
     3        k_y1_offset,k_y2_offset,l_i1_offset_3,l_i1_offset_4,
     4        l_i1_offset_5,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
            call tce_sort_6(dbl_mb(k_left),dbl_mb(k_left_sorted),
     1        int_mb(k_range+t_h1b-1),int_mb(k_range+t_h2b-1),
     2        int_mb(k_range+t_h3b-1),int_mb(k_range+t_p4b-1),
     3        int_mb(k_range+t_p5b-1),int_mb(k_range+t_p6b-1),
     4        4,5,6,1,2,3,1.0d0)
            if (restricted) then
              factor = 2.0d0
            else
              factor = 1.0d0
            endif
            if ((t_p4b.eq.t_p5b).and.(t_p5b.eq.t_p6b)) then
              factor = factor / 6.0d0
            else if ((t_p4b.eq.t_p5b).or.(t_p5b.eq.t_p6b)) then
              factor = factor / 2.0d0
            endif
            if ((t_h1b.eq.t_h2b).and.(t_h2b.eq.t_h3b)) then
              factor = factor / 6.0d0
            else if ((t_h1b.eq.t_h2b).or.(t_h2b.eq.t_h3b)) then
              factor = factor / 2.0d0
            endif
            i = 0
            do t_p4 = 1, int_mb(k_range+t_p4b-1)
             do t_p5 = 1, int_mb(k_range+t_p5b-1)
              do t_p6 = 1, int_mb(k_range+t_p6b-1)
               do t_h1 = 1, int_mb(k_range+t_h1b-1)
                do t_h2 = 1, int_mb(k_range+t_h2b-1)
                 do t_h3 = 1, int_mb(k_range+t_h3b-1)
                  i = i + 1
                  energy = energy + factor * dbl_mb(k_right+i-1)
     1                                     * dbl_mb(k_left_sorted+i-1)
     2      / (-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p4b-1)+t_p4-1)
     3         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     4         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     5         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     6         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     7         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1))
                 enddo
                enddo
               enddo
              enddo
             enddo
            enddo
            if (.not.MA_POP_STACK(l_left_sorted)) 
     1        call errquit('ccsd_t',6,MA_ERR)
            if (.not.MA_POP_STACK(l_left)) 
     1        call errquit('ccsd_t',6,MA_ERR)
            if (.not.MA_POP_STACK(l_right)) 
     1        call errquit('ccsd_t',6,MA_ERR)
            next = nxtask(nprocs,1)
            endif
            count = count + 1
            endif
            endif
            endif
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      next = nxtask(-nprocs,1)
      call ccsd2_t_left(dbl_mb(k_left),d_f1,d_i1_3,d_i1_4,d_i1_5,
     1  d_t1,d_v2,d_y1,d_y2,k_f1_offset,k_i1_offset_3,
     2  k_i1_offset_4,k_i1_offset_5,k_t1_offset,k_v2_offset,
     3  k_y1_offset,k_y2_offset,l_i1_offset_3,l_i1_offset_4,
     4  l_i1_offset_5,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,3)
      call cr_ccsd_t_N(dbl_mb(k_right),d_f1,d_i1_1,d_i1_2,
     1  k_t1_local,d_t2,d_v2,k_f1_offset,k_i1_offset_1,k_i1_offset_2,
     2  k_t1_offset,k_t2_offset,k_v2_offset,l_i1_offset_1,
     3  l_i1_offset_2,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,3)
      call ga_zero(g_energy)
      call ga_acc(g_energy,1,1,1,1,energy,1,1.0d0)
      call ga_sync()
      call ga_get(g_energy,1,1,1,1,energy,1)
      if (.not.ga_destroy(g_energy))
     1  call errquit('ccsd_t: GA problem',1,GA_ERR)
c - T1/X1 LOCALIZATION ------
         if(.not.MA_POP_STACK(l_t1_local))
     &      call errquit('l_t1_local',4,MA_ERR)
c ---------------------------
      return
      end
