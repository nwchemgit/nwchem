      SUBROUTINE ccsd2_q(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                   d_v2,k_v2_offset,d_y2,k_y2_offset,
     2                   d_e,k_e_offset,energy)
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
      integer d_v2
      integer k_v2_offset
      integer d_y2
      integer k_y2_offset
      integer d_e
      integer k_e_offset
      integer t_h1b, t_h1
      integer t_h2b, t_h2
      integer t_h3b, t_h3
      integer t_h4b, t_h4
      integer t_p5b, t_p5
      integer t_p6b, t_p6
      integer t_p7b, t_p7
      integer t_p8b, t_p8
      integer k_right,l_right
      integer k_left,l_left
      integer k_left_sorted,l_left_sorted
      integer size,i
      integer g_energy
      integer nxtask
      integer next
      integer nprocs
      integer count
      integer d_i1_1,d_i1_2
      integer k_i1_offset_1,l_i1_offset_1
      integer k_i1_offset_2,l_i1_offset_2
      double precision energy
      double precision factor
      external nxtask
c
c     Caution! k_right & k_left are not even allocated yet
c     but they will not be used.
c
      call ccsd2_q_right(dbl_mb(k_right),d_i1_1,d_i1_2,
     1  d_t1,d_t2,d_v2,k_i1_offset_1,k_i1_offset_2,
     2  k_t1_offset,k_t2_offset,k_v2_offset,l_i1_offset_1,
     3  l_i1_offset_2,t_h1b,t_h2b,t_h3b,t_h4b,
     4  t_p5b,t_p6b,t_p7b,t_p8b,1)
      call ccsd2_q_left(dbl_mb(k_left),d_v2,d_y2,k_v2_offset,
     1  k_y2_offset,t_h1b,t_h2b,t_h3b,t_h4b,
     2  t_p5b,t_p6b,t_p7b,t_p8b,1)
c
c     Accumulate
c
      if (.not.ga_create(mt_dbl,1,1,'perturbative',1,1,g_energy))
     1  call errquit('ccsd_q: GA problem',0,GA_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = nxtask(nprocs,1)
      energy=0.0d0
      do t_p5b = noab+1,noab+nvab
       do t_p6b = t_p5b,noab+nvab
        do t_p7b = t_p6b,noab+nvab
         do t_p8b = t_p7b,noab+nvab
          do t_h1b = 1,noab
           do t_h2b = t_h1b,noab
            do t_h3b = t_h2b,noab
             do t_h4b = t_h3b,noab
              if (next.eq.count) then
              if (int_mb(k_spin+t_p5b-1)
     1           +int_mb(k_spin+t_p6b-1)
     2           +int_mb(k_spin+t_p7b-1)
     3           +int_mb(k_spin+t_p8b-1)
     4        .eq.int_mb(k_spin+t_h1b-1)
     5           +int_mb(k_spin+t_h2b-1)
     6           +int_mb(k_spin+t_h3b-1)
     7           +int_mb(k_spin+t_h4b-1)) then
              if ((.not.restricted).or.
     1           (int_mb(k_spin+t_p5b-1)
     2           +int_mb(k_spin+t_p6b-1)
     3           +int_mb(k_spin+t_p7b-1)
     4           +int_mb(k_spin+t_p8b-1)
     5           +int_mb(k_spin+t_h1b-1)
     6           +int_mb(k_spin+t_h2b-1)
     7           +int_mb(k_spin+t_h3b-1)
     8           +int_mb(k_spin+t_h4b-1).le.12)) then
              if (ieor(int_mb(k_sym+t_p5b-1),
     1            ieor(int_mb(k_sym+t_p6b-1),
     2            ieor(int_mb(k_sym+t_p7b-1),
     3            ieor(int_mb(k_sym+t_p8b-1),
     4            ieor(int_mb(k_sym+t_h1b-1),
     5            ieor(int_mb(k_sym+t_h2b-1),
     6            ieor(int_mb(k_sym+t_h3b-1),
     7                 int_mb(k_sym+t_h4b-1)))))))).eq.0) then
              size = int_mb(k_range+t_p5b-1)
     1             * int_mb(k_range+t_p6b-1)
     2             * int_mb(k_range+t_p7b-1)
     3             * int_mb(k_range+t_p8b-1)
     4             * int_mb(k_range+t_h1b-1)
     5             * int_mb(k_range+t_h2b-1)
     6             * int_mb(k_range+t_h3b-1)
     7             * int_mb(k_range+t_h4b-1)
              if (.not.MA_PUSH_GET(mt_dbl,size,'right moment 2,3',
     1          l_right,k_right)) call errquit('ccsd_q',3,MA_ERR)
              if (.not.MA_PUSH_GET(mt_dbl,size,'left moment 2,3',
     1          l_left,k_left)) call errquit('ccsd_q',3,MA_ERR)
              if (.not.MA_PUSH_GET(mt_dbl,size,'left moment 2,3',
     1          l_left_sorted,k_left_sorted)) call errquit('ccsd_q',3,
     2          MA_ERR)
              do i = 1, size
               dbl_mb(k_right+i-1) = 0.0d0
               dbl_mb(k_left+i-1) = 0.0d0
              enddo
              call ccsd2_q_right(dbl_mb(k_right),d_i1_1,d_i1_2,
     1          d_t1,d_t2,d_v2,k_i1_offset_1,k_i1_offset_2,
     2          k_t1_offset,k_t2_offset,k_v2_offset,l_i1_offset_1,
     3          l_i1_offset_2,t_h1b,t_h2b,t_h3b,t_h4b,
     4          t_p5b,t_p6b,t_p7b,t_p8b,2)
              call ccsd2_q_left(dbl_mb(k_left),d_v2,d_y2,k_v2_offset,
     1          k_y2_offset,t_h1b,t_h2b,t_h3b,t_h4b,
     2          t_p5b,t_p6b,t_p7b,t_p8b,2)
              call tce_sort_8(dbl_mb(k_left),dbl_mb(k_left_sorted),
     1          int_mb(k_range+t_h1b-1),int_mb(k_range+t_h2b-1),
     2          int_mb(k_range+t_h3b-1),int_mb(k_range+t_h4b-1),
     3          int_mb(k_range+t_p5b-1),int_mb(k_range+t_p6b-1),
     4          int_mb(k_range+t_p7b-1),int_mb(k_range+t_p8b-1),
     5          5,6,7,8,1,2,3,4,1.0d0)
              if ((restricted).and.
     1           (int_mb(k_spin+t_p5b-1)
     2           +int_mb(k_spin+t_p6b-1)
     3           +int_mb(k_spin+t_p7b-1)
     4           +int_mb(k_spin+t_p8b-1)
     5           +int_mb(k_spin+t_h1b-1)
     6           +int_mb(k_spin+t_h2b-1)
     7           +int_mb(k_spin+t_h3b-1)
     8           +int_mb(k_spin+t_h4b-1).ne.12)) then
                factor = 2.0d0
              else
                factor = 1.0d0
              endif
              if ((t_p5b.eq.t_p6b).and.(t_p6b.eq.t_p7b)
     1                            .and.(t_p7b.eq.t_p8b)) then
                factor = factor / 24.0d0
              else if ((t_p5b.eq.t_p6b).and.(t_p6b.eq.t_p7b)) then
                factor = factor / 6.0d0
              else if ((t_p6b.eq.t_p7b).and.(t_p7b.eq.t_p8b)) then
                factor = factor / 6.0d0
              else if ((t_p5b.eq.t_p6b).and.(t_p7b.eq.t_p8b)) then
                factor = factor / 4.0d0
              else if (t_p5b.eq.t_p6b) then
                factor = factor / 2.0d0
              else if (t_p6b.eq.t_p7b) then
                factor = factor / 2.0d0
              else if (t_p7b.eq.t_p8b) then
                factor = factor / 2.0d0
              endif
              if ((t_h1b.eq.t_h2b).and.(t_h2b.eq.t_h3b)
     1                            .and.(t_h3b.eq.t_h4b)) then
                factor = factor / 24.0d0
              else if ((t_h1b.eq.t_h2b).and.(t_h2b.eq.t_h3b)) then
                factor = factor / 6.0d0
              else if ((t_h2b.eq.t_h3b).and.(t_h3b.eq.t_h4b)) then
                factor = factor / 6.0d0
              else if ((t_h1b.eq.t_h2b).and.(t_h3b.eq.t_h4b)) then
                factor = factor / 4.0d0
              else if (t_h1b.eq.t_h2b) then
                factor = factor / 2.0d0
              else if (t_h2b.eq.t_h3b) then
                factor = factor / 2.0d0
              else if (t_h3b.eq.t_h4b) then
                factor = factor / 2.0d0
              endif
              i = 0
              do t_p5 = 1, int_mb(k_range+t_p5b-1)
               do t_p6 = 1, int_mb(k_range+t_p6b-1)
                do t_p7 = 1, int_mb(k_range+t_p7b-1)
                 do t_p8 = 1, int_mb(k_range+t_p8b-1)
                  do t_h1 = 1, int_mb(k_range+t_h1b-1)
                   do t_h2 = 1, int_mb(k_range+t_h2b-1)
                    do t_h3 = 1, int_mb(k_range+t_h3b-1)
                     do t_h4 = 1, int_mb(k_range+t_h4b-1)
                    i = i + 1
                    energy = energy + factor * dbl_mb(k_right+i-1)
     1                                       * dbl_mb(k_left_sorted+i-1)
     2        / (-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     3           -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     4           -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p7b-1)+t_p7-1)
     5           -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p8b-1)+t_p8-1)
     6           +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     7           +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     8           +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1)
     9           +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h4b-1)+t_h4-1))
                     enddo
                    enddo
                   enddo
                  enddo
                 enddo
                enddo
               enddo
              enddo
              if (.not.MA_POP_STACK(l_left_sorted)) 
     1          call errquit('ccsd_q',6,MA_ERR)
              if (.not.MA_POP_STACK(l_left)) 
     1          call errquit('ccsd_q',6,MA_ERR)
              if (.not.MA_POP_STACK(l_right)) 
     1          call errquit('ccsd_q',6,MA_ERR)
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
      call ccsd2_q_left(dbl_mb(k_left),d_v2,d_y2,k_v2_offset,
     1  k_y2_offset,t_h1b,t_h2b,t_h3b,t_h4b,
     2  t_p5b,t_p6b,t_p7b,t_p8b,3)
      call ccsd2_q_right(dbl_mb(k_right),d_i1_1,d_i1_2,
     1  d_t1,d_t2,d_v2,k_i1_offset_1,k_i1_offset_2,
     2  k_t1_offset,k_t2_offset,k_v2_offset,l_i1_offset_1,
     3  l_i1_offset_2,t_h1b,t_h2b,t_h3b,t_h4b,
     4  t_p5b,t_p6b,t_p7b,t_p8b,3)
      call ga_zero(g_energy)
      call ga_acc(g_energy,1,1,1,1,energy,1,1.0d0)
      call ga_sync()
      call ga_get(g_energy,1,1,1,1,energy,1)
      if (.not.ga_destroy(g_energy))
     1  call errquit('ccsd_q: GA problem',1,GA_ERR)
      return
      end
