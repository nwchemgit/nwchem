      SUBROUTINE lr_ccsd_tq_1(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                   d_v2,k_v2_offset,
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
c
      integer l_hole,k_hole
      integer l_2hole,k_2hole
      integer hole_p_1
      integer hole_p_2
      double precision diffeorb,xlocden
c
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
      k_right = 0
      k_left = 0
c
c forming the nu(i) and mu(i,j) arrays (i,j refer to occupied indices)
c no offset for nu(i) and mu(i,j) is needed
      hole_p_1=nocc(1)-nfc(1)+nocc(ipol)-nfc(ipol)
      if (.not.ma_push_get(mt_dbl,hole_p_1,
     1  "one hole",l_hole,k_hole)) call errquit('lr_ccsd_t',1,MA_ERR)
      hole_p_2=hole_p_1*hole_p_1
      if (.not.ma_push_get(mt_dbl,hole_p_2,"two holes",
     2    l_2hole,k_2hole)) call errquit('lr_ccsd_t',2,MA_ERR)
c
c here goes the procedure to calculate the nu(i) and mu(i,j)
        do i=1,hole_p_1
         dbl_mb(k_hole+i-1)=0.0d0
        enddo
        do i=1,hole_p_2
         dbl_mb(k_2hole+i-1)=0.0d0
        enddo
        call tce_nu1(k_hole,d_t1,k_t1_offset)
        call tce_mu2(k_2hole,hole_p_1,d_t2,k_t2_offset)
c
c
c     Caution! k_right & k_left are not even allocated yet
c     but they won't be used.
c
c
c
c
      call ccsd2_q_right(dbl_mb(k_right),d_i1_1,d_i1_2,
     1  d_t1,d_t2,d_v2,k_i1_offset_1,k_i1_offset_2,
     2  k_t1_offset,k_t2_offset,k_v2_offset,l_i1_offset_1,
     3  l_i1_offset_2,t_h1b,t_h2b,t_h3b,t_h4b,
     4  t_p5b,t_p6b,t_p7b,t_p8b,1)
      call lr_t2_t2(dbl_mb(k_left),d_t2,k_t2_offset,
     &  t_h1b,t_h2b,t_h3b,t_h4b,
     &  t_p5b,t_p6b,t_p7b,t_p8b,1)
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
     1          l_right,k_right)) call errquit('lr_ccsd_tq',3,MA_ERR)
              if (.not.MA_PUSH_GET(mt_dbl,size,'T2 T2',
     1          l_left,k_left)) call errquit('lr-ccsd_tq',3,MA_ERR)
              do i = 1, size
               dbl_mb(k_right+i-1) = 0.0d0
               dbl_mb(k_left+i-1) = 0.0d0
              enddo
              call ccsd2_q_right(dbl_mb(k_right),d_i1_1,d_i1_2,
     1          d_t1,d_t2,d_v2,k_i1_offset_1,k_i1_offset_2,
     2          k_t1_offset,k_t2_offset,k_v2_offset,l_i1_offset_1,
     3          l_i1_offset_2,t_h1b,t_h2b,t_h3b,t_h4b,
     4          t_p5b,t_p6b,t_p7b,t_p8b,2)
              call lr_t2_t2(dbl_mb(k_left),d_t2,k_t2_offset,
     &          t_h1b,t_h2b,t_h3b,t_h4b,
     &          t_p5b,t_p6b,t_p7b,t_p8b,2)
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
c
        diffeorb=-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     3           -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     4           -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p7b-1)+t_p7-1)
     5           -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p8b-1)+t_p8-1)
     6           +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     7           +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     8           +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1)
     9           +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h4b-1)+t_h4-1)
c
        xlocden=1.0d0+dbl_mb(k_hole+int_mb(k_offset+t_h1b-1)+t_h1-1)
     &               +dbl_mb(k_hole+int_mb(k_offset+t_h2b-1)+t_h2-1)
     &               +dbl_mb(k_hole+int_mb(k_offset+t_h3b-1)+t_h3-1)
     &               +dbl_mb(k_hole+int_mb(k_offset+t_h4b-1)+t_h4-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h1b-1)+t_h1-1)
     &            +int_mb(k_offset+t_h2b-1)+t_h2-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h1b-1)+t_h1-1)
     &            +int_mb(k_offset+t_h3b-1)+t_h3-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h2b-1)+t_h2-1)
     &            +int_mb(k_offset+t_h3b-1)+t_h3-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h1b-1)+t_h1-1)
     &            +int_mb(k_offset+t_h4b-1)+t_h4-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h2b-1)+t_h2-1)
     &            +int_mb(k_offset+t_h4b-1)+t_h4-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h3b-1)+t_h3-1)
     &            +int_mb(k_offset+t_h4b-1)+t_h4-1)
c
                    energy = energy + factor * dbl_mb(k_right+i-1)
     1                                       * dbl_mb(k_left+i-1)/
     2                                         (xlocden)
c
                    energy = energy + factor * dbl_mb(k_right+i-1)
     1                                       * dbl_mb(k_right+i-1)/
     2                                         (diffeorb*xlocden)
c
                     enddo
                    enddo
                   enddo
                  enddo
                 enddo
                enddo
               enddo
              enddo
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
      call lr_t2_t2(dbl_mb(k_left),d_t2,k_t2_offset,
     &  t_h1b,t_h2b,t_h3b,t_h4b,
     &  t_p5b,t_p6b,t_p7b,t_p8b,3)
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
c
      if (.not.ma_pop_stack(l_2hole))
     1  call errquit("tce_energy: MA problem",19,MA_ERR)
      if (.not.ma_pop_stack(l_hole))
     1  call errquit("tce_energy: MA problem",19,MA_ERR)
c
      return
      end
c
c
c
c FULL VERSION 
c
c
c
      SUBROUTINE lr_ccsd_tq1p(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                   d_v2,k_v2_offset,
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
c local denominators
      integer l_hole,k_hole
      integer l_2hole,k_2hole
      integer hole_p_1
      integer hole_p_2
      double precision diffeorb,xlocden
c needed in lr_12prod_q call
      integer d_ii1_1,d_ii1_3
      integer k_ii1_offset_1,k_ii1_offset_3
      integer l_ii1_offset_1,l_ii1_offset_3
c
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
c
c forming the nu(i) and mu(i,j) arrays (i,j refer to occupied indices)
c no offset for nu(i) and mu(i,j) is needed
      hole_p_1=nocc(1)-nfc(1)+nocc(ipol)-nfc(ipol)
      if (.not.ma_push_get(mt_dbl,hole_p_1,
     1  "one hole",l_hole,k_hole)) call errquit('lr_ccsd_t',1,MA_ERR)
      hole_p_2=hole_p_1*hole_p_1
      if (.not.ma_push_get(mt_dbl,hole_p_2,"two holes",
     2    l_2hole,k_2hole)) call errquit('lr_ccsd_t',2,MA_ERR)
c
c here goes the procedure to calculate the nu(i) and mu(i,j)
        do i=1,hole_p_1
         dbl_mb(k_hole+i-1)=0.0d0
        enddo
        do i=1,hole_p_2
         dbl_mb(k_2hole+i-1)=0.0d0
        enddo
        call tce_nu1(k_hole,d_t1,k_t1_offset)
        call tce_mu2(k_2hole,hole_p_1,d_t2,k_t2_offset)
c
c
c
c
c     Caution! k_right & k_left are not even allocated yet
c     but they won't be used.
c
c
c
c
      call ccsd2_q_right(dbl_mb(k_right),d_i1_1,d_i1_2,
     1  d_t1,d_t2,d_v2,k_i1_offset_1,k_i1_offset_2,
     2  k_t1_offset,k_t2_offset,k_v2_offset,l_i1_offset_1,
     3  l_i1_offset_2,t_h1b,t_h2b,t_h3b,t_h4b,
     4  t_p5b,t_p6b,t_p7b,t_p8b,1)
      call lr_12prod_q(dbl_mb(k_left),d_ii1_1,d_ii1_3,d_t1,
     &d_t2,k_ii1_offset_1,k_ii1_offset_3,k_t1_offset,
     &k_t2_offset,l_ii1_offset_1,l_ii1_offset_3,
     &t_h1b,t_h2b,t_h3b,t_h4b,t_p5b,t_p6b,t_p7b,t_p8b
     &,1)
ccc      call lr_t1t3_q(dbl_mb(k_right2),d_t1,d_m3loc,
ccc     &k_t1_offset,k_m3loc_offset,t_h1b,
ccc     &t_h2b,t_h3b,t_h4b,t_p5b,t_p6b,t_p7b,t_p8b,1)
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
     1          l_right,k_right)) call errquit('lr_ccsd_tq',3,MA_ERR)
              if (.not.MA_PUSH_GET(mt_dbl,size,'T2 T2',
     1          l_left,k_left)) call errquit('lr-ccsd_tq',3,MA_ERR)
              do i = 1, size
               dbl_mb(k_right+i-1) = 0.0d0
               dbl_mb(k_left+i-1) = 0.0d0
              enddo
              call ccsd2_q_right(dbl_mb(k_right),d_i1_1,d_i1_2,
     1          d_t1,d_t2,d_v2,k_i1_offset_1,k_i1_offset_2,
     2          k_t1_offset,k_t2_offset,k_v2_offset,l_i1_offset_1,
     3          l_i1_offset_2,t_h1b,t_h2b,t_h3b,t_h4b,
     4          t_p5b,t_p6b,t_p7b,t_p8b,2)
              call lr_12prod_q(dbl_mb(k_left),d_ii1_1,d_ii1_3,d_t1,
     &           d_t2,k_ii1_offset_1,k_ii1_offset_3,k_t1_offset,
     &           k_t2_offset,l_ii1_offset_1,l_ii1_offset_3,
     &           t_h1b,t_h2b,t_h3b,t_h4b,t_p5b,t_p6b,t_p7b,t_p8b
     &           ,2)
c              call lr_t1t3_q(dbl_mb(k_right2),d_t1,d_m3loc,
c     &           k_t1_offset,k_m3loc_offset,t_h1b,
c     &           t_h2b,t_h3b,t_h4b,t_p5b,t_p6b,t_p7b,t_p8b,2)
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
c
        diffeorb=-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     3           -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     4           -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p7b-1)+t_p7-1)
     5           -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p8b-1)+t_p8-1)
     6           +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     7           +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     8           +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1)
     9           +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h4b-1)+t_h4-1)
c
        xlocden=1.0d0+dbl_mb(k_hole+int_mb(k_offset+t_h1b-1)+t_h1-1)
     &               +dbl_mb(k_hole+int_mb(k_offset+t_h2b-1)+t_h2-1)
     &               +dbl_mb(k_hole+int_mb(k_offset+t_h3b-1)+t_h3-1)
     &               +dbl_mb(k_hole+int_mb(k_offset+t_h4b-1)+t_h4-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h1b-1)+t_h1-1)
     &            +int_mb(k_offset+t_h2b-1)+t_h2-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h1b-1)+t_h1-1)
     &            +int_mb(k_offset+t_h3b-1)+t_h3-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h2b-1)+t_h2-1)
     &            +int_mb(k_offset+t_h3b-1)+t_h3-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h1b-1)+t_h1-1)
     &            +int_mb(k_offset+t_h4b-1)+t_h4-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h2b-1)+t_h2-1)
     &            +int_mb(k_offset+t_h4b-1)+t_h4-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h3b-1)+t_h3-1)
     &            +int_mb(k_offset+t_h4b-1)+t_h4-1)
c <(all_12products)+*M4/local_den>
                    energy = energy + factor * dbl_mb(k_right+i-1)
     1                                       * dbl_mb(k_left+i-1)/
     2                                         (xlocden)
c <T4+*(M4/local_den)>
                    energy = energy + factor * dbl_mb(k_right+i-1)
     1                                       * dbl_mb(k_right+i-1)/
     2                                         (diffeorb*xlocden)
c
                     enddo
                    enddo
                   enddo
                  enddo
                 enddo
                enddo
               enddo
              enddo
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
      call lr_12prod_q(dbl_mb(k_left),d_ii1_1,d_ii1_3,d_t1,
     &   d_t2,k_ii1_offset_1,k_ii1_offset_3,k_t1_offset,
     &   k_t2_offset,l_ii1_offset_1,l_ii1_offset_3,
     &   t_h1b,t_h2b,t_h3b,t_h4b,t_p5b,t_p6b,t_p7b,t_p8b
     &   ,3)
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
c
      if (.not.ma_pop_stack(l_2hole))
     1  call errquit("tce_energy: MA problem",19,MA_ERR)
      if (.not.ma_pop_stack(l_hole))
     1  call errquit("tce_energy: MA problem",19,MA_ERR)
c
      return
      end
c
c
c
c
c  PRODUCT 1/2T2*T2
c
c
c
      SUBROUTINE lr_t2_t2(a_i0,d_t2,k_t2_offset,t_h1b,t_h2b,t_h3b,t_h4b,
     &t_p5b,t_p6b,t_p7b,t_p8b,toggle)
C     $Id$
C     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
C     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
C     i0 ( p5 p6 p7 p8 h1 h2 h3 h4 )_tt + = 1/2 * P( 36 ) * t ( p5 p6 h1 h2 )_t * t ( p7 p8 h3 h4 )_t
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
      INTEGER t_p5b
      INTEGER t_p6b
      INTEGER t_p7b
      INTEGER t_p8b
      INTEGER t_h1b
      INTEGER t_h2b
      INTEGER t_h3b
      INTEGER t_h4b
      INTEGER toggle
      INTEGER d_t2
      INTEGER k_t2_offset
      DOUBLE PRECISION a_i0(*)
      IF (toggle .eq. 2) CALL lr_t2_t2_1(d_t2,k_t2_offset,d_t2,k_t2_offs
     &et,a_i0,t_p5b,t_p6b,t_p7b,t_p8b,t_h1b,t_h2b,t_h3b,t_h4b)
      RETURN
      END
      SUBROUTINE lr_t2_t2_1(d_a,k_a_offset,d_b,k_b_offset,a_c,t_p5b,t_p6
     &b,t_p7b,t_p8b,t_h1b,t_h2b,t_h3b,t_h4b)
C     $Id$
C     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
C     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
C     i0 ( p5 p6 p7 p8 h1 h2 h3 h4 )_tt + = 1/2 * P( 36 ) * t ( p5 p6 h1 h2 )_t * t ( p7 p8 h3 h4 )_t
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
      INTEGER d_a
      INTEGER k_a_offset
      INTEGER d_b
      INTEGER k_b_offset
      INTEGER t_p5b
      INTEGER t_p6b
      INTEGER t_p7b
      INTEGER t_p8b
      INTEGER t_h1b
      INTEGER t_h2b
      INTEGER t_h3b
      INTEGER t_h4b
      INTEGER p5b
      INTEGER p6b
      INTEGER p7b
      INTEGER p8b
      INTEGER h1b
      INTEGER h2b
      INTEGER h3b
      INTEGER h4b
      INTEGER dimc
      INTEGER l_c_sort
      INTEGER k_c_sort
      INTEGER p5b_1
      INTEGER p6b_1
      INTEGER h1b_1
      INTEGER h2b_1
      INTEGER p7b_2
      INTEGER p8b_2
      INTEGER h3b_2
      INTEGER h4b_2
      INTEGER dim_common
      INTEGER dima_sort
      INTEGER dima
      INTEGER dimb_sort
      INTEGER dimb
      INTEGER l_a_sort
      INTEGER k_a_sort
      INTEGER l_a
      INTEGER k_a
      INTEGER l_b_sort
      INTEGER k_b_sort
      INTEGER l_b
      INTEGER k_b
      DOUBLE PRECISION a_c(*)
      LOGICAL skipped
      DO p5b = noab+1,noab+nvab
      DO p6b = p5b,noab+nvab
      DO p7b = noab+1,noab+nvab
      DO p8b = p7b,noab+nvab
      DO h1b = 1,noab
      DO h2b = h1b,noab
      DO h3b = 1,noab
      DO h4b = h3b,noab
      skipped = .true.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &4b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &4b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &4b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &4b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &4b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p8b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p8b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p8b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p8b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p8b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p8b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &4b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF (.not.skipped) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1
     &)+int_mb(k_spin+p7b-1)+int_mb(k_spin+p8b-1)+int_mb(k_spin+h1b-1)+i
     &nt_mb(k_spin+h2b-1)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.1
     &6)) THEN
      IF (int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1)+int_mb(k_spin+p7b-1)
     &+int_mb(k_spin+p8b-1) .eq. int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-
     &1)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),ieor(int_mb(k_sym+p6b-1),ieor(int_mb(
     &k_sym+p7b-1),ieor(int_mb(k_sym+p8b-1),ieor(int_mb(k_sym+h1b-1),ieo
     &r(int_mb(k_sym+h2b-1),ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1)
     &))))))) .eq. ieor(irrep_t,irrep_t)) THEN
      dimc = int_mb(k_range+p5b-1) * int_mb(k_range+p6b-1) * int_mb(k_ra
     &nge+p7b-1) * int_mb(k_range+p8b-1) * int_mb(k_range+h1b-1) * int_m
     &b(k_range+h2b-1) * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c_sort,k_c_sort)) CALL
     & ERRQUIT('lr_t2_t2_1',0,MA_ERR)
      CALL DFILL(dimc,0.0d0,dbl_mb(k_c_sort),1)
      IF (int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1) .eq. int_mb(k_spin+h
     &1b-1)+int_mb(k_spin+h2b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),ieor(int_mb(k_sym+p6b-1),ieor(int_mb(
     &k_sym+h1b-1),int_mb(k_sym+h2b-1)))) .eq. irrep_t) THEN
      CALL TCE_RESTRICTED_4(p5b,p6b,h1b,h2b,p5b_1,p6b_1,h1b_1,h2b_1)
      CALL TCE_RESTRICTED_4(p7b,p8b,h3b,h4b,p7b_2,p8b_2,h3b_2,h4b_2)
      dim_common = 1
      dima_sort = int_mb(k_range+p5b-1) * int_mb(k_range+p6b-1) * int_mb
     &(k_range+h1b-1) * int_mb(k_range+h2b-1)
      dima = dim_common * dima_sort
      dimb_sort = int_mb(k_range+p7b-1) * int_mb(k_range+p8b-1) * int_mb
     &(k_range+h3b-1) * int_mb(k_range+h4b-1)
      dimb = dim_common * dimb_sort
      IF ((dima .gt. 0) .and. (dimb .gt. 0)) THEN
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a_sort,k_a_sort)) CALL
     & ERRQUIT('lr_t2_t2_1',1,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a,k_a)) CALL ERRQUIT('
     &lr_t2_t2_1',2,MA_ERR)
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h2b_1
     & - 1 + noab * (h1b_1 - 1 + noab * (p6b_1 - noab - 1 + nvab * (p5b_
     &1 - noab - 1)))))
      CALL TCE_SORT_4(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p5b-1)
     &,int_mb(k_range+p6b-1),int_mb(k_range+h1b-1),int_mb(k_range+h2b-1)
     &,4,3,2,1,1.0d0)
      IF (.not.MA_POP_STACK(l_a)) CALL ERRQUIT('lr_t2_t2_1',3,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b_sort,k_b_sort)) CALL
     & ERRQUIT('lr_t2_t2_1',4,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b,k_b)) CALL ERRQUIT('
     &lr_t2_t2_1',5,MA_ERR)
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),(h4b_2
     & - 1 + noab * (h3b_2 - 1 + noab * (p8b_2 - noab - 1 + nvab * (p7b_
     &2 - noab - 1)))))
      CALL TCE_SORT_4(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+p7b-1)
     &,int_mb(k_range+p8b-1),int_mb(k_range+h3b-1),int_mb(k_range+h4b-1)
     &,4,3,2,1,1.0d0)
      IF (.not.MA_POP_STACK(l_b)) CALL ERRQUIT('lr_t2_t2_1',6,MA_ERR)
      CALL DGEMM('T','N',dima_sort,dimb_sort,dim_common,1.0d0,dbl_mb(k_a
     &_sort),dim_common,dbl_mb(k_b_sort),dim_common,1.0d0,dbl_mb(k_c_sor
     &t),dima_sort)
      IF (.not.MA_POP_STACK(l_b_sort)) CALL ERRQUIT('lr_t2_t2_1',7,MA_ER
     &R)
      IF (.not.MA_POP_STACK(l_a_sort)) CALL ERRQUIT('lr_t2_t2_1',8,MA_ER
     &R)
      END IF
      END IF
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,7,4,3,6,5,2,1,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,7,4,3,2,6,5,1,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,7,4,3,2,6,1,5,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,7,4,3,6,2,5,1,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,7,4,3,6,2,1,5,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &4b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,7,4,3,2,1,6,5,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,8,7,3,6,5,2,1,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,8,7,3,2,6,5,1,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,8,7,3,2,6,1,5,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,8,7,3,6,2,5,1,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,8,7,3,6,2,1,5,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &4b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,8,7,3,2,1,6,5,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,8,3,7,6,5,2,1,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,8,3,7,2,6,5,1,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,8,3,7,2,6,1,5,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,8,3,7,6,2,5,1,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,8,3,7,6,2,1,5,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &4b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,8,3,7,2,1,6,5,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,4,7,3,6,5,2,1,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,4,7,3,2,6,5,1,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,4,7,3,2,6,1,5,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,4,7,3,6,2,5,1,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,4,7,3,6,2,1,5,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &4b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,4,7,3,2,1,6,5,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,4,3,7,6,5,2,1,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,4,3,7,2,6,5,1,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,4,3,7,2,6,1,5,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,4,3,7,6,2,5,1,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,4,3,7,6,2,1,5,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &4b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,4,3,7,2,1,6,5,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p8b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,3,8,7,6,5,2,1,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p8b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,3,8,7,2,6,5,1,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p8b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,3,8,7,2,6,1,5,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p8b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,3,8,7,6,2,5,1,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p8b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,3,8,7,6,2,1,5,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p8b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &4b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,3,8,7,2,1,6,5,1.0d0/2.0d0)
      END IF
      IF (.not.MA_POP_STACK(l_c_sort)) CALL ERRQUIT('lr_t2_t2_1',9,MA_ER
     &R)
      END IF
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      RETURN
      END
c
c
c
c
c ALL PRODUCTS 1/2T2*T2+1/2T1*T1*T2+1/24T1*T1*T1*T1
c
c
c
      SUBROUTINE lr_12prod_q(a_i0,d_i1_1,d_i1_2,d_t1,d_t2,k_i1_offset_1,
     &k_i1_offset_2,k_t1_offset,k_t2_offset,l_i1_offset_1,l_i1_offset_2,
     &t_h1b,t_h2b,t_h3b,t_h4b,t_p5b,t_p6b,t_p7b,t_p8b,toggle)
C     $Id$
C     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
C     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
C     i0 ( p5 p6 p7 p8 h1 h2 h3 h4 )_tt + = 1/2 * P( 36 ) * t ( p5 p6 h1 h2 )_t * i1 ( p7 p8 h3 h4 )_t
C         i1 ( p5 p6 h1 h2 )_t + = 1 * t ( p5 p6 h1 h2 )_t
C         i1 ( p5 p6 h1 h2 )_tt + = 1 * P( 4 ) * t ( p5 h1 )_t * t ( p6 h2 )_t
C     i0 ( p5 p6 p7 p8 h1 h2 h3 h4 )_tttt + = -3/2 * P( 16 ) * t ( p5 h1 )_t * i1 ( p6 p7 p8 h2 h3 h4 )_ttt
C         i1 ( p5 p6 p7 h1 h2 h3 )_ttt + = 1/9 * P( 9 ) * t ( p5 h1 )_t * i2 ( p6 p7 h2 h3 )_tt
C             i2 ( p5 p6 h1 h2 )_tt + = -1/4 * P( 4 ) * t ( p5 h1 )_t * t ( p6 h2 )_t
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
      INTEGER t_p5b
      INTEGER t_p6b
      INTEGER t_p7b
      INTEGER t_p8b
      INTEGER t_h1b
      INTEGER t_h2b
      INTEGER t_h3b
      INTEGER t_h4b
      INTEGER toggle
      INTEGER d_t2
      INTEGER k_t2_offset
      INTEGER d_i1_1
      INTEGER k_i1_offset_1
      INTEGER l_i1_offset_1
      INTEGER d_t1
      INTEGER k_t1_offset
      INTEGER d_i1_2
      INTEGER k_i1_offset_2
      INTEGER l_i1_offset_2
      INTEGER size_i1_1
      INTEGER d_i2
      INTEGER k_i2_offset
      INTEGER size_i1_2
      INTEGER l_i2_offset
      INTEGER size_i2
      DOUBLE PRECISION a_i0(*)
      CHARACTER*255 filename
      IF (toggle .eq. 3) THEN
      CALL DELETEFILE(d_i1_2)
      IF (.not.MA_POP_STACK(l_i1_offset_2)) CALL ERRQUIT('lr_12prod_q',-
     &1,MA_ERR)
      END IF
      IF (toggle .eq. 3) THEN
      CALL DELETEFILE(d_i1_1)
      IF (.not.MA_POP_STACK(l_i1_offset_1)) CALL ERRQUIT('lr_12prod_q',-
     &1,MA_ERR)
      END IF
      IF (toggle .eq. 1) CALL OFFSET_lr_12prod_q_1_1(l_i1_offset_1,k_i1_
     &offset_1,size_i1_1)
      IF (toggle .eq. 1) CALL TCE_FILENAME('lr_12prod_q_1_1_i1',filename
     &)
      IF (toggle .eq. 1) CALL CREATEFILE(filename,d_i1_1,size_i1_1)
      IF (toggle .eq. 1) CALL lr_12prod_q_1_1(d_t2,k_t2_offset,d_i1_1,k_
     &i1_offset_1)
      IF (toggle .eq. 1) CALL lr_12prod_q_1_2(d_t1,k_t1_offset,d_t1,k_t1
     &_offset,d_i1_1,k_i1_offset_1)
      IF (toggle .eq. 1) CALL RECONCILEFILE(d_i1_1,size_i1_1)
      IF (toggle .eq. 2) CALL lr_12prod_q_1(d_t2,k_t2_offset,d_i1_1,k_i1
     &_offset_1,a_i0,t_p5b,t_p6b,t_p7b,t_p8b,t_h1b,t_h2b,t_h3b,t_h4b)
      IF (toggle .eq. 1) CALL OFFSET_lr_12prod_q_2_1(l_i1_offset_2,k_i1_
     &offset_2,size_i1_2)
      IF (toggle .eq. 1) CALL TCE_FILENAME('lr_12prod_q_2_1_i1',filename
     &)
      IF (toggle .eq. 1) CALL CREATEFILE(filename,d_i1_2,size_i1_2)
      IF (toggle .eq. 1) CALL OFFSET_lr_12prod_q_2_1_1(l_i2_offset,k_i2_
     &offset,size_i2)
      IF (toggle .eq. 1) CALL TCE_FILENAME('lr_12prod_q_2_1_1_i2',filena
     &me)
      IF (toggle .eq. 1) CALL CREATEFILE(filename,d_i2,size_i2)
      IF (toggle .eq. 1) CALL lr_12prod_q_2_1_1(d_t1,k_t1_offset,d_t1,k_
     &t1_offset,d_i2,k_i2_offset)
      IF (toggle .eq. 1) CALL RECONCILEFILE(d_i2,size_i2)
      IF (toggle .eq. 1) CALL lr_12prod_q_2_1(d_t1,k_t1_offset,d_i2,k_i2
     &_offset,d_i1_2,k_i1_offset_2)
      IF (toggle .eq. 1) THEN
      CALL DELETEFILE(d_i2)
      IF (.not.MA_POP_STACK(l_i2_offset)) CALL ERRQUIT('lr_12prod_q',-1,
     &MA_ERR)
      END IF
      IF (toggle .eq. 1) CALL RECONCILEFILE(d_i1_2,size_i1_2)
      IF (toggle .eq. 2) CALL lr_12prod_q_2(d_t1,k_t1_offset,d_i1_2,k_i1
     &_offset_2,a_i0,t_p5b,t_p6b,t_p7b,t_p8b,t_h1b,t_h2b,t_h3b,t_h4b)
      RETURN
      END
      SUBROUTINE lr_12prod_q_1(d_a,k_a_offset,d_b,k_b_offset,a_c,t_p5b,t
     &_p6b,t_p7b,t_p8b,t_h1b,t_h2b,t_h3b,t_h4b)
C     $Id$
C     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
C     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
C     i0 ( p5 p6 p7 p8 h1 h2 h3 h4 )_tt + = 1/2 * P( 36 ) * t ( p5 p6 h1 h2 )_t * i1 ( p7 p8 h3 h4 )_t
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
      INTEGER d_a
      INTEGER k_a_offset
      INTEGER d_b
      INTEGER k_b_offset
      INTEGER t_p5b
      INTEGER t_p6b
      INTEGER t_p7b
      INTEGER t_p8b
      INTEGER t_h1b
      INTEGER t_h2b
      INTEGER t_h3b
      INTEGER t_h4b
      INTEGER p5b
      INTEGER p6b
      INTEGER p7b
      INTEGER p8b
      INTEGER h1b
      INTEGER h2b
      INTEGER h3b
      INTEGER h4b
      INTEGER dimc
      INTEGER l_c_sort
      INTEGER k_c_sort
      INTEGER p5b_1
      INTEGER p6b_1
      INTEGER h1b_1
      INTEGER h2b_1
      INTEGER p7b_2
      INTEGER p8b_2
      INTEGER h3b_2
      INTEGER h4b_2
      INTEGER dim_common
      INTEGER dima_sort
      INTEGER dima
      INTEGER dimb_sort
      INTEGER dimb
      INTEGER l_a_sort
      INTEGER k_a_sort
      INTEGER l_a
      INTEGER k_a
      INTEGER l_b_sort
      INTEGER k_b_sort
      INTEGER l_b
      INTEGER k_b
      DOUBLE PRECISION a_c(*)
      LOGICAL skipped
      DO p5b = noab+1,noab+nvab
      DO p6b = p5b,noab+nvab
      DO p7b = noab+1,noab+nvab
      DO p8b = p7b,noab+nvab
      DO h1b = 1,noab
      DO h2b = h1b,noab
      DO h3b = 1,noab
      DO h4b = h3b,noab
      skipped = .true.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &4b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &4b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &4b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &4b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &4b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p8b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p8b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p8b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p8b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p8b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p8b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &4b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h2b)) skipped = .fals
     &e.
      IF (.not.skipped) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1
     &)+int_mb(k_spin+p7b-1)+int_mb(k_spin+p8b-1)+int_mb(k_spin+h1b-1)+i
     &nt_mb(k_spin+h2b-1)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.1
     &6)) THEN
      IF (int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1)+int_mb(k_spin+p7b-1)
     &+int_mb(k_spin+p8b-1) .eq. int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-
     &1)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),ieor(int_mb(k_sym+p6b-1),ieor(int_mb(
     &k_sym+p7b-1),ieor(int_mb(k_sym+p8b-1),ieor(int_mb(k_sym+h1b-1),ieo
     &r(int_mb(k_sym+h2b-1),ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1)
     &))))))) .eq. ieor(irrep_t,irrep_t)) THEN
      dimc = int_mb(k_range+p5b-1) * int_mb(k_range+p6b-1) * int_mb(k_ra
     &nge+p7b-1) * int_mb(k_range+p8b-1) * int_mb(k_range+h1b-1) * int_m
     &b(k_range+h2b-1) * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c_sort,k_c_sort)) CALL
     & ERRQUIT('lr_12prod_q_1',0,MA_ERR)
      CALL DFILL(dimc,0.0d0,dbl_mb(k_c_sort),1)
      IF (int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1) .eq. int_mb(k_spin+h
     &1b-1)+int_mb(k_spin+h2b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),ieor(int_mb(k_sym+p6b-1),ieor(int_mb(
     &k_sym+h1b-1),int_mb(k_sym+h2b-1)))) .eq. irrep_t) THEN
      CALL TCE_RESTRICTED_4(p5b,p6b,h1b,h2b,p5b_1,p6b_1,h1b_1,h2b_1)
      CALL TCE_RESTRICTED_4(p7b,p8b,h3b,h4b,p7b_2,p8b_2,h3b_2,h4b_2)
      dim_common = 1
      dima_sort = int_mb(k_range+p5b-1) * int_mb(k_range+p6b-1) * int_mb
     &(k_range+h1b-1) * int_mb(k_range+h2b-1)
      dima = dim_common * dima_sort
      dimb_sort = int_mb(k_range+p7b-1) * int_mb(k_range+p8b-1) * int_mb
     &(k_range+h3b-1) * int_mb(k_range+h4b-1)
      dimb = dim_common * dimb_sort
      IF ((dima .gt. 0) .and. (dimb .gt. 0)) THEN
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a_sort,k_a_sort)) CALL
     & ERRQUIT('lr_12prod_q_1',1,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a,k_a)) CALL ERRQUIT('
     &lr_12prod_q_1',2,MA_ERR)
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h2b_1
     & - 1 + noab * (h1b_1 - 1 + noab * (p6b_1 - noab - 1 + nvab * (p5b_
     &1 - noab - 1)))))
      CALL TCE_SORT_4(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p5b-1)
     &,int_mb(k_range+p6b-1),int_mb(k_range+h1b-1),int_mb(k_range+h2b-1)
     &,4,3,2,1,1.0d0)
      IF (.not.MA_POP_STACK(l_a)) CALL ERRQUIT('lr_12prod_q_1',3,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b_sort,k_b_sort)) CALL
     & ERRQUIT('lr_12prod_q_1',4,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b,k_b)) CALL ERRQUIT('
     &lr_12prod_q_1',5,MA_ERR)
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),(h4b_2
     & - 1 + noab * (h3b_2 - 1 + noab * (p8b_2 - noab - 1 + nvab * (p7b_
     &2 - noab - 1)))))
      CALL TCE_SORT_4(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+p7b-1)
     &,int_mb(k_range+p8b-1),int_mb(k_range+h3b-1),int_mb(k_range+h4b-1)
     &,4,3,2,1,1.0d0)
      IF (.not.MA_POP_STACK(l_b)) CALL ERRQUIT('lr_12prod_q_1',6,MA_ERR)
      CALL DGEMM('T','N',dima_sort,dimb_sort,dim_common,1.0d0,dbl_mb(k_a
     &_sort),dim_common,dbl_mb(k_b_sort),dim_common,1.0d0,dbl_mb(k_c_sor
     &t),dima_sort)
      IF (.not.MA_POP_STACK(l_b_sort)) CALL ERRQUIT('lr_12prod_q_1',7,MA
     &_ERR)
      IF (.not.MA_POP_STACK(l_a_sort)) CALL ERRQUIT('lr_12prod_q_1',8,MA
     &_ERR)
      END IF
      END IF
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,7,4,3,6,5,2,1,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,7,4,3,2,6,5,1,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,7,4,3,2,6,1,5,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,7,4,3,6,2,5,1,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,7,4,3,6,2,1,5,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &4b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,7,4,3,2,1,6,5,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,8,7,3,6,5,2,1,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,8,7,3,2,6,5,1,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,8,7,3,2,6,1,5,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,8,7,3,6,2,5,1,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,8,7,3,6,2,1,5,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &4b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,8,7,3,2,1,6,5,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,8,3,7,6,5,2,1,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,8,3,7,2,6,5,1,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,8,3,7,2,6,1,5,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,8,3,7,6,2,5,1,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,8,3,7,6,2,1,5,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &4b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,8,3,7,2,1,6,5,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,4,7,3,6,5,2,1,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,4,7,3,2,6,5,1,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,4,7,3,2,6,1,5,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,4,7,3,6,2,5,1,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,4,7,3,6,2,1,5,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p6b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &4b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,4,7,3,2,1,6,5,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,4,3,7,6,5,2,1,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,4,3,7,2,6,5,1,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,4,3,7,2,6,1,5,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,4,3,7,6,2,5,1,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,4,3,7,6,2,1,5,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &4b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),8,4,3,7,2,1,6,5,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p8b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,3,8,7,6,5,2,1,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p8b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,3,8,7,2,6,5,1,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p8b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,3,8,7,2,6,1,5,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p8b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h2b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,3,8,7,6,2,5,1,-1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p8b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,3,8,7,6,2,1,5,1.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p7b) .and. (t_p6b .eq. p8b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p6b) .and. (t_h1b .eq. h3b) .and. (t_h2b .eq. h
     &4b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h2b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+p8b-1),int_mb(k_range+p7b-1),int_
     &mb(k_range+h2b-1),int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_
     &mb(k_range+p5b-1),4,3,8,7,2,1,6,5,1.0d0/2.0d0)
      END IF
      IF (.not.MA_POP_STACK(l_c_sort)) CALL ERRQUIT('lr_12prod_q_1',9,MA
     &_ERR)
      END IF
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      RETURN
      END
      SUBROUTINE lr_12prod_q_1_1(d_a,k_a_offset,d_c,k_c_offset)
C     $Id$
C     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
C     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
C     i1 ( p5 p6 h1 h2 )_t + = 1 * t ( p5 p6 h1 h2 )_t
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
      INTEGER d_a
      INTEGER k_a_offset
      INTEGER d_c
      INTEGER k_c_offset
      INTEGER nxtask
      INTEGER next
      INTEGER nprocs
      INTEGER count
      INTEGER p5b
      INTEGER p6b
      INTEGER h1b
      INTEGER h2b
      INTEGER dimc
      INTEGER p5b_1
      INTEGER p6b_1
      INTEGER h1b_1
      INTEGER h2b_1
      INTEGER dim_common
      INTEGER dima_sort
      INTEGER dima
      INTEGER l_a_sort
      INTEGER k_a_sort
      INTEGER l_a
      INTEGER k_a
      INTEGER l_c
      INTEGER k_c
      EXTERNAL nxtask
      nprocs = GA_NNODES()
      count = 0
      next = nxtask(nprocs,1)
      DO p5b = noab+1,noab+nvab
      DO p6b = p5b,noab+nvab
      DO h1b = 1,noab
      DO h2b = h1b,noab
      IF (next.eq.count) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1
     &)+int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1).ne.8)) THEN
      IF (int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1) .eq. int_mb(k_spin+h
     &1b-1)+int_mb(k_spin+h2b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),ieor(int_mb(k_sym+p6b-1),ieor(int_mb(
     &k_sym+h1b-1),int_mb(k_sym+h2b-1)))) .eq. irrep_t) THEN
      dimc = int_mb(k_range+p5b-1) * int_mb(k_range+p6b-1) * int_mb(k_ra
     &nge+h1b-1) * int_mb(k_range+h2b-1)
      CALL TCE_RESTRICTED_4(p5b,p6b,h1b,h2b,p5b_1,p6b_1,h1b_1,h2b_1)
      dim_common = 1
      dima_sort = int_mb(k_range+p5b-1) * int_mb(k_range+p6b-1) * int_mb
     &(k_range+h1b-1) * int_mb(k_range+h2b-1)
      dima = dim_common * dima_sort
      IF (dima .gt. 0) THEN
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a_sort,k_a_sort)) CALL
     & ERRQUIT('lr_12prod_q_1_1',0,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a,k_a)) CALL ERRQUIT('
     &lr_12prod_q_1_1',1,MA_ERR)
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h2b_1
     & - 1 + noab * (h1b_1 - 1 + noab * (p6b_1 - noab - 1 + nvab * (p5b_
     &1 - noab - 1)))))
      CALL TCE_SORT_4(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p5b-1)
     &,int_mb(k_range+p6b-1),int_mb(k_range+h1b-1),int_mb(k_range+h2b-1)
     &,4,3,2,1,1.0d0)
      IF (.not.MA_POP_STACK(l_a)) CALL ERRQUIT('lr_12prod_q_1_1',2,MA_ER
     &R)
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c,k_c)) CALL ERRQUIT('
     &lr_12prod_q_1_1',3,MA_ERR)
      CALL TCE_SORT_4(dbl_mb(k_a_sort),dbl_mb(k_c),int_mb(k_range+h2b-1)
     &,int_mb(k_range+h1b-1),int_mb(k_range+p6b-1),int_mb(k_range+p5b-1)
     &,4,3,2,1,1.0d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(h2b -
     & 1 + noab * (h1b - 1 + noab * (p6b - noab - 1 + nvab * (p5b - noab
     & - 1)))))
      IF (.not.MA_POP_STACK(l_c)) CALL ERRQUIT('lr_12prod_q_1_1',4,MA_ER
     &R)
      IF (.not.MA_POP_STACK(l_a_sort)) CALL ERRQUIT('lr_12prod_q_1_1',5,
     &MA_ERR)
      END IF
      END IF
      END IF
      END IF
      next = nxtask(nprocs,1)
      END IF
      count = count + 1
      END DO
      END DO
      END DO
      END DO
      next = nxtask(-nprocs,1)
      call GA_SYNC()
      RETURN
      END
      SUBROUTINE OFFSET_lr_12prod_q_1_1(l_a_offset,k_a_offset,size)
C     $Id$
C     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
C     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
C     i1 ( p5 p6 h1 h2 )_t
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
      INTEGER l_a_offset
      INTEGER k_a_offset
      INTEGER size
      INTEGER length
      INTEGER addr
      INTEGER p5b
      INTEGER p6b
      INTEGER h1b
      INTEGER h2b
      length = 0
      DO p5b = noab+1,noab+nvab
      DO p6b = p5b,noab+nvab
      DO h1b = 1,noab
      DO h2b = h1b,noab
      IF (int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1) .eq. int_mb(k_spin+h
     &1b-1)+int_mb(k_spin+h2b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),ieor(int_mb(k_sym+p6b-1),ieor(int_mb(
     &k_sym+h1b-1),int_mb(k_sym+h2b-1)))) .eq. irrep_t) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1
     &)+int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1).ne.8)) THEN
      length = length + 1
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
      IF (.not.MA_PUSH_GET(mt_int,2*length+1,'noname',l_a_offset,k_a_off
     &set)) CALL ERRQUIT('lr_12prod_q_1_1',0,MA_ERR)
      int_mb(k_a_offset) = length
      addr = 0
      size = 0
      DO p5b = noab+1,noab+nvab
      DO p6b = p5b,noab+nvab
      DO h1b = 1,noab
      DO h2b = h1b,noab
      IF (int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1) .eq. int_mb(k_spin+h
     &1b-1)+int_mb(k_spin+h2b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),ieor(int_mb(k_sym+p6b-1),ieor(int_mb(
     &k_sym+h1b-1),int_mb(k_sym+h2b-1)))) .eq. irrep_t) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1
     &)+int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1).ne.8)) THEN
      addr = addr + 1
      int_mb(k_a_offset+addr) = h2b - 1 + noab * (h1b - 1 + noab * (p6b 
     &- noab - 1 + nvab * (p5b - noab - 1)))
      int_mb(k_a_offset+length+addr) = size
      size = size + int_mb(k_range+p5b-1) * int_mb(k_range+p6b-1) * int_
     &mb(k_range+h1b-1) * int_mb(k_range+h2b-1)
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
      RETURN
      END
      SUBROUTINE lr_12prod_q_1_2(d_a,k_a_offset,d_b,k_b_offset,d_c,k_c_o
     &ffset)
C     $Id$
C     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
C     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
C     i1 ( p5 p6 h1 h2 )_tt + = 1 * P( 4 ) * t ( p5 h1 )_t * t ( p6 h2 )_t
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
      INTEGER d_a
      INTEGER k_a_offset
      INTEGER d_b
      INTEGER k_b_offset
      INTEGER d_c
      INTEGER k_c_offset
      INTEGER nxtask
      INTEGER next
      INTEGER nprocs
      INTEGER count
      INTEGER p5b
      INTEGER p6b
      INTEGER h1b
      INTEGER h2b
      INTEGER dimc
      INTEGER l_c_sort
      INTEGER k_c_sort
      INTEGER p5b_1
      INTEGER h1b_1
      INTEGER p6b_2
      INTEGER h2b_2
      INTEGER dim_common
      INTEGER dima_sort
      INTEGER dima
      INTEGER dimb_sort
      INTEGER dimb
      INTEGER l_a_sort
      INTEGER k_a_sort
      INTEGER l_a
      INTEGER k_a
      INTEGER l_b_sort
      INTEGER k_b_sort
      INTEGER l_b
      INTEGER k_b
      INTEGER l_c
      INTEGER k_c
      EXTERNAL nxtask
      nprocs = GA_NNODES()
      count = 0
      next = nxtask(nprocs,1)
      DO p5b = noab+1,noab+nvab
      DO p6b = noab+1,noab+nvab
      DO h1b = 1,noab
      DO h2b = 1,noab
      IF (next.eq.count) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1
     &)+int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1).ne.8)) THEN
      IF (int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1) .eq. int_mb(k_spin+h
     &1b-1)+int_mb(k_spin+h2b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),ieor(int_mb(k_sym+p6b-1),ieor(int_mb(
     &k_sym+h1b-1),int_mb(k_sym+h2b-1)))) .eq. ieor(irrep_t,irrep_t)) TH
     &EN
      dimc = int_mb(k_range+p5b-1) * int_mb(k_range+p6b-1) * int_mb(k_ra
     &nge+h1b-1) * int_mb(k_range+h2b-1)
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c_sort,k_c_sort)) CALL
     & ERRQUIT('lr_12prod_q_1_2',0,MA_ERR)
      CALL DFILL(dimc,0.0d0,dbl_mb(k_c_sort),1)
      IF (int_mb(k_spin+p5b-1) .eq. int_mb(k_spin+h1b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),int_mb(k_sym+h1b-1)) .eq. irrep_t) TH
     &EN
      CALL TCE_RESTRICTED_2(p5b,h1b,p5b_1,h1b_1)
      CALL TCE_RESTRICTED_2(p6b,h2b,p6b_2,h2b_2)
      dim_common = 1
      dima_sort = int_mb(k_range+p5b-1) * int_mb(k_range+h1b-1)
      dima = dim_common * dima_sort
      dimb_sort = int_mb(k_range+p6b-1) * int_mb(k_range+h2b-1)
      dimb = dim_common * dimb_sort
      IF ((dima .gt. 0) .and. (dimb .gt. 0)) THEN
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a_sort,k_a_sort)) CALL
     & ERRQUIT('lr_12prod_q_1_2',1,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a,k_a)) CALL ERRQUIT('
     &lr_12prod_q_1_2',2,MA_ERR)
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h1b_1
     & - 1 + noab * (p5b_1 - noab - 1)))
      CALL TCE_SORT_2(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p5b-1)
     &,int_mb(k_range+h1b-1),2,1,1.0d0)
      IF (.not.MA_POP_STACK(l_a)) CALL ERRQUIT('lr_12prod_q_1_2',3,MA_ER
     &R)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b_sort,k_b_sort)) CALL
     & ERRQUIT('lr_12prod_q_1_2',4,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b,k_b)) CALL ERRQUIT('
     &lr_12prod_q_1_2',5,MA_ERR)
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),(h2b_2
     & - 1 + noab * (p6b_2 - noab - 1)))
      CALL TCE_SORT_2(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+p6b-1)
     &,int_mb(k_range+h2b-1),2,1,1.0d0)
      IF (.not.MA_POP_STACK(l_b)) CALL ERRQUIT('lr_12prod_q_1_2',6,MA_ER
     &R)
      CALL DGEMM('T','N',dima_sort,dimb_sort,dim_common,1.0d0,dbl_mb(k_a
     &_sort),dim_common,dbl_mb(k_b_sort),dim_common,1.0d0,dbl_mb(k_c_sor
     &t),dima_sort)
      IF (.not.MA_POP_STACK(l_b_sort)) CALL ERRQUIT('lr_12prod_q_1_2',7,
     &MA_ERR)
      IF (.not.MA_POP_STACK(l_a_sort)) CALL ERRQUIT('lr_12prod_q_1_2',8,
     &MA_ERR)
      END IF
      END IF
      END IF
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c,k_c)) CALL ERRQUIT('
     &lr_12prod_q_1_2',9,MA_ERR)
      IF ((p5b .le. p6b) .and. (h1b .le. h2b)) THEN
      CALL TCE_SORT_4(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+h2b-1)
     &,int_mb(k_range+p6b-1),int_mb(k_range+h1b-1),int_mb(k_range+p5b-1)
     &,4,2,3,1,1.0d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(h2b -
     & 1 + noab * (h1b - 1 + noab * (p6b - noab - 1 + nvab * (p5b - noab
     & - 1)))))
      END IF
      IF ((p5b .le. p6b) .and. (h2b .le. h1b)) THEN
      CALL TCE_SORT_4(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+h2b-1)
     &,int_mb(k_range+p6b-1),int_mb(k_range+h1b-1),int_mb(k_range+p5b-1)
     &,4,2,1,3,-1.0d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(h1b -
     & 1 + noab * (h2b - 1 + noab * (p6b - noab - 1 + nvab * (p5b - noab
     & - 1)))))
      END IF
      IF ((p6b .le. p5b) .and. (h1b .le. h2b)) THEN
      CALL TCE_SORT_4(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+h2b-1)
     &,int_mb(k_range+p6b-1),int_mb(k_range+h1b-1),int_mb(k_range+p5b-1)
     &,2,4,3,1,-1.0d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(h2b -
     & 1 + noab * (h1b - 1 + noab * (p5b - noab - 1 + nvab * (p6b - noab
     & - 1)))))
      END IF
      IF ((p6b .le. p5b) .and. (h2b .le. h1b)) THEN
      CALL TCE_SORT_4(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+h2b-1)
     &,int_mb(k_range+p6b-1),int_mb(k_range+h1b-1),int_mb(k_range+p5b-1)
     &,2,4,1,3,1.0d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(h1b -
     & 1 + noab * (h2b - 1 + noab * (p5b - noab - 1 + nvab * (p6b - noab
     & - 1)))))
      END IF
      IF (.not.MA_POP_STACK(l_c)) CALL ERRQUIT('lr_12prod_q_1_2',10,MA_E
     &RR)
      IF (.not.MA_POP_STACK(l_c_sort)) CALL ERRQUIT('lr_12prod_q_1_2',11
     &,MA_ERR)
      END IF
      END IF
      END IF
      next = nxtask(nprocs,1)
      END IF
      count = count + 1
      END DO
      END DO
      END DO
      END DO
      next = nxtask(-nprocs,1)
      call GA_SYNC()
      RETURN
      END
      SUBROUTINE lr_12prod_q_2(d_a,k_a_offset,d_b,k_b_offset,a_c,t_p5b,t
     &_p6b,t_p7b,t_p8b,t_h1b,t_h2b,t_h3b,t_h4b)
C     $Id$
C     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
C     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
C     i0 ( p5 p6 p7 p8 h1 h2 h3 h4 )_tttt + = -3/2 * P( 16 ) * t ( p5 h1 )_t * i1 ( p6 p7 p8 h2 h3 h4 )_ttt
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
      INTEGER d_a
      INTEGER k_a_offset
      INTEGER d_b
      INTEGER k_b_offset
      INTEGER t_p5b
      INTEGER t_p6b
      INTEGER t_p7b
      INTEGER t_p8b
      INTEGER t_h1b
      INTEGER t_h2b
      INTEGER t_h3b
      INTEGER t_h4b
      INTEGER p5b
      INTEGER p6b
      INTEGER p7b
      INTEGER p8b
      INTEGER h1b
      INTEGER h2b
      INTEGER h3b
      INTEGER h4b
      INTEGER dimc
      INTEGER l_c_sort
      INTEGER k_c_sort
      INTEGER p5b_1
      INTEGER h1b_1
      INTEGER p6b_2
      INTEGER p7b_2
      INTEGER p8b_2
      INTEGER h2b_2
      INTEGER h3b_2
      INTEGER h4b_2
      INTEGER dim_common
      INTEGER dima_sort
      INTEGER dima
      INTEGER dimb_sort
      INTEGER dimb
      INTEGER l_a_sort
      INTEGER k_a_sort
      INTEGER l_a
      INTEGER k_a
      INTEGER l_b_sort
      INTEGER k_b_sort
      INTEGER l_b
      INTEGER k_b
      DOUBLE PRECISION a_c(*)
      LOGICAL skipped
      DO p5b = noab+1,noab+nvab
      DO p6b = noab+1,noab+nvab
      DO p7b = p6b,noab+nvab
      DO p8b = p7b,noab+nvab
      DO h1b = 1,noab
      DO h2b = 1,noab
      DO h3b = h2b,noab
      DO h4b = h3b,noab
      skipped = .true.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h1b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p6b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p6b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p6b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p6b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h1b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p6b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p6b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p6b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p6b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h1b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p6b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p5b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p6b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p5b) .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p6b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p5b) .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h4b)) skipped = .fals
     &e.
      IF ((t_p5b .eq. p6b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p5b) .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h1b)) skipped = .fals
     &e.
      IF (.not.skipped) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1
     &)+int_mb(k_spin+p7b-1)+int_mb(k_spin+p8b-1)+int_mb(k_spin+h1b-1)+i
     &nt_mb(k_spin+h2b-1)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.1
     &6)) THEN
      IF (int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1)+int_mb(k_spin+p7b-1)
     &+int_mb(k_spin+p8b-1) .eq. int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-
     &1)+int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),ieor(int_mb(k_sym+p6b-1),ieor(int_mb(
     &k_sym+p7b-1),ieor(int_mb(k_sym+p8b-1),ieor(int_mb(k_sym+h1b-1),ieo
     &r(int_mb(k_sym+h2b-1),ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1)
     &))))))) .eq. ieor(irrep_t,ieor(irrep_t,ieor(irrep_t,irrep_t)))) TH
     &EN
      dimc = int_mb(k_range+p5b-1) * int_mb(k_range+p6b-1) * int_mb(k_ra
     &nge+p7b-1) * int_mb(k_range+p8b-1) * int_mb(k_range+h1b-1) * int_m
     &b(k_range+h2b-1) * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c_sort,k_c_sort)) CALL
     & ERRQUIT('lr_12prod_q_2',0,MA_ERR)
      CALL DFILL(dimc,0.0d0,dbl_mb(k_c_sort),1)
      IF (int_mb(k_spin+p5b-1) .eq. int_mb(k_spin+h1b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),int_mb(k_sym+h1b-1)) .eq. irrep_t) TH
     &EN
      CALL TCE_RESTRICTED_2(p5b,h1b,p5b_1,h1b_1)
      CALL TCE_RESTRICTED_6(p6b,p7b,p8b,h2b,h3b,h4b,p6b_2,p7b_2,p8b_2,h2
     &b_2,h3b_2,h4b_2)
      dim_common = 1
      dima_sort = int_mb(k_range+p5b-1) * int_mb(k_range+h1b-1)
      dima = dim_common * dima_sort
      dimb_sort = int_mb(k_range+p6b-1) * int_mb(k_range+p7b-1) * int_mb
     &(k_range+p8b-1) * int_mb(k_range+h2b-1) * int_mb(k_range+h3b-1) * 
     &int_mb(k_range+h4b-1)
      dimb = dim_common * dimb_sort
      IF ((dima .gt. 0) .and. (dimb .gt. 0)) THEN
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a_sort,k_a_sort)) CALL
     & ERRQUIT('lr_12prod_q_2',1,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a,k_a)) CALL ERRQUIT('
     &lr_12prod_q_2',2,MA_ERR)
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h1b_1
     & - 1 + noab * (p5b_1 - noab - 1)))
      CALL TCE_SORT_2(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p5b-1)
     &,int_mb(k_range+h1b-1),2,1,1.0d0)
      IF (.not.MA_POP_STACK(l_a)) CALL ERRQUIT('lr_12prod_q_2',3,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b_sort,k_b_sort)) CALL
     & ERRQUIT('lr_12prod_q_2',4,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b,k_b)) CALL ERRQUIT('
     &lr_12prod_q_2',5,MA_ERR)
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),(h4b_2
     & - 1 + noab * (h3b_2 - 1 + noab * (h2b_2 - 1 + noab * (p8b_2 - noa
     &b - 1 + nvab * (p7b_2 - noab - 1 + nvab * (p6b_2 - noab - 1)))))))
      CALL TCE_SORT_6(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+p6b-1)
     &,int_mb(k_range+p7b-1),int_mb(k_range+p8b-1),int_mb(k_range+h2b-1)
     &,int_mb(k_range+h3b-1),int_mb(k_range+h4b-1),6,5,4,3,2,1,1.0d0)
      IF (.not.MA_POP_STACK(l_b)) CALL ERRQUIT('lr_12prod_q_2',6,MA_ERR)
      CALL DGEMM('T','N',dima_sort,dimb_sort,dim_common,1.0d0,dbl_mb(k_a
     &_sort),dim_common,dbl_mb(k_b_sort),dim_common,1.0d0,dbl_mb(k_c_sor
     &t),dima_sort)
      IF (.not.MA_POP_STACK(l_b_sort)) CALL ERRQUIT('lr_12prod_q_2',7,MA
     &_ERR)
      IF (.not.MA_POP_STACK(l_a_sort)) CALL ERRQUIT('lr_12prod_q_2',8,MA
     &_ERR)
      END IF
      END IF
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+h2b-1),int_mb(k_range+p8b-1),int_
     &mb(k_range+p7b-1),int_mb(k_range+p6b-1),int_mb(k_range+h1b-1),int_
     &mb(k_range+p5b-1),8,6,5,4,7,3,2,1,-3.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+h2b-1),int_mb(k_range+p8b-1),int_
     &mb(k_range+p7b-1),int_mb(k_range+p6b-1),int_mb(k_range+h1b-1),int_
     &mb(k_range+p5b-1),8,6,5,4,3,7,2,1,3.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+h2b-1),int_mb(k_range+p8b-1),int_
     &mb(k_range+p7b-1),int_mb(k_range+p6b-1),int_mb(k_range+h1b-1),int_
     &mb(k_range+p5b-1),8,6,5,4,3,2,7,1,-3.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p5b) .and. (t_p6b .eq. p6b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h1b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+h2b-1),int_mb(k_range+p8b-1),int_
     &mb(k_range+p7b-1),int_mb(k_range+p6b-1),int_mb(k_range+h1b-1),int_
     &mb(k_range+p5b-1),8,6,5,4,3,2,1,7,3.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p6b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+h2b-1),int_mb(k_range+p8b-1),int_
     &mb(k_range+p7b-1),int_mb(k_range+p6b-1),int_mb(k_range+h1b-1),int_
     &mb(k_range+p5b-1),6,8,5,4,7,3,2,1,3.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p6b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+h2b-1),int_mb(k_range+p8b-1),int_
     &mb(k_range+p7b-1),int_mb(k_range+p6b-1),int_mb(k_range+h1b-1),int_
     &mb(k_range+p5b-1),6,8,5,4,3,7,2,1,-3.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p6b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+h2b-1),int_mb(k_range+p8b-1),int_
     &mb(k_range+p7b-1),int_mb(k_range+p6b-1),int_mb(k_range+h1b-1),int_
     &mb(k_range+p5b-1),6,8,5,4,3,2,7,1,3.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p6b) .and. (t_p6b .eq. p5b) .and. (t_p7b .eq. p7b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h1b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+h2b-1),int_mb(k_range+p8b-1),int_
     &mb(k_range+p7b-1),int_mb(k_range+p6b-1),int_mb(k_range+h1b-1),int_
     &mb(k_range+p5b-1),6,8,5,4,3,2,1,7,-3.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p6b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+h2b-1),int_mb(k_range+p8b-1),int_
     &mb(k_range+p7b-1),int_mb(k_range+p6b-1),int_mb(k_range+h1b-1),int_
     &mb(k_range+p5b-1),6,5,8,4,7,3,2,1,-3.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p6b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+h2b-1),int_mb(k_range+p8b-1),int_
     &mb(k_range+p7b-1),int_mb(k_range+p6b-1),int_mb(k_range+h1b-1),int_
     &mb(k_range+p5b-1),6,5,8,4,3,7,2,1,3.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p6b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+h2b-1),int_mb(k_range+p8b-1),int_
     &mb(k_range+p7b-1),int_mb(k_range+p6b-1),int_mb(k_range+h1b-1),int_
     &mb(k_range+p5b-1),6,5,8,4,3,2,7,1,-3.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p6b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p5b)
     & .and. (t_p8b .eq. p8b) .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h1b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+h2b-1),int_mb(k_range+p8b-1),int_
     &mb(k_range+p7b-1),int_mb(k_range+p6b-1),int_mb(k_range+h1b-1),int_
     &mb(k_range+p5b-1),6,5,8,4,3,2,1,7,3.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p6b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p5b) .and. (t_h1b .eq. h1b) .and. (t_h2b .eq. h
     &2b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+h2b-1),int_mb(k_range+p8b-1),int_
     &mb(k_range+p7b-1),int_mb(k_range+p6b-1),int_mb(k_range+h1b-1),int_
     &mb(k_range+p5b-1),6,5,4,8,7,3,2,1,3.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p6b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p5b) .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h
     &1b) .and. (t_h3b .eq. h3b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+h2b-1),int_mb(k_range+p8b-1),int_
     &mb(k_range+p7b-1),int_mb(k_range+p6b-1),int_mb(k_range+h1b-1),int_
     &mb(k_range+p5b-1),6,5,4,8,3,7,2,1,-3.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p6b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p5b) .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h1b) .and. (t_h4b .eq. h4b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+h2b-1),int_mb(k_range+p8b-1),int_
     &mb(k_range+p7b-1),int_mb(k_range+p6b-1),int_mb(k_range+h1b-1),int_
     &mb(k_range+p5b-1),6,5,4,8,3,2,7,1,3.0d0/2.0d0)
      END IF
      IF ((t_p5b .eq. p6b) .and. (t_p6b .eq. p7b) .and. (t_p7b .eq. p8b)
     & .and. (t_p8b .eq. p5b) .and. (t_h1b .eq. h2b) .and. (t_h2b .eq. h
     &3b) .and. (t_h3b .eq. h4b) .and. (t_h4b .eq. h1b)) THEN
      CALL TCE_SORTACC_8(dbl_mb(k_c_sort),a_c,int_mb(k_range+h4b-1),int_
     &mb(k_range+h3b-1),int_mb(k_range+h2b-1),int_mb(k_range+p8b-1),int_
     &mb(k_range+p7b-1),int_mb(k_range+p6b-1),int_mb(k_range+h1b-1),int_
     &mb(k_range+p5b-1),6,5,4,8,3,2,1,7,-3.0d0/2.0d0)
      END IF
      IF (.not.MA_POP_STACK(l_c_sort)) CALL ERRQUIT('lr_12prod_q_2',9,MA
     &_ERR)
      END IF
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      RETURN
      END
      SUBROUTINE lr_12prod_q_2_1(d_a,k_a_offset,d_b,k_b_offset,d_c,k_c_o
     &ffset)
C     $Id$
C     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
C     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
C     i1 ( p5 p6 p7 h1 h2 h3 )_ttt + = 1/9 * P( 9 ) * t ( p5 h1 )_t * i2 ( p6 p7 h2 h3 )_tt
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
      INTEGER d_a
      INTEGER k_a_offset
      INTEGER d_b
      INTEGER k_b_offset
      INTEGER d_c
      INTEGER k_c_offset
      INTEGER nxtask
      INTEGER next
      INTEGER nprocs
      INTEGER count
      INTEGER p5b
      INTEGER p6b
      INTEGER p7b
      INTEGER h1b
      INTEGER h2b
      INTEGER h3b
      INTEGER dimc
      INTEGER l_c_sort
      INTEGER k_c_sort
      INTEGER p5b_1
      INTEGER h1b_1
      INTEGER p6b_2
      INTEGER p7b_2
      INTEGER h2b_2
      INTEGER h3b_2
      INTEGER dim_common
      INTEGER dima_sort
      INTEGER dima
      INTEGER dimb_sort
      INTEGER dimb
      INTEGER l_a_sort
      INTEGER k_a_sort
      INTEGER l_a
      INTEGER k_a
      INTEGER l_b_sort
      INTEGER k_b_sort
      INTEGER l_b
      INTEGER k_b
      INTEGER l_c
      INTEGER k_c
      EXTERNAL nxtask
      nprocs = GA_NNODES()
      count = 0
      next = nxtask(nprocs,1)
      DO p5b = noab+1,noab+nvab
      DO p6b = noab+1,noab+nvab
      DO p7b = p6b,noab+nvab
      DO h1b = 1,noab
      DO h2b = 1,noab
      DO h3b = h2b,noab
      IF (next.eq.count) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1
     &)+int_mb(k_spin+p7b-1)+int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1)+i
     &nt_mb(k_spin+h3b-1).ne.12)) THEN
      IF (int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1)+int_mb(k_spin+p7b-1)
     & .eq. int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1)+int_mb(k_spin+h3b-
     &1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),ieor(int_mb(k_sym+p6b-1),ieor(int_mb(
     &k_sym+p7b-1),ieor(int_mb(k_sym+h1b-1),ieor(int_mb(k_sym+h2b-1),int
     &_mb(k_sym+h3b-1)))))) .eq. ieor(irrep_t,ieor(irrep_t,irrep_t))) TH
     &EN
      dimc = int_mb(k_range+p5b-1) * int_mb(k_range+p6b-1) * int_mb(k_ra
     &nge+p7b-1) * int_mb(k_range+h1b-1) * int_mb(k_range+h2b-1) * int_m
     &b(k_range+h3b-1)
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c_sort,k_c_sort)) CALL
     & ERRQUIT('lr_12prod_q_2_1',0,MA_ERR)
      CALL DFILL(dimc,0.0d0,dbl_mb(k_c_sort),1)
      IF (int_mb(k_spin+p5b-1) .eq. int_mb(k_spin+h1b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),int_mb(k_sym+h1b-1)) .eq. irrep_t) TH
     &EN
      CALL TCE_RESTRICTED_2(p5b,h1b,p5b_1,h1b_1)
      CALL TCE_RESTRICTED_4(p6b,p7b,h2b,h3b,p6b_2,p7b_2,h2b_2,h3b_2)
      dim_common = 1
      dima_sort = int_mb(k_range+p5b-1) * int_mb(k_range+h1b-1)
      dima = dim_common * dima_sort
      dimb_sort = int_mb(k_range+p6b-1) * int_mb(k_range+p7b-1) * int_mb
     &(k_range+h2b-1) * int_mb(k_range+h3b-1)
      dimb = dim_common * dimb_sort
      IF ((dima .gt. 0) .and. (dimb .gt. 0)) THEN
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a_sort,k_a_sort)) CALL
     & ERRQUIT('lr_12prod_q_2_1',1,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a,k_a)) CALL ERRQUIT('
     &lr_12prod_q_2_1',2,MA_ERR)
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h1b_1
     & - 1 + noab * (p5b_1 - noab - 1)))
      CALL TCE_SORT_2(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p5b-1)
     &,int_mb(k_range+h1b-1),2,1,1.0d0)
      IF (.not.MA_POP_STACK(l_a)) CALL ERRQUIT('lr_12prod_q_2_1',3,MA_ER
     &R)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b_sort,k_b_sort)) CALL
     & ERRQUIT('lr_12prod_q_2_1',4,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b,k_b)) CALL ERRQUIT('
     &lr_12prod_q_2_1',5,MA_ERR)
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),(h3b_2
     & - 1 + noab * (h2b_2 - 1 + noab * (p7b_2 - noab - 1 + nvab * (p6b_
     &2 - noab - 1)))))
      CALL TCE_SORT_4(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+p6b-1)
     &,int_mb(k_range+p7b-1),int_mb(k_range+h2b-1),int_mb(k_range+h3b-1)
     &,4,3,2,1,1.0d0)
      IF (.not.MA_POP_STACK(l_b)) CALL ERRQUIT('lr_12prod_q_2_1',6,MA_ER
     &R)
      CALL DGEMM('T','N',dima_sort,dimb_sort,dim_common,1.0d0,dbl_mb(k_a
     &_sort),dim_common,dbl_mb(k_b_sort),dim_common,1.0d0,dbl_mb(k_c_sor
     &t),dima_sort)
      IF (.not.MA_POP_STACK(l_b_sort)) CALL ERRQUIT('lr_12prod_q_2_1',7,
     &MA_ERR)
      IF (.not.MA_POP_STACK(l_a_sort)) CALL ERRQUIT('lr_12prod_q_2_1',8,
     &MA_ERR)
      END IF
      END IF
      END IF
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c,k_c)) CALL ERRQUIT('
     &lr_12prod_q_2_1',9,MA_ERR)
      IF ((p5b .le. p6b) .and. (h1b .le. h2b)) THEN
      CALL TCE_SORT_6(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+h3b-1)
     &,int_mb(k_range+h2b-1),int_mb(k_range+p7b-1),int_mb(k_range+p6b-1)
     &,int_mb(k_range+h1b-1),int_mb(k_range+p5b-1),6,4,3,5,2,1,1.0d0/9.0
     &d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(h3b -
     & 1 + noab * (h2b - 1 + noab * (h1b - 1 + noab * (p7b - noab - 1 + 
     &nvab * (p6b - noab - 1 + nvab * (p5b - noab - 1)))))))
      END IF
      IF ((p5b .le. p6b) .and. (h2b .le. h1b) .and. (h1b .le. h3b)) THEN
      CALL TCE_SORT_6(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+h3b-1)
     &,int_mb(k_range+h2b-1),int_mb(k_range+p7b-1),int_mb(k_range+p6b-1)
     &,int_mb(k_range+h1b-1),int_mb(k_range+p5b-1),6,4,3,2,5,1,-1.0d0/9.
     &0d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(h3b -
     & 1 + noab * (h1b - 1 + noab * (h2b - 1 + noab * (p7b - noab - 1 + 
     &nvab * (p6b - noab - 1 + nvab * (p5b - noab - 1)))))))
      END IF
      IF ((p5b .le. p6b) .and. (h3b .le. h1b)) THEN
      CALL TCE_SORT_6(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+h3b-1)
     &,int_mb(k_range+h2b-1),int_mb(k_range+p7b-1),int_mb(k_range+p6b-1)
     &,int_mb(k_range+h1b-1),int_mb(k_range+p5b-1),6,4,3,2,1,5,1.0d0/9.0
     &d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(h1b -
     & 1 + noab * (h3b - 1 + noab * (h2b - 1 + noab * (p7b - noab - 1 + 
     &nvab * (p6b - noab - 1 + nvab * (p5b - noab - 1)))))))
      END IF
      IF ((p6b .le. p5b) .and. (p5b .le. p7b) .and. (h1b .le. h2b)) THEN
      CALL TCE_SORT_6(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+h3b-1)
     &,int_mb(k_range+h2b-1),int_mb(k_range+p7b-1),int_mb(k_range+p6b-1)
     &,int_mb(k_range+h1b-1),int_mb(k_range+p5b-1),4,6,3,5,2,1,-1.0d0/9.
     &0d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(h3b -
     & 1 + noab * (h2b - 1 + noab * (h1b - 1 + noab * (p7b - noab - 1 + 
     &nvab * (p5b - noab - 1 + nvab * (p6b - noab - 1)))))))
      END IF
      IF ((p6b .le. p5b) .and. (p5b .le. p7b) .and. (h2b .le. h1b) .and.
     & (h1b .le. h3b)) THEN
      CALL TCE_SORT_6(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+h3b-1)
     &,int_mb(k_range+h2b-1),int_mb(k_range+p7b-1),int_mb(k_range+p6b-1)
     &,int_mb(k_range+h1b-1),int_mb(k_range+p5b-1),4,6,3,2,5,1,1.0d0/9.0
     &d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(h3b -
     & 1 + noab * (h1b - 1 + noab * (h2b - 1 + noab * (p7b - noab - 1 + 
     &nvab * (p5b - noab - 1 + nvab * (p6b - noab - 1)))))))
      END IF
      IF ((p6b .le. p5b) .and. (p5b .le. p7b) .and. (h3b .le. h1b)) THEN
      CALL TCE_SORT_6(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+h3b-1)
     &,int_mb(k_range+h2b-1),int_mb(k_range+p7b-1),int_mb(k_range+p6b-1)
     &,int_mb(k_range+h1b-1),int_mb(k_range+p5b-1),4,6,3,2,1,5,-1.0d0/9.
     &0d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(h1b -
     & 1 + noab * (h3b - 1 + noab * (h2b - 1 + noab * (p7b - noab - 1 + 
     &nvab * (p5b - noab - 1 + nvab * (p6b - noab - 1)))))))
      END IF
      IF ((p7b .le. p5b) .and. (h1b .le. h2b)) THEN
      CALL TCE_SORT_6(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+h3b-1)
     &,int_mb(k_range+h2b-1),int_mb(k_range+p7b-1),int_mb(k_range+p6b-1)
     &,int_mb(k_range+h1b-1),int_mb(k_range+p5b-1),4,3,6,5,2,1,1.0d0/9.0
     &d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(h3b -
     & 1 + noab * (h2b - 1 + noab * (h1b - 1 + noab * (p5b - noab - 1 + 
     &nvab * (p7b - noab - 1 + nvab * (p6b - noab - 1)))))))
      END IF
      IF ((p7b .le. p5b) .and. (h2b .le. h1b) .and. (h1b .le. h3b)) THEN
      CALL TCE_SORT_6(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+h3b-1)
     &,int_mb(k_range+h2b-1),int_mb(k_range+p7b-1),int_mb(k_range+p6b-1)
     &,int_mb(k_range+h1b-1),int_mb(k_range+p5b-1),4,3,6,2,5,1,-1.0d0/9.
     &0d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(h3b -
     & 1 + noab * (h1b - 1 + noab * (h2b - 1 + noab * (p5b - noab - 1 + 
     &nvab * (p7b - noab - 1 + nvab * (p6b - noab - 1)))))))
      END IF
      IF ((p7b .le. p5b) .and. (h3b .le. h1b)) THEN
      CALL TCE_SORT_6(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+h3b-1)
     &,int_mb(k_range+h2b-1),int_mb(k_range+p7b-1),int_mb(k_range+p6b-1)
     &,int_mb(k_range+h1b-1),int_mb(k_range+p5b-1),4,3,6,2,1,5,1.0d0/9.0
     &d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(h1b -
     & 1 + noab * (h3b - 1 + noab * (h2b - 1 + noab * (p5b - noab - 1 + 
     &nvab * (p7b - noab - 1 + nvab * (p6b - noab - 1)))))))
      END IF
      IF (.not.MA_POP_STACK(l_c)) CALL ERRQUIT('lr_12prod_q_2_1',10,MA_E
     &RR)
      IF (.not.MA_POP_STACK(l_c_sort)) CALL ERRQUIT('lr_12prod_q_2_1',11
     &,MA_ERR)
      END IF
      END IF
      END IF
      next = nxtask(nprocs,1)
      END IF
      count = count + 1
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      next = nxtask(-nprocs,1)
      call GA_SYNC()
      RETURN
      END
      SUBROUTINE OFFSET_lr_12prod_q_2_1(l_a_offset,k_a_offset,size)
C     $Id$
C     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
C     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
C     i1 ( p5 p6 p7 h1 h2 h3 )_ttt
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
      INTEGER l_a_offset
      INTEGER k_a_offset
      INTEGER size
      INTEGER length
      INTEGER addr
      INTEGER p5b
      INTEGER p6b
      INTEGER p7b
      INTEGER h1b
      INTEGER h2b
      INTEGER h3b
      length = 0
      DO p5b = noab+1,noab+nvab
      DO p6b = p5b,noab+nvab
      DO p7b = p6b,noab+nvab
      DO h1b = 1,noab
      DO h2b = h1b,noab
      DO h3b = h2b,noab
      IF (int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1)+int_mb(k_spin+p7b-1)
     & .eq. int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1)+int_mb(k_spin+h3b-
     &1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),ieor(int_mb(k_sym+p6b-1),ieor(int_mb(
     &k_sym+p7b-1),ieor(int_mb(k_sym+h1b-1),ieor(int_mb(k_sym+h2b-1),int
     &_mb(k_sym+h3b-1)))))) .eq. ieor(irrep_t,ieor(irrep_t,irrep_t))) TH
     &EN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1
     &)+int_mb(k_spin+p7b-1)+int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1)+i
     &nt_mb(k_spin+h3b-1).ne.12)) THEN
      length = length + 1
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      IF (.not.MA_PUSH_GET(mt_int,2*length+1,'noname',l_a_offset,k_a_off
     &set)) CALL ERRQUIT('lr_12prod_q_2_1',0,MA_ERR)
      int_mb(k_a_offset) = length
      addr = 0
      size = 0
      DO p5b = noab+1,noab+nvab
      DO p6b = p5b,noab+nvab
      DO p7b = p6b,noab+nvab
      DO h1b = 1,noab
      DO h2b = h1b,noab
      DO h3b = h2b,noab
      IF (int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1)+int_mb(k_spin+p7b-1)
     & .eq. int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1)+int_mb(k_spin+h3b-
     &1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),ieor(int_mb(k_sym+p6b-1),ieor(int_mb(
     &k_sym+p7b-1),ieor(int_mb(k_sym+h1b-1),ieor(int_mb(k_sym+h2b-1),int
     &_mb(k_sym+h3b-1)))))) .eq. ieor(irrep_t,ieor(irrep_t,irrep_t))) TH
     &EN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1
     &)+int_mb(k_spin+p7b-1)+int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1)+i
     &nt_mb(k_spin+h3b-1).ne.12)) THEN
      addr = addr + 1
      int_mb(k_a_offset+addr) = h3b - 1 + noab * (h2b - 1 + noab * (h1b 
     &- 1 + noab * (p7b - noab - 1 + nvab * (p6b - noab - 1 + nvab * (p5
     &b - noab - 1)))))
      int_mb(k_a_offset+length+addr) = size
      size = size + int_mb(k_range+p5b-1) * int_mb(k_range+p6b-1) * int_
     &mb(k_range+p7b-1) * int_mb(k_range+h1b-1) * int_mb(k_range+h2b-1) 
     &* int_mb(k_range+h3b-1)
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      RETURN
      END
      SUBROUTINE lr_12prod_q_2_1_1(d_a,k_a_offset,d_b,k_b_offset,d_c,k_c
     &_offset)
C     $Id$
C     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
C     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
C     i2 ( p5 p6 h1 h2 )_tt + = -1/4 * P( 4 ) * t ( p5 h1 )_t * t ( p6 h2 )_t
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
      INTEGER d_a
      INTEGER k_a_offset
      INTEGER d_b
      INTEGER k_b_offset
      INTEGER d_c
      INTEGER k_c_offset
      INTEGER nxtask
      INTEGER next
      INTEGER nprocs
      INTEGER count
      INTEGER p5b
      INTEGER p6b
      INTEGER h1b
      INTEGER h2b
      INTEGER dimc
      INTEGER l_c_sort
      INTEGER k_c_sort
      INTEGER p5b_1
      INTEGER h1b_1
      INTEGER p6b_2
      INTEGER h2b_2
      INTEGER dim_common
      INTEGER dima_sort
      INTEGER dima
      INTEGER dimb_sort
      INTEGER dimb
      INTEGER l_a_sort
      INTEGER k_a_sort
      INTEGER l_a
      INTEGER k_a
      INTEGER l_b_sort
      INTEGER k_b_sort
      INTEGER l_b
      INTEGER k_b
      INTEGER l_c
      INTEGER k_c
      EXTERNAL nxtask
      nprocs = GA_NNODES()
      count = 0
      next = nxtask(nprocs,1)
      DO p5b = noab+1,noab+nvab
      DO p6b = noab+1,noab+nvab
      DO h1b = 1,noab
      DO h2b = 1,noab
      IF (next.eq.count) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1
     &)+int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1).ne.8)) THEN
      IF (int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1) .eq. int_mb(k_spin+h
     &1b-1)+int_mb(k_spin+h2b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),ieor(int_mb(k_sym+p6b-1),ieor(int_mb(
     &k_sym+h1b-1),int_mb(k_sym+h2b-1)))) .eq. ieor(irrep_t,irrep_t)) TH
     &EN
      dimc = int_mb(k_range+p5b-1) * int_mb(k_range+p6b-1) * int_mb(k_ra
     &nge+h1b-1) * int_mb(k_range+h2b-1)
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c_sort,k_c_sort)) CALL
     & ERRQUIT('lr_12prod_q_2_1_1',0,MA_ERR)
      CALL DFILL(dimc,0.0d0,dbl_mb(k_c_sort),1)
      IF (int_mb(k_spin+p5b-1) .eq. int_mb(k_spin+h1b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),int_mb(k_sym+h1b-1)) .eq. irrep_t) TH
     &EN
      CALL TCE_RESTRICTED_2(p5b,h1b,p5b_1,h1b_1)
      CALL TCE_RESTRICTED_2(p6b,h2b,p6b_2,h2b_2)
      dim_common = 1
      dima_sort = int_mb(k_range+p5b-1) * int_mb(k_range+h1b-1)
      dima = dim_common * dima_sort
      dimb_sort = int_mb(k_range+p6b-1) * int_mb(k_range+h2b-1)
      dimb = dim_common * dimb_sort
      IF ((dima .gt. 0) .and. (dimb .gt. 0)) THEN
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a_sort,k_a_sort)) CALL
     & ERRQUIT('lr_12prod_q_2_1_1',1,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a,k_a)) CALL ERRQUIT('
     &lr_12prod_q_2_1_1',2,MA_ERR)
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h1b_1
     & - 1 + noab * (p5b_1 - noab - 1)))
      CALL TCE_SORT_2(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p5b-1)
     &,int_mb(k_range+h1b-1),2,1,1.0d0)
      IF (.not.MA_POP_STACK(l_a)) CALL ERRQUIT('lr_12prod_q_2_1_1',3,MA_
     &ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b_sort,k_b_sort)) CALL
     & ERRQUIT('lr_12prod_q_2_1_1',4,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b,k_b)) CALL ERRQUIT('
     &lr_12prod_q_2_1_1',5,MA_ERR)
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),(h2b_2
     & - 1 + noab * (p6b_2 - noab - 1)))
      CALL TCE_SORT_2(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+p6b-1)
     &,int_mb(k_range+h2b-1),2,1,1.0d0)
      IF (.not.MA_POP_STACK(l_b)) CALL ERRQUIT('lr_12prod_q_2_1_1',6,MA_
     &ERR)
      CALL DGEMM('T','N',dima_sort,dimb_sort,dim_common,1.0d0,dbl_mb(k_a
     &_sort),dim_common,dbl_mb(k_b_sort),dim_common,1.0d0,dbl_mb(k_c_sor
     &t),dima_sort)
      IF (.not.MA_POP_STACK(l_b_sort)) CALL ERRQUIT('lr_12prod_q_2_1_1',
     &7,MA_ERR)
      IF (.not.MA_POP_STACK(l_a_sort)) CALL ERRQUIT('lr_12prod_q_2_1_1',
     &8,MA_ERR)
      END IF
      END IF
      END IF
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c,k_c)) CALL ERRQUIT('
     &lr_12prod_q_2_1_1',9,MA_ERR)
      IF ((p5b .le. p6b) .and. (h1b .le. h2b)) THEN
      CALL TCE_SORT_4(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+h2b-1)
     &,int_mb(k_range+p6b-1),int_mb(k_range+h1b-1),int_mb(k_range+p5b-1)
     &,4,2,3,1,-1.0d0/4.0d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(h2b -
     & 1 + noab * (h1b - 1 + noab * (p6b - noab - 1 + nvab * (p5b - noab
     & - 1)))))
      END IF
      IF ((p5b .le. p6b) .and. (h2b .le. h1b)) THEN
      CALL TCE_SORT_4(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+h2b-1)
     &,int_mb(k_range+p6b-1),int_mb(k_range+h1b-1),int_mb(k_range+p5b-1)
     &,4,2,1,3,1.0d0/4.0d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(h1b -
     & 1 + noab * (h2b - 1 + noab * (p6b - noab - 1 + nvab * (p5b - noab
     & - 1)))))
      END IF
      IF ((p6b .le. p5b) .and. (h1b .le. h2b)) THEN
      CALL TCE_SORT_4(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+h2b-1)
     &,int_mb(k_range+p6b-1),int_mb(k_range+h1b-1),int_mb(k_range+p5b-1)
     &,2,4,3,1,1.0d0/4.0d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(h2b -
     & 1 + noab * (h1b - 1 + noab * (p5b - noab - 1 + nvab * (p6b - noab
     & - 1)))))
      END IF
      IF ((p6b .le. p5b) .and. (h2b .le. h1b)) THEN
      CALL TCE_SORT_4(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+h2b-1)
     &,int_mb(k_range+p6b-1),int_mb(k_range+h1b-1),int_mb(k_range+p5b-1)
     &,2,4,1,3,-1.0d0/4.0d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(h1b -
     & 1 + noab * (h2b - 1 + noab * (p5b - noab - 1 + nvab * (p6b - noab
     & - 1)))))
      END IF
      IF (.not.MA_POP_STACK(l_c)) CALL ERRQUIT('lr_12prod_q_2_1_1',10,MA
     &_ERR)
      IF (.not.MA_POP_STACK(l_c_sort)) CALL ERRQUIT('lr_12prod_q_2_1_1',
     &11,MA_ERR)
      END IF
      END IF
      END IF
      next = nxtask(nprocs,1)
      END IF
      count = count + 1
      END DO
      END DO
      END DO
      END DO
      next = nxtask(-nprocs,1)
      call GA_SYNC()
      RETURN
      END
      SUBROUTINE OFFSET_lr_12prod_q_2_1_1(l_a_offset,k_a_offset,size)
C     $Id$
C     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
C     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
C     i2 ( p5 p6 h1 h2 )_tt
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
      INTEGER l_a_offset
      INTEGER k_a_offset
      INTEGER size
      INTEGER length
      INTEGER addr
      INTEGER p5b
      INTEGER p6b
      INTEGER h1b
      INTEGER h2b
      length = 0
      DO p5b = noab+1,noab+nvab
      DO p6b = p5b,noab+nvab
      DO h1b = 1,noab
      DO h2b = h1b,noab
      IF (int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1) .eq. int_mb(k_spin+h
     &1b-1)+int_mb(k_spin+h2b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),ieor(int_mb(k_sym+p6b-1),ieor(int_mb(
     &k_sym+h1b-1),int_mb(k_sym+h2b-1)))) .eq. ieor(irrep_t,irrep_t)) TH
     &EN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1
     &)+int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1).ne.8)) THEN
      length = length + 1
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
      IF (.not.MA_PUSH_GET(mt_int,2*length+1,'noname',l_a_offset,k_a_off
     &set)) CALL ERRQUIT('lr_12prod_q_2_1_1',0,MA_ERR)
      int_mb(k_a_offset) = length
      addr = 0
      size = 0
      DO p5b = noab+1,noab+nvab
      DO p6b = p5b,noab+nvab
      DO h1b = 1,noab
      DO h2b = h1b,noab
      IF (int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1) .eq. int_mb(k_spin+h
     &1b-1)+int_mb(k_spin+h2b-1)) THEN
      IF (ieor(int_mb(k_sym+p5b-1),ieor(int_mb(k_sym+p6b-1),ieor(int_mb(
     &k_sym+h1b-1),int_mb(k_sym+h2b-1)))) .eq. ieor(irrep_t,irrep_t)) TH
     &EN
      IF ((.not.restricted).or.(int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1
     &)+int_mb(k_spin+h1b-1)+int_mb(k_spin+h2b-1).ne.8)) THEN
      addr = addr + 1
      int_mb(k_a_offset+addr) = h2b - 1 + noab * (h1b - 1 + noab * (p6b 
     &- noab - 1 + nvab * (p5b - noab - 1)))
      int_mb(k_a_offset+length+addr) = size
      size = size + int_mb(k_range+p5b-1) * int_mb(k_range+p6b-1) * int_
     &mb(k_range+h1b-1) * int_mb(k_range+h2b-1)
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
      RETURN
      END
