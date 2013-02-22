      SUBROUTINE lr_ccsd_t(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                     d_f1,k_f1_offset,d_v2,k_v2_offset,
     2                     d_e,k_e_offset,
     3                     energy1,energy2,energy3,
     4                     energy4,energy5,energy6,
     5                     size_t1)
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
      integer d_e
      integer k_e_offset
      integer t_h1b, t_h1
      integer t_h2b, t_h2
      integer t_h3b, t_h3
      integer t_p4b, t_p4
      integer t_p5b, t_p5
      integer t_p6b, t_p6
      integer k_doubles,l_doubles
      integer k_den,l_den
      integer k_right,l_right
      integer size,i
      integer g_energy
      integer l_hole,k_hole
      integer l_2hole,k_2hole
      integer hole_p_1
      integer hole_p_2
c
      integer nxtask
      integer next
      integer nprocs
      integer count
      integer d_i1_1,d_i1_2,d_i1_3
      integer k_i1_offset_1,l_i1_offset_1
      integer k_i1_offset_2,l_i1_offset_2
      integer k_i1_offset_3,l_i1_offset_3
c - T1/X1 LOCALIZATION -------------------
      integer l_t1_local,k_t1_local
      integer size_t1
c ---------------------------------------
      double precision energy1,energy2,energy3,energy4
      double precision energy5,energy6
      double precision factor
      double precision den0,den1,den2,den3,num1,num2
      external nxtask
c
      k_right = 0
      k_den = 0
c
c - T1/X1 LOCALIZATION ----------
c    opening l_t1_local and l_x1_local
        if (.not.MA_PUSH_GET(mt_dbl,size_t1,'t1_local',
     1      l_t1_local,k_t1_local))
     1      call errquit('t1_local',1,MA_ERR)
        call ma_zero(dbl_mb(k_t1_local),size_t1)
c    copy d_t1 ==> l_t1_local
cc        call ga_get(d_t1,1,size_t1,1,1,dbl_mb(k_t1_local),1)
      call get_block(d_t1,dbl_mb(k_t1_local),size_t1,0)
c -------------------------------
c
c     Get singles & doubles part of the denominator
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
c     Caution! k_right & k_den are not even allocated yet
c     but they won't be used.
c
      call cr_ccsd_t_N(dbl_mb(k_right),d_f1,d_i1_1,d_i1_2,
     1  k_t1_local,d_t2,d_v2,k_f1_offset,k_i1_offset_1,k_i1_offset_2,
     2  k_t1_offset,k_t2_offset,k_v2_offset,l_i1_offset_1,
     3  l_i1_offset_2,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,1)
      call cr_ccsd_t_E(dbl_mb(k_den),d_i1_3,
     1  k_t1_local,d_t2,k_i1_offset_3,k_t1_offset,k_t2_offset,
     2  l_i1_offset_3,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,1)
c
c     Get the numerator
c
      if (.not.ga_create(mt_dbl,1,1,'perturbative',1,1,g_energy))
     1  call errquit('ccsd_t: GA problem',0,GA_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = nxtask(nprocs,1)
      num1=0.0d0 !for <S3(M3)+ M3>/(1+L(J3))
      num2=0.0d0 !for <(S3(M3)+S1S2)+M3>/(1+L(J3))
      den0=0.0d0  !for <T3[2]+M3>/(1+L(J3))
      den1=0.0d0 !for <(T3[2]+S1S2)+M3>/(1+L(J3))
      den2=0.0d0 !for <(T3[2])+R03 M3[2]>/(1+L(J3))
      den3=0.0d0 !for <(T3[2]+S1S2)+R03 M3[2]>/(1+L(J3))
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
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) doubles',
     1        l_doubles,k_doubles)) call errquit('ccsd_t',2,MA_ERR)
            if (.not.MA_PUSH_GET(mt_dbl,size,'moment 2,3',
     1        l_right,k_right)) call errquit('ccsd_t',3,MA_ERR)
            if (.not.MA_PUSH_GET(mt_dbl,size,'denominator',
     1        l_den,k_den)) call errquit('ccsd_t',3,MA_ERR)
            do i = 1, size
             dbl_mb(k_doubles+i-1) = 0.0d0
             dbl_mb(k_right+i-1) = 0.0d0
             dbl_mb(k_den+i-1) = 0.0d0
            enddo
            call ccsd_t_doubles(dbl_mb(k_doubles),d_t2,d_v2,k_t2_offset,
     1        k_v2_offset,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
            call cr_ccsd_t_N(dbl_mb(k_right),d_f1,d_i1_1,d_i1_2,
     1        k_t1_local,d_t2,d_v2,
     &        k_f1_offset,k_i1_offset_1,k_i1_offset_2,
     2        k_t1_offset,k_t2_offset,k_v2_offset,l_i1_offset_1,
     3        l_i1_offset_2,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
            call cr_ccsd_t_E(dbl_mb(k_den),d_i1_3,
     1        k_t1_local,d_t2,k_i1_offset_3,k_t1_offset,k_t2_offset,
     2        l_i1_offset_3,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
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
c num1 => <S3(M3)+ M3>/(1+L(J3))
                  num1 = num1 + factor * dbl_mb(k_right+i-1)
     &                           * dbl_mb(k_right+i-1)/
     &        ((-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p4b-1)+t_p4-1)
     &          -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     &          -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     &          +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     &          +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     &          +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1))
     &         *(1.0d0+dbl_mb(k_hole+int_mb(k_offset+t_h1b-1)+t_h1-1)
     &                +dbl_mb(k_hole+int_mb(k_offset+t_h2b-1)+t_h2-1)
     &                +dbl_mb(k_hole+int_mb(k_offset+t_h3b-1)+t_h3-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h1b-1)+t_h1-1)
     &            +int_mb(k_offset+t_h2b-1)+t_h2-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h1b-1)+t_h1-1)
     &            +int_mb(k_offset+t_h3b-1)+t_h3-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h2b-1)+t_h2-1)
     &            +int_mb(k_offset+t_h3b-1)+t_h3-1)))
c num2 => <(S3(M3)+S1S2)+M3>/(1+L(J3))
                  num2 = num2 + factor * dbl_mb(k_right+i-1)
     &                           * dbl_mb(k_right+i-1)/
     &        ((-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p4b-1)+t_p4-1)
     &          -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     &          -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     &          +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     &          +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     &          +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1))
     &         *(1.0d0+dbl_mb(k_hole+int_mb(k_offset+t_h1b-1)+t_h1-1)
     &                +dbl_mb(k_hole+int_mb(k_offset+t_h2b-1)+t_h2-1)
     &                +dbl_mb(k_hole+int_mb(k_offset+t_h3b-1)+t_h3-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h1b-1)+t_h1-1)
     &            +int_mb(k_offset+t_h2b-1)+t_h2-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h1b-1)+t_h1-1)
     &            +int_mb(k_offset+t_h3b-1)+t_h3-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h2b-1)+t_h2-1)
     &            +int_mb(k_offset+t_h3b-1)+t_h3-1)))  
     &          +factor*dbl_mb(k_den+i-1)*dbl_mb(k_right+i-1)/
     &          (1.0d0+dbl_mb(k_hole+int_mb(k_offset+t_h1b-1)+t_h1-1)
     &                +dbl_mb(k_hole+int_mb(k_offset+t_h2b-1)+t_h2-1)
     &                +dbl_mb(k_hole+int_mb(k_offset+t_h3b-1)+t_h3-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h1b-1)+t_h1-1)
     &            +int_mb(k_offset+t_h2b-1)+t_h2-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h1b-1)+t_h1-1)
     &            +int_mb(k_offset+t_h3b-1)+t_h3-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h2b-1)+t_h2-1)
     &            +int_mb(k_offset+t_h3b-1)+t_h3-1))
c den0 => <T3[2]+M3>/(1+L(J3))
                  den0 = den0 + factor * dbl_mb(k_doubles+i-1)
     &                           * dbl_mb(k_right+i-1)/
     &        ((-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p4b-1)+t_p4-1)
     &          -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     &          -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     &          +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     &          +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     &          +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1))
     &         *(1.0d0+dbl_mb(k_hole+int_mb(k_offset+t_h1b-1)+t_h1-1)
     &                +dbl_mb(k_hole+int_mb(k_offset+t_h2b-1)+t_h2-1)
     &                +dbl_mb(k_hole+int_mb(k_offset+t_h3b-1)+t_h3-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h1b-1)+t_h1-1)
     &            +int_mb(k_offset+t_h2b-1)+t_h2-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h1b-1)+t_h1-1)
     &            +int_mb(k_offset+t_h3b-1)+t_h3-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h2b-1)+t_h2-1)
     &            +int_mb(k_offset+t_h3b-1)+t_h3-1)))
c den1 => <(T3[2]+S1S2)+M3>/(1+L(J3))
                  den1 = den1 + factor * dbl_mb(k_doubles+i-1)
     &                           * dbl_mb(k_right+i-1)/
     &        ((-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p4b-1)+t_p4-1)
     &          -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     &          -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     &          +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     &          +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     &          +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1))
     &         *(1.0d0+dbl_mb(k_hole+int_mb(k_offset+t_h1b-1)+t_h1-1)
     &                +dbl_mb(k_hole+int_mb(k_offset+t_h2b-1)+t_h2-1)
     &                +dbl_mb(k_hole+int_mb(k_offset+t_h3b-1)+t_h3-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h1b-1)+t_h1-1)
     &            +int_mb(k_offset+t_h2b-1)+t_h2-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h1b-1)+t_h1-1)
     &            +int_mb(k_offset+t_h3b-1)+t_h3-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h2b-1)+t_h2-1)
     &            +int_mb(k_offset+t_h3b-1)+t_h3-1)))
     &          +factor*dbl_mb(k_den+i-1)*dbl_mb(k_right+i-1)/
     &          (1.0d0+dbl_mb(k_hole+int_mb(k_offset+t_h1b-1)+t_h1-1)
     &                +dbl_mb(k_hole+int_mb(k_offset+t_h2b-1)+t_h2-1)
     &                +dbl_mb(k_hole+int_mb(k_offset+t_h3b-1)+t_h3-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h1b-1)+t_h1-1)
     &            +int_mb(k_offset+t_h2b-1)+t_h2-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h1b-1)+t_h1-1)
     &            +int_mb(k_offset+t_h3b-1)+t_h3-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h2b-1)+t_h2-1)
     &            +int_mb(k_offset+t_h3b-1)+t_h3-1))
c den2 => <(T3[2])+R03 M3[2]>/(1+L(J3))
                  den2 = den2 + factor * dbl_mb(k_doubles+i-1)
     &                           * dbl_mb(k_doubles+i-1)/
     &        ((-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p4b-1)+t_p4-1)
     &          -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     &          -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     &          +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     &          +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     &          +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1))
     &         *(1.0d0+dbl_mb(k_hole+int_mb(k_offset+t_h1b-1)+t_h1-1)
     &                +dbl_mb(k_hole+int_mb(k_offset+t_h2b-1)+t_h2-1)
     &                +dbl_mb(k_hole+int_mb(k_offset+t_h3b-1)+t_h3-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h1b-1)+t_h1-1)
     &            +int_mb(k_offset+t_h2b-1)+t_h2-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h1b-1)+t_h1-1)
     &            +int_mb(k_offset+t_h3b-1)+t_h3-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h2b-1)+t_h2-1)
     &            +int_mb(k_offset+t_h3b-1)+t_h3-1)))
c den3 => <(T3[2]+S1S2)+R03 M3[2]>/(1+L(J3))
                  den3 = den3 + factor * dbl_mb(k_doubles+i-1)
     &                           * dbl_mb(k_doubles+i-1)/
     &        ((-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p4b-1)+t_p4-1)
     &          -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     &          -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     &          +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     &          +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     &          +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1))
     &         *(1.0d0+dbl_mb(k_hole+int_mb(k_offset+t_h1b-1)+t_h1-1)
     &                +dbl_mb(k_hole+int_mb(k_offset+t_h2b-1)+t_h2-1)
     &                +dbl_mb(k_hole+int_mb(k_offset+t_h3b-1)+t_h3-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h1b-1)+t_h1-1)
     &            +int_mb(k_offset+t_h2b-1)+t_h2-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h1b-1)+t_h1-1)
     &            +int_mb(k_offset+t_h3b-1)+t_h3-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h2b-1)+t_h2-1)
     &            +int_mb(k_offset+t_h3b-1)+t_h3-1)))
     &          +factor*dbl_mb(k_den+i-1)*dbl_mb(k_doubles+i-1)/
     &          (1.0d0+dbl_mb(k_hole+int_mb(k_offset+t_h1b-1)+t_h1-1)
     &                +dbl_mb(k_hole+int_mb(k_offset+t_h2b-1)+t_h2-1)
     &                +dbl_mb(k_hole+int_mb(k_offset+t_h3b-1)+t_h3-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h1b-1)+t_h1-1)
     &            +int_mb(k_offset+t_h2b-1)+t_h2-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h1b-1)+t_h1-1)
     &            +int_mb(k_offset+t_h3b-1)+t_h3-1)
     &     +dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+t_h2b-1)+t_h2-1)
     &            +int_mb(k_offset+t_h3b-1)+t_h3-1))
                 enddo
                enddo
               enddo
              enddo
             enddo
            enddo
            if (.not.MA_POP_STACK(l_den)) 
     1        call errquit('ccsd_t',6,MA_ERR)
            if (.not.MA_POP_STACK(l_right)) 
     1        call errquit('ccsd_t',6,MA_ERR)
            if (.not.MA_POP_STACK(l_doubles)) 
     1        call errquit('ccsd_t',7,MA_ERR)
cccx            if (.not.MA_POP_STACK(l_singles)) 
cccx     1        call errquit('ccsd_t',8,MA_ERR)
c
            next = nxtask(nprocs,1)
            endif
            count = count + 1
c
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
      call cr_ccsd_t_E(dbl_mb(k_den),d_i1_3,
     1  k_t1_local,d_t2,k_i1_offset_3,k_t1_offset,k_t2_offset,
     2  l_i1_offset_3,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,3)
      call cr_ccsd_t_N(dbl_mb(k_right),d_f1,d_i1_1,d_i1_2,
     1  k_t1_local,d_t2,d_v2,k_f1_offset,k_i1_offset_1,k_i1_offset_2,
     2  k_t1_offset,k_t2_offset,k_v2_offset,l_i1_offset_1,
     3  l_i1_offset_2,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,3)
      call ga_zero(g_energy)
      call ga_acc(g_energy,1,1,1,1,num1,1,1.0d0)
      call ga_sync()
      call ga_get(g_energy,1,1,1,1,num1,1)
      call ga_zero(g_energy)
      call ga_acc(g_energy,1,1,1,1,num2,1,1.0d0)
      call ga_sync()
      call ga_get(g_energy,1,1,1,1,num2,1)
      call ga_zero(g_energy)
      call ga_acc(g_energy,1,1,1,1,den0,1,1.0d0)
      call ga_sync()
      call ga_get(g_energy,1,1,1,1,den0,1)
      call ga_zero(g_energy)
      call ga_acc(g_energy,1,1,1,1,den1,1,1.0d0)
      call ga_sync()
      call ga_get(g_energy,1,1,1,1,den1,1)
      call ga_zero(g_energy)
      call ga_acc(g_energy,1,1,1,1,den2,1,1.0d0)
      call ga_sync()
      call ga_get(g_energy,1,1,1,1,den2,1)
      call ga_zero(g_energy)
      call ga_acc(g_energy,1,1,1,1,den3,1,1.0d0)
      call ga_sync()
      call ga_get(g_energy,1,1,1,1,den3,1)
c
      if (.not.ga_destroy(g_energy))
     1  call errquit('ccsd_t: GA problem',1,GA_ERR)
      energy1 = num1
      energy2 = num2
      energy3 = den0
      energy4 = den1
      energy5 = den2
      energy6 = den3
c
      if (.not.ma_pop_stack(l_2hole))
     1  call errquit("tce_energy: MA problem",19,MA_ERR)
      if (.not.ma_pop_stack(l_hole))
     1  call errquit("tce_energy: MA problem",19,MA_ERR)
c
c - T1/X1 LOCALIZATION ------
         if(.not.MA_POP_STACK(l_t1_local))
     &      call errquit('l_t1_local',4,MA_ERR)
c ---------------------------
      return
      end
c
c
c
c
c
c
      subroutine tce_nu1(k_hole,d_t1,k_t1_offset)
c
c $Id$
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer d_t1
      integer p1b
      integer h2b
      integer p1
      integer h2
      integer k_t1_offset
      integer size
      integer l_t1,k_t1
      integer i
      integer k_hole
      integer pp1b,hh2b
      integer spinsum
c
c     ================
c     Loop over blocks
c     ================
c
      do p1b = noab+1,noab+nvab
        do h2b = 1,noab
            if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h2b-1)) then
            SPINSUM=int_mb(k_spin+p1b-1)+int_mb(k_spin+h2b-1)
            if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h2b-1)).eq.0) then
              size = int_mb(k_range+p1b-1) * int_mb(k_range+h2b-1)
              if (.not.ma_push_get(mt_dbl,size,'t1',l_t1,k_t1))
     1          call errquit('tce_guess_t1: MA problem',0,MA_ERR)
              IF(.NOT.RESTRICTED.OR.SPINSUM.NE.4) THEN
              call get_hash_block(d_t1,dbl_mb(k_t1),size,
     1          int_mb(k_t1_offset),((p1b-noab-1)*noab+h2b-1))
              i = 0
              do p1 = 1,int_mb(k_range+p1b-1)
                do h2 = 1,int_mb(k_range+h2b-1)
                  i = i + 1
                  dbl_mb(k_hole+int_mb(k_offset+h2b-1)+h2-1) =
     1            dbl_mb(k_hole+int_mb(k_offset+h2b-1)+h2-1) +
     1            dbl_mb(k_t1+i-1)*dbl_mb(k_t1+i-1)
                enddo
              enddo
              if (.not.ma_pop_stack(l_t1))
     1          call errquit('tce_guess_t1: MA problem',1,MA_ERR)
              END IF
              IF(RESTRICTED.AND.SPINSUM.EQ.4) THEN  ! 2beta+restr
                pp1b=int_mb(k_alpha+p1b-1)
                hh2b=int_mb(k_alpha+h2b-1)
              call get_hash_block(d_t1,dbl_mb(k_t1),size,
     1          int_mb(k_t1_offset),((pp1b-noab-1)*noab+hh2b-1))
              i = 0
              do p1 = 1,int_mb(k_range+p1b-1)
                do h2 = 1,int_mb(k_range+h2b-1)
                  i = i + 1
                  dbl_mb(k_hole+int_mb(k_offset+h2b-1)+h2-1) =
     1            dbl_mb(k_hole+int_mb(k_offset+h2b-1)+h2-1) +
     1            dbl_mb(k_t1+i-1)*dbl_mb(k_t1+i-1)
                enddo
              enddo
              if (.not.ma_pop_stack(l_t1))
     1          call errquit('tce_guess_t1: MA problem',1,MA_ERR)
              END IF
            endif
            endif
        enddo
      enddo
      return
      end
c
c
c
c
      subroutine tce_mu2(k_2hole,hole_p_1,d_t2,k_t2_offset)
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
      integer d_t2
      integer p1b
      integer p2b
      integer h3b
      integer h4b
      integer p1
      integer p2
      integer h3
      integer h4
      integer k_2hole,hole_p_1
      integer k_t2_offset
      integer size
      integer l_t2,k_t2
      integer i
      integer ipa1,ipa2,iha3,iha4
      integer spinsum
      integer pp1b,pp2b,hh3b,hh4b
c
c     ================
c     Loop over blocks
c     ================
c
      do p1b = noab+1,noab+nvab
        do p2b = p1b,noab+nvab
          do h3b = 1,noab
            do h4b = h3b,noab
                if (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)
     1            .eq. int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) then
                SPINSUM = int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+
     1                    int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)
                if (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),
     1            ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1))))
     2            .eq. 0) then
                 size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
     1                * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
                 if (.not.ma_push_get(mt_dbl,size,'t2',l_t2,k_t2))
     1             call errquit('tce_guess_t2: MA problem',0,MA_ERR)
                  IF(.NOT.RESTRICTED.OR.SPINSUM.NE.8) THEN
                  call get_hash_block(d_t2,dbl_mb(k_t2),size,
     1              int_mb(k_t2_offset),((((p1b-noab-1)*nvab+p2b-noab-1)
     2              *noab+h3b-1)*noab+h4b-1))
                 i = 0
                 do p1 = 1,int_mb(k_range+p1b-1)
                   do p2 = 1,int_mb(k_range+p2b-1)
                     do h3 = 1,int_mb(k_range+h3b-1)
                       do h4 = 1,int_mb(k_range+h4b-1)
                         i = i + 1
          ipa1=int_mb(k_offset+p1b-1)+p1
          ipa2=int_mb(k_offset+p2b-1)+p2
          iha3=int_mb(k_offset+h3b-1)+h3
          iha4=int_mb(k_offset+h4b-1)+h4
          if(ipa1.lt.ipa2.AND.iha3.lt.iha4) THEN
          dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+h3b-1)+h3-1)
     1          +int_mb(k_offset+h4b-1)+h4-1) =
     2    dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+h3b-1)+h3-1)
     3          +int_mb(k_offset+h4b-1)+h4-1)
     4          +dbl_mb(k_t2+i-1)*dbl_mb(k_t2+i-1)
          dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+h4b-1)+h4-1)
     1          +int_mb(k_offset+h3b-1)+h3-1) =
     2    dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+h4b-1)+h4-1)
     3          +int_mb(k_offset+h3b-1)+h3-1) +
     4           dbl_mb(k_t2+i-1)*dbl_mb(k_t2+i-1)
         end if
                       enddo
                     enddo
                   enddo
                 enddo
                 if (.not.ma_pop_stack(l_t2))
     1             call errquit('tce_guess_t2: MA problem',1,MA_ERR)
                 END IF
                 IF(RESTRICTED.AND.SPINSUM.EQ.8) THEN  ! 4beta+restr
                  pp1b=int_mb(k_alpha+p1b-1)
                  pp2b=int_mb(k_alpha+p2b-1)
                  hh3b=int_mb(k_alpha+h3b-1)
                  hh4b=int_mb(k_alpha+h4b-1)
                call get_hash_block(d_t2,dbl_mb(k_t2),size,
     1          int_mb(k_t2_offset),((((pp1b-noab-1)*nvab+pp2b-noab-1)
     2         *noab+hh3b-1)*noab+hh4b-1))
                 i = 0
                 do p1 = 1,int_mb(k_range+p1b-1)
                   do p2 = 1,int_mb(k_range+p2b-1)
                     do h3 = 1,int_mb(k_range+h3b-1)
                       do h4 = 1,int_mb(k_range+h4b-1)
                         i = i + 1
          ipa1=int_mb(k_offset+p1b-1)+p1
          ipa2=int_mb(k_offset+p2b-1)+p2
          iha3=int_mb(k_offset+h3b-1)+h3
          iha4=int_mb(k_offset+h4b-1)+h4
          if(ipa1.lt.ipa2.AND.iha3.lt.iha4) THEN
          dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+h3b-1)+h3-1)
     1          +int_mb(k_offset+h4b-1)+h4-1) =
     2    dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+h3b-1)+h3-1)
     3          +int_mb(k_offset+h4b-1)+h4-1)
     4          +dbl_mb(k_t2+i-1)*dbl_mb(k_t2+i-1)
          dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+h4b-1)+h4-1)
     1          +int_mb(k_offset+h3b-1)+h3-1) =
     2    dbl_mb(k_2hole+hole_p_1*(int_mb(k_offset+h4b-1)+h4-1)
     3          +int_mb(k_offset+h3b-1)+h3-1) +
     4           dbl_mb(k_t2+i-1)*dbl_mb(k_t2+i-1)
         end if
                       enddo
                     enddo
                   enddo
                 enddo
                 if (.not.ma_pop_stack(l_t2))
     1             call errquit('tce_guess_t2: MA problem',1,MA_ERR)
                  END IF   ! 4beta+restr
                endif
                endif
            enddo
          enddo
        enddo
      enddo
      return
      end

