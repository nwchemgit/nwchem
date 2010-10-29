      SUBROUTINE cr_ccsd_t(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                     d_f1,k_f1_offset,d_v2,k_v2_offset,
     2                     d_e,k_e_offset,energy1,energy2,size_t1)
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
      integer k_singles,l_singles
      integer k_doubles,l_doubles
      integer k_den,l_den
      integer k_right,l_right
      integer size,i
      integer g_energy
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
      double precision energy1,energy2
      double precision factor
      double precision den0,den1,den2,num1,num2
      external nxtask
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
c (here t1 - on GA !)
      call cr_ccsd_t_D(d_t1,d_t2,d_e,d_t1,d_t2,k_t1_offset,
     1  k_t2_offset,k_e_offset,k_t1_offset,k_t2_offset)
      call reconcilefile(d_e,1)
      call get_block(d_e,den0,1,0)
      den1 = 0.0d0
      den2 = 0.0d0
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
      num1=0.0d0
      num2=0.0d0
      do t_p4b = noab+1,noab+nvab
       do t_p5b = t_p4b,noab+nvab
        do t_p6b = t_p5b,noab+nvab
         do t_h1b = 1,noab
          do t_h2b = t_h1b,noab
           do t_h3b = t_h2b,noab
            if (next.eq.count) then
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
            size = int_mb(k_range+t_p4b-1)
     1           * int_mb(k_range+t_p5b-1)
     2           * int_mb(k_range+t_p6b-1)
     3           * int_mb(k_range+t_h1b-1)
     4           * int_mb(k_range+t_h2b-1)
     5           * int_mb(k_range+t_h3b-1)
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) singles',
     1        l_singles,k_singles)) call errquit('ccsd_t',1,MA_ERR)
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) doubles',
     1        l_doubles,k_doubles)) call errquit('ccsd_t',2,MA_ERR)
            if (.not.MA_PUSH_GET(mt_dbl,size,'moment 2,3',
     1        l_right,k_right)) call errquit('ccsd_t',3,MA_ERR)
            if (.not.MA_PUSH_GET(mt_dbl,size,'denominator',
     1        l_den,k_den)) call errquit('ccsd_t',3,MA_ERR)
            do i = 1, size
             dbl_mb(k_singles+i-1) = 0.0d0
             dbl_mb(k_doubles+i-1) = 0.0d0
             dbl_mb(k_right+i-1) = 0.0d0
             dbl_mb(k_den+i-1) = 0.0d0
            enddo
            call ccsd_t_singles_l(dbl_mb(k_singles),k_t1_local,
     1        d_v2,k_t1_offset,
     1        k_v2_offset,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
            call ccsd_t_doubles_l(dbl_mb(k_doubles),
     1        d_t2,d_v2,k_t2_offset,
     1        k_v2_offset,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
            call cr_ccsd_t_N(dbl_mb(k_right),d_f1,d_i1_1,d_i1_2,
     1        k_t1_local,d_t2,d_v2,k_f1_offset,
     1        k_i1_offset_1,k_i1_offset_2,
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
                  num1 = num1 + factor * dbl_mb(k_right+i-1)
     1                                 * dbl_mb(k_doubles+i-1)
     2      / (-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p4b-1)+t_p4-1)
     3         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     4         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     5         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     6         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     7         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1))
                  num2 = num2 + factor * dbl_mb(k_right+i-1)
     1             * (dbl_mb(k_singles+i-1) + dbl_mb(k_doubles+i-1))
     2      / (-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p4b-1)+t_p4-1)
     3         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     4         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     5         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     6         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     7         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1))
                  den1 = den1 + factor * dbl_mb(k_den+i-1)
     1                                 * dbl_mb(k_doubles+i-1)
     2      / (-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p4b-1)+t_p4-1)
     3         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     4         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     5         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     6         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     7         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1))
                  den2 = den2 + factor * dbl_mb(k_den+i-1)
     1             * (dbl_mb(k_singles+i-1) + dbl_mb(k_doubles+i-1))
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
            if (.not.MA_POP_STACK(l_den)) 
     1        call errquit('ccsd_t',6,MA_ERR)
            if (.not.MA_POP_STACK(l_right)) 
     1        call errquit('ccsd_t',6,MA_ERR)
            if (.not.MA_POP_STACK(l_doubles)) 
     1        call errquit('ccsd_t',7,MA_ERR)
            if (.not.MA_POP_STACK(l_singles)) 
     1        call errquit('ccsd_t',8,MA_ERR)
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
      call ga_acc(g_energy,1,1,1,1,den1,1,1.0d0)
      call ga_sync()
      call ga_get(g_energy,1,1,1,1,den1,1)
      call ga_zero(g_energy)
      call ga_acc(g_energy,1,1,1,1,den2,1,1.0d0)
      call ga_sync()
      call ga_get(g_energy,1,1,1,1,den2,1)
      if (.not.ga_destroy(g_energy))
     1  call errquit('ccsd_t: GA problem',1,GA_ERR)
      den1 = den1 + den0
      den2 = den2 + den0
      energy1 = num1/(1.0d0+den1)
      energy2 = num2/(1.0d0+den2)
c - T1/X1 LOCALIZATION ------
         if(.not.MA_POP_STACK(l_t1_local))
     &      call errquit('l_t1_local',4,MA_ERR)
c ---------------------------
      return
      end
