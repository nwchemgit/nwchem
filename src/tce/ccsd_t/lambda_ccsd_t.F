      SUBROUTINE lambda_ccsd_t(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                  d_y1,k_y1_offset,d_y2,k_y2_offset,
     2                  d_f1,k_f1_offset,d_v2,k_v2_offset,
     3                  energy1,energy2)
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
      integer d_y1
      integer k_y1_offset
      integer d_y2
      integer k_y2_offset
      integer d_f1
      integer k_f1_offset
      integer d_v2
      integer k_v2_offset
      integer t_h1b, t_h1
      integer t_h2b, t_h2
      integer t_h3b, t_h3
      integer t_p4b, t_p4
      integer t_p5b, t_p5
      integer t_p6b, t_p6
      integer k_tdoubles,l_tdoubles
      integer k_ysingles,l_ysingles
      integer k_ydoublesA,l_ydoublesA ! L3 form, before sorting
      integer k_ydoubles,l_ydoubles   ! T3 form, after sorting
      integer size,i
      integer g_energy
      integer nxtask
      integer next
      integer nprocs
      integer count
      double precision energy1,energy2
      double precision factor
      external nxtask
C
#ifdef DEBUG_PRINT
      print*,"entering lambda_ccsd_t"
#endif
      if (.not.ga_create(mt_dbl,1,1,'perturbative',1,1,g_energy))
     1  call errquit('lambda_ccsd_t: GA problem',0,GA_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = nxtask(nprocs,1)
      energy1=0.0d0
      energy2=0.0d0
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
c
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) right doubles',
     1        l_tdoubles,k_tdoubles)) 
     2        call errquit('lambda_ccsd_t: MA error',301,MA_ERR)
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) left singles',
     1        l_ysingles,k_ysingles)) 
     2        call errquit('lambda_ccsd_t: MA error',302,MA_ERR)
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) left doubles',
     1        l_ydoubles,k_ydoubles)) 
     2        call errquit('lambda_ccsd_t: MA error',303,MA_ERR)
c
            do i = 1, size
             dbl_mb(k_tdoubles+i-1) = 0.0d0
             dbl_mb(k_ysingles+i-1) = 0.0d0
             dbl_mb(k_ydoubles+i-1) = 0.0d0
            enddo
c
#ifdef DEBUG_PRINT
      print*,"before ccsd_t_doubles"
#endif
c
c           <T|{W*T2}|0> : i0(p4 p5 p6 h1 h2 h3)_vt += -P(9)*Sum(h7)*t(p4 p5 h1 h7)_t*v(h7 p6 h2 h3)_v
c                        : i0(p4 p5 p6 h1 h2 h3)_vt += -P(9)*Sum(p7)*t(p4 p7 h1 h2)_t*v(p5 p6 h3 p7)_v
c
            call ccsd_t_doubles(dbl_mb(k_tdoubles),d_t2,d_v2,
     1           k_t2_offset,k_v2_offset,
     2           t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
c
#ifdef DEBUG_PRINT
      print*,"before lambda_ccsd_t_left"
#endif
c
c           <0|L(F+W)|T> : i0(h4 h5 h6 p1 p2 p3)_yv += P(9)*y(h4 p1)_y*v(h5 h6 p2 p3)_v
c
            call lambda_ccsd_t_left(dbl_mb(k_ysingles),d_f1,d_v2,
     1           d_y1,d_y2,k_f1_offset,k_v2_offset,
     2           k_y1_offset,k_y2_offset,
     3           t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,1)
c
c           <0|L(F+W)|T> : i0(h4 h5 h6 p1 p2 p3)_yf += P(9)*y(h4 h5 p1 p2)_y*f(h6 p3)_f
c                        : i0(h4 h5 h6 p1 p2 p3)_yv += -P(9)*Sum(h7)*y(h4 h7 p1 p2)_y*v(h5 h6 h7 p3)_v
c                        : i0(h4 h5 h6 p1 p2 p3)_yv += -P(9)*Sum(p7)*y(h4 h5 p1 p7)_y*v(h6 p7 p2 p3)_v
c
            call lambda_ccsd_t_left(dbl_mb(k_ydoubles),d_f1,d_v2,
     1           d_y1,d_y2,k_f1_offset,k_v2_offset,
     2           k_y1_offset,k_y2_offset,
     3           t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
c
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
                  energy1 = energy1 + factor * dbl_mb(k_tdoubles+i-1)
     1                                       * dbl_mb(k_ydoubles+i-1)
     2      / (-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p4b-1)+t_p4-1)
     3         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     4         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     5         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     6         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     7         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1))
                  energy2 = energy2 + factor * dbl_mb(k_tdoubles+i-1)
     1             * (dbl_mb(k_ysingles+i-1) + dbl_mb(k_ydoubles+i-1))
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
            if (.not.MA_POP_STACK(l_ydoubles)) 
     1        call errquit('lambda_ccsd_t',314,MA_ERR)
            if (.not.MA_POP_STACK(l_ysingles)) 
     1        call errquit('lambda_ccsd_t',312,MA_ERR)
            if (.not.MA_POP_STACK(l_tdoubles)) 
     1        call errquit('lambda_ccsd_t',311,MA_ERR)
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
      call ga_zero(g_energy)
      call ga_acc(g_energy,1,1,1,1,energy1,1,1.0d0)
      call ga_sync()
      call ga_get(g_energy,1,1,1,1,energy1,1)
      call ga_zero(g_energy)
      call ga_acc(g_energy,1,1,1,1,energy2,1,1.0d0)
      call ga_sync()
      call ga_get(g_energy,1,1,1,1,energy2,1)
      if (.not.ga_destroy(g_energy))
     1  call errquit('lambda_ccsd_t: GA problem',1,GA_ERR)
      return
      end
