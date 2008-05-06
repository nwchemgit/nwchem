      SUBROUTINE lambda_ccsd_t(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                  d_y1,k_y1_offset,d_y2,k_y2_offset,
     2                  d_f1,k_f1_offset,d_v2,k_v2_offset,
     3                  energy1,energy2)
C
C     $Id: lambda_ccsd_t.F,v 1.1 2008-05-06 16:20:59 jhammond Exp $
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
c      integer k_tsingles,l_tsingles
      integer k_tdoubles,l_tdoubles
      integer k_ysingles,l_ysingles
      integer k_ydoubles,l_ydoubles
      integer size,i
      integer g_energy
      integer NXTVAL
      integer next
      integer nprocs
      integer count
      double precision energy1,energy2
      double precision factor
      external NXTVAL
C
#ifdef DEBUG_PRINT
      print*,"entering lambda_ccsd_t"
#endif
      if (.not.ga_create(mt_dbl,1,1,'perturbative',1,1,g_energy))
     1  call errquit('lambda_ccsd_t: GA problem',0,GA_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTVAL(nprocs)
      energy1=0.0d0
      energy2=0.0d0
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
c
c            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) singles',l_tsingles,
c     1        k_tsingles)) 
c     2        call errquit('lambda_ccsd_t: MA error',11,MA_ERR)
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) doubles',l_tdoubles,
     1        k_tdoubles)) 
     2        call errquit('lambda_ccsd_t: MA error',12,MA_ERR)
            if (.not.MA_PUSH_GET(mt_dbl,size,'Lambda 1',l_ysingles,
     1        k_ysingles)) 
     2        call errquit('lambda_ccsd_t: MA error',21,MA_ERR)
            if (.not.MA_PUSH_GET(mt_dbl,size,'Lambda 2',l_ydoubles,
     1        k_ydoubles)) 
     2        call errquit('lambda_ccsd_t: MA error',22,MA_ERR)
c
c
            do i = 1, size
c             dbl_mb(k_tsingles+i-1) = 0.0d0
             dbl_mb(k_tdoubles+i-1) = 0.0d0
             dbl_mb(k_ysingles+i-1) = 0.0d0
             dbl_mb(k_ydoubles+i-1) = 0.0d0
            enddo
c
c            call ccsd_t_singles(dbl_mb(k_tsingles),d_t1,d_v2,
c     1           k_t1_offset,k_v2_offset,
c     2           t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
c
#ifdef DEBUG_PRINT
      print*,"before ccsd_t_doubles"
#endif
            call ccsd_t_doubles(dbl_mb(k_tdoubles),d_t2,d_v2,
     1           k_t2_offset,k_v2_offset,
     2           t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
c
#ifdef DEBUG_PRINT
      print*,"before lambda_ccsd_t_singles"
#endif
            call lambda_ccsd_t_singles(dbl_mb(k_ysingles),
     1           d_f1,d_v2,d_y1,d_y2,k_f1_offset,
     2           k_v2_offset,k_y1_offset,k_y2_offset,
     3           t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
c
#ifdef DEBUG_PRINT
      print*,"before lambda_ccsd_t_doubles"
#endif
            call lambda_ccsd_t_doubles(dbl_mb(k_ydoubles),
     1           d_y2,d_v2,k_y2_offset,k_v2_offset,
     2           t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
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
     1        call errquit('lambda_ccsd_t',32,MA_ERR)
            if (.not.MA_POP_STACK(l_ysingles)) 
     1        call errquit('lambda_ccsd_t',31,MA_ERR)
            if (.not.MA_POP_STACK(l_tdoubles)) 
     1        call errquit('lambda_ccsd_t',42,MA_ERR)
c            if (.not.MA_POP_STACK(l_tsingles)) 
c     1        call errquit('ccsd_t',41,MA_ERR)
            endif
            endif
            endif
            next = NXTVAL(nprocs)
            endif
            count = count + 1
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      next = NXTVAL(-nprocs)
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
