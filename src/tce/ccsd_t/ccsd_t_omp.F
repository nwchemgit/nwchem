      SUBROUTINE ccsd_t_omp(d_t1,k_t1_offset,
     &                      d_t2,k_t2_offset,
     &                      d_v2,k_v2_offset,
     &                      energy1,energy2,size_t1)
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "hybrid.fh"
#include "stdio.fh"
      integer, intent(in) :: d_t1, d_t2, d_v2
      integer, intent(in) :: k_t1_offset, k_t2_offset, k_v2_offset
      integer t_h1b, t_h1
      integer t_h2b, t_h2
      integer t_h3b, t_h3
      integer t_p4b, t_p4
      integer t_p5b, t_p5
      integer t_p6b, t_p6
      integer k_singles,l_singles
      integer k_doubles,l_doubles
      integer size,i
      integer nxtask
      integer next
      integer nprocs
      integer count
      integer l_t1_local,k_t1_local
      integer size_t1
      double precision energy1,energy2
      double precision factor
      external nxtask
      logical nodezero
#if 1
      integer range_p4,range_p5,range_p6
      integer range_h1,range_h2,range_h3
      integer offset_p4,offset_p5,offset_p6
      integer offset_h1,offset_h2,offset_h3
#endif
      nodezero=(ga_nodeid().eq.0)
      if (.not.MA_PUSH_GET(mt_dbl,size_t1,'t1_local',
     1                     l_t1_local,k_t1_local))
     1    call errquit('t1_local',1,MA_ERR)
      call ma_zero(dbl_mb(k_t1_local),size_t1)
      call ga_get(d_t1,1,size_t1,1,1,dbl_mb(k_t1_local),size_t1)
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
            if (next.eq.count) then
            size = int_mb(k_range+t_p4b-1)
     1           * int_mb(k_range+t_p5b-1)
     2           * int_mb(k_range+t_p6b-1)
     3           * int_mb(k_range+t_h1b-1)
     4           * int_mb(k_range+t_h2b-1)
     5           * int_mb(k_range+t_h3b-1)
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) singles',
     &                           l_singles,k_singles))
     &          call errquit('ccsd_t: MA error',1,MA_ERR)
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) doubles',
     &                           l_doubles,k_doubles))
     &          call errquit('ccsd_t: MA error',2,MA_ERR)
            call dfill(size, 0.0d0, dbl_mb(k_singles), 1)
            call dfill(size, 0.0d0, dbl_mb(k_doubles), 1)
            call ccsd_t_singles_omp(dbl_mb(k_singles),
     1        k_t1_local,d_v2,k_t1_offset,k_v2_offset,
     2        t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
            call ccsd_t_doubles_omp(dbl_mb(k_doubles),d_t2,d_v2,
     1        k_t2_offset,k_v2_offset,
     2        t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
#if 1
            range_p4  = int_mb(k_range+t_p4b-1)
            range_p5  = int_mb(k_range+t_p5b-1)
            range_p6  = int_mb(k_range+t_p6b-1)
            range_h1  = int_mb(k_range+t_h1b-1)
            range_h2  = int_mb(k_range+t_h2b-1)
            range_h3  = int_mb(k_range+t_h3b-1)
            offset_p4 = k_evl_sorted+int_mb(k_offset+t_p4b-1)+t_p4
            offset_p5 = k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5
            offset_p6 = k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6
            offset_h1 = k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1
            offset_h2 = k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2
            offset_h3 = k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3
            call ccsd_t_dot(dbl_mb(k_singles),dbl_mb(k_doubles),
     &                      restricted,
     &                      t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,
     &                      dbl_mb(offset_h1),dbl_mb(offset_h2),
     &                      dbl_mb(offset_h3),dbl_mb(offset_p4),
     &                      dbl_mb(offset_p5),dbl_mb(offset_p6),
     &                      range_h1,range_h2,range_h3,
     &                      range_p4,range_p5,range_p6,
     &                      energy1,energy2)
#else
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
                  energy1 = energy1 + factor * dbl_mb(k_doubles+i-1)
     1                                       * dbl_mb(k_doubles+i-1)
     2      / (-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p4b-1)+t_p4-1)
     3         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     4         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     5         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     6         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     7         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1))
                  energy2 =  energy2 + factor * dbl_mb(k_doubles+i-1)
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
#endif
            if (.not.MA_POP_STACK(l_doubles))
     1        call errquit('ccsd_t_omp',3,MA_ERR)
            if (.not.MA_POP_STACK(l_singles))
     1        call errquit('ccsd_t_omp',4,MA_ERR)
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
      call ga_dgop(1975,energy1,1,'+')
      call ga_dgop(1976,energy2,1,'+')
      if(.not.MA_POP_STACK(l_t1_local))
     &   call errquit('l_t1_local',4,MA_ERR)
      return
      end
