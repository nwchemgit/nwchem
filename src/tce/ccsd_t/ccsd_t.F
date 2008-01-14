      SUBROUTINE ccsd_t(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                  d_v2,k_v2_offset,energy1,energy2)
C
C     $Id: ccsd_t.F,v 1.2 2008-01-14 15:20:39 jhammond Exp $
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
      integer t_h1b, t_h1
      integer t_h2b, t_h2
      integer t_h3b, t_h3
      integer t_p4b, t_p4
      integer t_p5b, t_p5
      integer t_p6b, t_p6
      integer k_singles,l_singles
      integer k_doubles,l_doubles
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
      if (.not.ga_create(mt_dbl,1,1,'perturbative',1,1,g_energy))
     1  call errquit('ccsd_t: GA problem',0,GA_ERR)
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
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) singles',
     1        l_singles,k_singles)) call errquit('ccsd_t',1,MA_ERR)
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) doubles',
     1        l_doubles,k_doubles)) call errquit('ccsd_t',2,MA_ERR)
            do i = 1, size
             dbl_mb(k_singles+i-1) = 0.0d0
             dbl_mb(k_doubles+i-1) = 0.0d0
            enddo
            call ccsd_t_singles(dbl_mb(k_singles),d_t1,d_v2,k_t1_offset,
     1        k_v2_offset,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
            call ccsd_t_doubles(dbl_mb(k_doubles),d_t2,d_v2,k_t2_offset,
     1        k_v2_offset,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
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
                  energy2 = energy2 + factor * dbl_mb(k_doubles+i-1)
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
            if (.not.MA_POP_STACK(l_doubles)) 
     1        call errquit('ccsd_t',3,MA_ERR)
            if (.not.MA_POP_STACK(l_singles)) 
     1        call errquit('ccsd_t',4,MA_ERR)
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
     1  call errquit('ccsd_t: GA problem',1,GA_ERR)
      return
      end
