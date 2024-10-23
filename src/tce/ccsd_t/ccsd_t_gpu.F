c#define NUM_DEVICE 1
      SUBROUTINE ccsd_t_gpu(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                  d_v2,k_v2_offset,energy1,energy2,size_t1,icuda)
C
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "hybrid.fh"
#include "stdio.fh"
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
      integer nxtask
      integer next
      integer nprocs
      integer count
c - T1/X1 LOCALIZATION -------------------
      integer l_t1_local,k_t1_local
      integer size_t1
c ---------------------------------------

      double precision energy1,energy2,energy2_t
      double precision factor
      double precision factor_l(1)
      double precision energy_l(2),total_d(1),total_s(1)
      external nxtask
      external device_init
      integer has_GPU
      external check_device
      logical nodezero

      integer icuda
      integer cuda_device_number
      cuda_device_number = 0

      nodezero=(ga_nodeid().eq.0)
      has_GPU = check_device(icuda)
      if (has_GPU.eq.1) then
       call device_init(icuda,cuda_device_number)
       if (cuda_device_number .eq. 30) then
         call errquit("cuda",30,INPUT_ERR)
       endif
      endif
      if (nodezero) then
        write(*,'(A,I3,A)') "Using ",icuda," device per node"
      endif
      if (nodezero) call util_flush(LuOut)

ckbn In a large cluster it is better to get and broadcast
        if (.not.MA_PUSH_GET(mt_dbl,size_t1,'t1_local',
     1      l_t1_local,k_t1_local))
     1      call errquit('t1_local',1,MA_ERR)
        call ma_zero(dbl_mb(k_t1_local),size_t1)
c    copy d_t1 ==> l_t1_local
! copied pattern from cr-eomccsd_t/cr_eomccsd_t.F
#if 1
        call util_mygabcast(d_t1,size_t1,1,dbl_mb(k_t1_local),size_t1)
#else
        call ga_get(d_t1,1,size_t1,1,1,dbl_mb(k_t1_local),size_t1)
#endif
c -------------------------------
c
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
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) singles',
     &                           l_singles,k_singles)) then
              call errquit('ccsd_t_gpu: MA error - singles',size,MA_ERR)
            endif
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) doubles',
     &                           l_doubles,k_doubles)) then
              call errquit('ccsd_t_gpu: MA error - doubles',size,MA_ERR)
            endif
            call dfill(size, 0.0d0, dbl_mb(k_singles), 1)
            call dfill(size, 0.0d0, dbl_mb(k_doubles), 1)

            has_GPU = check_device(icuda)
            if (has_GPU.eq.1) then
              call initmemmodule()
            endif

            has_GPU = check_device(icuda)
            call ccsd_t_singles_gpu(dbl_mb(k_singles),
     1        k_t1_local,d_v2,k_t1_offset,
     2        k_v2_offset,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2,
     3        has_GPU)

            call ccsd_t_doubles_gpu(dbl_mb(k_doubles),d_t2,d_v2,
     +        k_t2_offset,
     1        k_v2_offset,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2, 
     2        has_GPU)
c            else
c            call ccsd_t_doubles(dbl_mb(k_doubles),d_t2,d_v2,k_t2_offset,
c     1        k_v2_offset,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2, 0)
c            endif

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
            has_GPU = check_device(icuda)
            if (has_GPU.eq.0) then

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

            else
c     GPU process
              factor_l(1) = factor
              call compute_en(factor_l, energy_l,
     1                  dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)),
     2                  dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)),
     3                  dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)),
     4                  dbl_mb(k_evl_sorted+int_mb(k_offset+t_p4b-1)),
     5                  dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)),
     6                  dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)),
     7             int_mb(k_range+t_h1b-1),int_mb(k_range+t_h2b-1),
     8             int_mb(k_range+t_h3b-1),int_mb(k_range+t_p4b-1),
     9             int_mb(k_range+t_p5b-1),int_mb(k_range+t_p6b-1),
     1            dbl_mb(k_doubles),dbl_mb(k_singles),total_d,total_s)
              energy1 = energy1 + energy_l(1)
              energy2 = energy2 + energy_l(2)
c    release GPU memory
              call dev_release()
              call finalizememmodule()
           endif

            if (.not.MA_POP_STACK(l_doubles)) 
     1        call errquit('ccsd_t',3,MA_ERR)

            if (.not.MA_POP_STACK(l_singles)) 
     1        call errquit('ccsd_t',4,MA_ERR)
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
      call ga_dgop(1975,energy1,1,'+')
      call ga_dgop(1976,energy2,1,'+')

c - T1/X1 LOCALIZATION ------
         if(.not.MA_POP_STACK(l_t1_local))
     &      call errquit('l_t1_local',4,MA_ERR)
c ---------------------------
      return
      end
c $Id$
