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
      integer g_total_d, g_total_s
      integer nxtask
      integer next
      integer nprocs
      integer count
      integer armci_master
c - T1/X1 LOCALIZATION -------------------
      integer l_t1_local,k_t1_local
      integer size_t1
c ---------------------------------------

      double precision energy1,energy2,energy2_t
      double precision factor
      double precision factor_l(1)
      double precision energy_l(2),total_d(1),total_s(1),total_all_d
      double precision total_all_s
      external nxtask
c Wenjing
c for getting device information
      external integer armci_master
      external device_init
      double precision time1
ckbn -2
c      double precision sum_s,sum_d
c    NEW...
      integer has_GPU
      external check_device
      logical nodezero
c      static int device_id=-1
c
c - T1/X1 LOCALIZATION ----------
c    opening l_t1_local and l_x1_local
c     NEW...

cTCE_CUDA 
      integer icuda
      integer cuda_device_number
      cuda_device_number = 0

ckbn -2
      nodezero=(ga_nodeid().eq.0)
ckbn      sum_s = 0.0d0
ckbn      sum_d = 0.0d0
      has_GPU = check_device(icuda)
      if (has_GPU.eq.1) then
       call device_init(icuda,cuda_device_number)
       if(cuda_device_number .eq. 30 ) call errquit("cuda",30,INPUT_ERR)
      endif
      if(nodezero) 
     + write(*,'(A,I3,A)') "Using ",icuda, " device per node"
      if (nodezero) call util_flush(LuOut)

ckbn In a large cluster it is better to get and broadcast
        if (.not.MA_PUSH_GET(mt_dbl,size_t1,'t1_local',
     1      l_t1_local,k_t1_local))
     1      call errquit('t1_local',1,MA_ERR)
        call ma_zero(dbl_mb(k_t1_local),size_t1)
c    copy d_t1 ==> l_t1_local
        call ga_get(d_t1,1,size_t1,1,1,dbl_mb(k_t1_local),size_t1)
c -------------------------------
c


c      if (.not.ga_create(mt_dbl,1,1,'total_d',1,1,g_total_d))
c     1  call errquit('ccsd_t: GA problem',0,GA_ERR)
c      if (.not.ga_create(mt_dbl,1,1,'total_s',1,1,g_total_s))
c     1  call errquit('ccsd_t: GA problem',0,GA_ERR)

      nprocs = GA_NNODES()
      count = 0
      next = nxtask(nprocs,1)
c      total_all_d = 0.0d0
c      total_all_s = 0.0d0

      energy1=0.0d0
      energy2=0.0d0
      do t_p4b = noab+1,noab+nvab
       do t_p5b = t_p4b,noab+nvab
        do t_p6b = t_p5b,noab+nvab
         do t_h1b = 1,noab
          do t_h2b = t_h1b,noab
           do t_h3b = t_h2b,noab
ccx            if (next.eq.count) then                        


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
            time1 = - util_wallsec()
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) singles',l_singles,
     1        k_singles)) call errquit('ccsd_t: MA error',1,MA_ERR)

            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) doubles',l_doubles,
     1        k_doubles)) call errquit('ccsd_t: MA error',2,MA_ERR)
            time1=time1+ util_wallsec()
c            write (*,*) 'time for MA_PUSH_GET ', time1
ccx            do i = 1, size
ccx             dbl_mb(k_singles+i-1) = 0.0d0
ccx             dbl_mb(k_doubles+i-1) = 0.0d0
ccx            enddo
c zeroing ---
            time1 = - util_wallsec()
        call dfill(size, 0.0d0, dbl_mb(k_singles), 1)


            time1=time1+ util_wallsec()
c            write (*,*) 'time for dfill MA_PUSH_GET ', time1
        call dfill(size, 0.0d0, dbl_mb(k_doubles), 1)
c -----------
c      call device_init()
c            device_me = get_device_id()
c        NEW..
c      init GPU mem
            has_GPU = check_device(icuda)
            if (has_GPU.eq.1) then
              call initmemmodule()
            endif

c       NEW
            has_GPU = check_device(icuda)
            call ccsd_t_singles_gpu(dbl_mb(k_singles),
     1        k_t1_local,d_v2,k_t1_offset,
     2        k_v2_offset,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2,
     3        has_GPU)

c            device_me = get_device_id()
c            if (device_me<NUM_DEVICE) then
            call ccsd_t_doubles_gpu(dbl_mb(k_doubles),d_t2,d_v2,
     +        k_t2_offset,
     1        k_v2_offset,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2, 
     2        has_GPU)
c            else
c            call ccsd_t_doubles(dbl_mb(k_doubles),d_t2,d_v2,k_t2_offset,
c     1        k_v2_offset,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2, 0)
c            endif

ckbn            do i=1,size
ckbn             sum_s = sum_s + dbl_mb(k_singles+i-1)
ckbn             sum_d = sum_d + dbl_mb(k_doubles+i-1)
ckbn            enddo

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
c             write (*,*) 't_p4 is ', int_mb(k_offset+t_p4b-1),
c     1       int_mb(k_offset+t_p5b-1),int_mb(k_offset+t_p6b-1),
c     2       int_mb(k_offset+t_h1b-1),int_mb(k_offset+t_h2b-1),
c     3       int_mb(k_offset+t_h3b-1)


c       NEW
c            device_me = get_device_id()
            has_GPU = check_device(icuda)
            if (has_GPU.eq.0) then
c      CPU process
            do t_p4 = 1, int_mb(k_range+t_p4b-1)
             do t_p5 = 1, int_mb(k_range+t_p5b-1)
              do t_p6 = 1, int_mb(k_range+t_p6b-1)
               do t_h1 = 1, int_mb(k_range+t_h1b-1)
                do t_h2 = 1, int_mb(k_range+t_h2b-1)
                 do t_h3 = 1, int_mb(k_range+t_h3b-1)
                  i = i + 1
c                  energy1 = energy1 +  dbl_mb(k_singles+i-1)

                  energy1 = energy1 + factor * dbl_mb(k_doubles+i-1)
     1                                       * dbl_mb(k_doubles+i-1)
     2      / (-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p4b-1)+t_p4-1)
     3         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     4         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     5         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     6         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     7         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1))
c                  energy2_t = factor * dbl_mb(k_doubles+i-1)
c                  energy2_t = dbl_mb(k_doubles+i-1)
c                  energy2 =  energy2 + factor * dbl_mb(k_doubles+i-1)
                  energy2 =  energy2 + factor * dbl_mb(k_doubles+i-1)
     1             * (dbl_mb(k_singles+i-1) + dbl_mb(k_doubles+i-1))
     2      / (-dbl_mb(k_evl_sorted+int_mb(k_offset+t_p4b-1)+t_p4-1)
     3         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p5b-1)+t_p5-1)
     4         -dbl_mb(k_evl_sorted+int_mb(k_offset+t_p6b-1)+t_p6-1)
     5         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)+t_h1-1)
     6         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h2b-1)+t_h2-1)
     7         +dbl_mb(k_evl_sorted+int_mb(k_offset+t_h3b-1)+t_h3-1))
c                  write (*,*) 'CPU energy ', energy2, energy2_t
c                  energy2 = energy2+energy2_t
c                  total_all_d = total_all_d + dbl_mb(k_doubles+1)
c     dbl_mb(k_evl_sorted+int_mb
c     &(k_offset+t_p4b-1)+t_p4-1)
c     dbl_mb(k_doubles+i-1)
c                  total_all_s = total_all_s + dbl_mb(k_singles+i-1)
                enddo
               enddo
              enddo
              enddo
             enddo
            enddo
c            total_all_d = total_all_d + dbl_mb(k_doubles+1)
c            total_all_s = total_all_s + dbl_mb(k_singles+1)

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
c          write (*,*) 'F paras ',dbl_mb(k_doubles),dbl_mb(k_singles+1),
c     1        dbl_mb(k_evl_sorted+int_mb(k_offset+t_h1b-1)),factor_l(1)
c              write (*,*) 'GPU energy ', energy1, energy_l(1), energy2, 
c     1             energy_l(2)
c              write (*,*) 'double and single', total_d(1), total_s(1)
c              total_all_d = total_all_d + total_d(1)
c              total_all_s = total_all_s + total_s(1)
              energy1 = energy1 + energy_l(1)
              energy2 = energy2 + energy_l(2)
c    release GPU memory
              call dev_release()
              call finalizememmodule()
c            endif
           endif

c
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
ccx            next = nxtask(nprocs,1)
ccx            endif
ccx            count = count + 1
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      next = nxtask(-nprocs,1)
#if 0
      if (.not.ga_create(mt_dbl,1,1,'perturbative',1,1,g_energy))
     1  call errquit('ccsd_t: GA problem',0,GA_ERR)
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
#else
      call ga_dgop(1975,energy1,1,'+')
      call ga_dgop(1976,energy2,1,'+')
#endif



c - T1/X1 LOCALIZATION ------
         if(.not.MA_POP_STACK(l_t1_local))
     &      call errquit('l_t1_local',4,MA_ERR)
c ---------------------------
      return
      end
c $Id$
