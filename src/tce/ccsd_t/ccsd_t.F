c     Modifications (c) 2021 Advanced Micro Devices, Inc. All Rights Reserved.

C #ifndef OFFLOAD_DEBUG_DATA
C #define OFFLOAD_DEBUG_DATA
C #endif

#ifdef USE_OMP_TEAMS_DISTRIBUTE
#define TEAMS_DISTRIBUTE teams distribute
#else
#define TEAMS_DISTRIBUTE
#endif

#ifdef USE_OMP_SCHEDULE_STATIC_1
#define OMP_SCHEDULE schedule(static,1)
#else
#define OMP_SCHEDULE
#endif

      SUBROUTINE ccsd_t(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1            d_v2,k_v2_offset,energy1,energy2,size_t1)
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
#include "offl.fh"
#include "ccsd_t_ps.fh"
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
#ifdef USE_F90_ALLOCATABLE
      double precision, allocatable :: f_singles(:),f_doubles(:)
#ifdef USE_FASTMEM
      !dec$ attributes fastmem :: f_singles, f_doubles
#endif
      integer alloc_error
#else
      integer k_singles,l_singles
      integer k_doubles,l_doubles
#endif
      integer size,i
      integer nxtask0
      integer next
      integer nprocs
      integer count
      integer offset_p4,offset_p5,offset_p6
      integer offset_h1,offset_h2,offset_h3
      integer range_p4,range_p5,range_p6
      integer range_h1,range_h2,range_h3
c - T1/X1 LOCALIZATION -------------------
      integer l_t1_local,k_t1_local
      integer size_t1
c ---------------------------------------
      double precision energy(2)
      double precision energy1,energy2
      integer k_aux,k_list,l_aux,l_list
      integer ccsd_t_6tasks,tot_task,w_thresh,task_thresh
      integer k,iptr,k_grain
      double precision wall
      integer tasks_skip
      integer n_progr,pct_progr
      parameter(n_progr=20)
      logical i_progr(n_progr)
      external nxtask0,ccsd_t_6tasks
      logical offload_enabled
      external offload_enabled
C
c
c - T1/X1 LOCALIZATION ----------
c    opening l_t1_local and l_x1_local
        if (.not.MA_PUSH_GET(mt_dbl,size_t1,'t1_local',
     1      l_t1_local,k_t1_local))
     1      call errquit('ccsd_t: t1_local size=',size_t1,MA_ERR)
        call ma_zero(dbl_mb(k_t1_local),size_t1)
c    copy d_t1 ==> l_t1_local
ccx        call ga_get(d_t1,1,size_t1,1,1,dbl_mb(k_t1_local),1)
      call get_block(d_t1,dbl_mb(k_t1_local),size_t1,0)
c -------------------------------
C
      nprocs = GA_NNODES()
      energy(1)=0.0d0
      energy(2)=0.0d0
      energy1 = 0.0d0
      energy2 = 0.0d0
c     estimate triplesx size
      range_p4=0
      do t_p4b = noab+1,noab+nvab
         range_p4 = max(range_p4,int_mb(k_range+t_p4b-1))
      enddo
      range_h1=0
      do t_h1b = 1,noab
         range_h1 = max(range_h1,int_mb(k_range+t_h1b-1))
      enddo
      size=(range_p4**3)*(range_h1**3)

      call util_align64(size)
      triplesx_mxlgth=size
#ifdef USE_F90_ALLOCATABLE
      allocate( f_singles(1:size), stat=alloc_error)
      if (alloc_error.ne.0) then
        call errquit('ccsd_t: MA error sgl',0,MA_ERR)
      endif
      allocate( f_doubles(1:size), stat=alloc_error)
      if (alloc_error.ne.0) then
        call errquit('ccsd_t: MA error dbl',0,MA_ERR)
      endif
#ifdef TCE_OPENACC
!$acc enter data create(f_singles(1:size), f_doubles(1:size))
#endif

#if USE_OPENMP
#if USE_OFFLOAD
      if (offload_enabled()) then
!$omp target enter data
!$omp&map(alloc:f_singles)
!$omp&map(alloc:f_doubles)
#ifdef OFFLOAD_DEBUG_DATA
        write (*,*) ga_nodeid(), 'alloc f_singles', size, 'elements'
        write (*,*) ga_nodeid(), 'alloc f_doubles', size, 'elements'
#endif
      endif
#endif
#endif
#else
#if USE_OFFLOAD
      call errquit('ccsd_t: offload not supported with MA',0,0)
#endif
      if (.not.MA_PUSH_GET(mt_dbl,size,'(T) singles',l_singles,
     1     k_singles)) call errquit('ccsd_t: MA error sgl',
     2     size,MA_ERR)
      if (.not.MA_PUSH_GET(mt_dbl,size,'(T) doubles',l_doubles,
     1     k_doubles)) call errquit('ccsd_t: MA error dbl',
     2     size,MA_ERR)
#endif
      wall=-util_wallsec()
      tot_task= ccsd_t_6tasks(restricted,noab,nvab,
     1                        int_mb(k_spin),int_mb(k_sym))
      if (.not.ma_push_get(mt_int,7*tot_task,"list.task",
     1  l_list,k_list)) call errquit("k_list",1,MA_ERR)
      if (.not.ma_push_get(mt_int,7*tot_task,"auxtask",
     1  l_aux,k_aux)) call errquit("k_aux",2,MA_ERR)
c
c     get first task with weight lt ? 8
c
c      w_thresh=10
      w_thresh=0
      w_thresh=w_thresh**6
      call ccsd_t_neword(tot_task, w_thresh,task_thresh,
     R     restricted,noab,nvab,
     K     int_mb(k_spin),int_mb(k_sym),
     K     int_mb(k_range),
     A     int_mb(k_aux),int_mb(k_list))


      if (.not.MA_POP_STACK(l_aux))
     1     call errquit('ordering',3,MA_ERR)
      count = 0
      k_grain=1
      tasks_skip=0
      next = nxtask0(nprocs,k_grain,tasks_skip)
      do k=1,n_progr
         i_progr(k)=.true.
      enddo
c     stagger start of loop
      call util_mpinap(100)
      if(task_thresh.gt.1) then
      do k=1,task_thresh-1
c
         if (next.eq.count) then
#ifdef USE_F90_ALLOCATABLE
            call ccsd_t_loop(k,energy1,energy2,
     &              int_mb(k_list),int_mb(k_range),int_mb(k_offset),
     &              f_singles,f_doubles,
     &              k_t1_local,k_t1_offset,
     &              d_t2,d_v2,k_t2_offset,k_v2_offset,
     &              restricted,k_evl_sorted,size)
#else
            call ccsd_t_loop(k,energy1,energy2,
     &              int_mb(k_list),int_mb(k_range),int_mb(k_offset),
     &              dbl_mb(k_singles),dbl_mb(k_doubles),
     &              k_t1_local,k_t1_offset,
     &              d_t2,d_v2,k_t2_offset,k_v2_offset,
     &              restricted,k_evl_sorted,size)
#endif
c
            if(ga_nodeid().eq.2) then
               pct_progr=(k*n_progr)/tot_task
               if(i_progr(pct_progr)) then
                  i_progr(pct_progr)=.false.
                  write(6,'(a,i5,a,i4,a,f15.1,a,f9.1)')
     &                  '0task ',k,'  done ',
     &                  int((k*100d0)/tot_task),'%  at',
     &                  wall+util_wallsec(),' sec, (size)^1/6= ',
     &                  (size)**(1d0/6d0)
                  call util_flush(6)
               endif
            endif
            next = nxtask0(nprocs,k_grain,tasks_skip)
         endif
         count = count + 1
      enddo
      endif
      if(task_thresh.le.tot_task) then
      next = nxtask0(-nprocs,k_grain,tasks_skip)
      if(.true.) then
         next = nxtask0(nprocs,k_grain,tasks_skip)
         count=0
         do k=task_thresh,tot_task
            if (next.eq.count) then
#ifdef USE_F90_ALLOCATABLE
               call ccsd_t_loop(k,energy1,energy2,
     &                 int_mb(k_list),int_mb(k_range),int_mb(k_offset),
     &                 f_singles,f_doubles,
     &                 k_t1_local,k_t1_offset,
     &                 d_t2,d_v2,k_t2_offset,k_v2_offset,
     &                 restricted,k_evl_sorted,size)
#else
               call ccsd_t_loop(k,energy1,energy2,
     &                 int_mb(k_list),int_mb(k_range),int_mb(k_offset),
     &                 dbl_mb(k_singles),dbl_mb(k_doubles),
     &                 k_t1_local,k_t1_offset,
     &                 d_t2,d_v2,k_t2_offset,k_v2_offset,
     &                 restricted,k_evl_sorted,size)
#endif
c
               next = nxtask0(nprocs,k_grain,tasks_skip)
               if(ga_nodeid().eq.2) then
                  pct_progr=(k*n_progr)/tot_task
                  if(i_progr(pct_progr)) then
                     i_progr(pct_progr)=.false.
                     write(6,'(a,i8,a,i4,a,f15.1,a,f9.1)')
     &                     ' task',k,'  done ',
     &                     int((k*100d0)/tot_task),'%  at',
     &                     wall+util_wallsec(),' sec, (size)^1/6= ',
     &                     (size)**(1d0/6d0)
                     call util_flush(6)
                  endif
               endif
            endif
            count = count + 1
         enddo
      endif
      endif

      if (.not.MA_POP_STACK(l_list))
     1     call errquit('ordering',3,MA_ERR)
#ifdef USE_F90_ALLOCATABLE
#ifdef TCE_OPENACC
!$acc exit data delete(f_singles(1:size), f_doubles(1:size))
#endif
#if USE_OPENMP
#if USE_OFFLOAD
      if (offload_enabled()) then
!$omp target exit data
!$omp&map(release:f_singles)
!$omp&map(release:f_doubles)
#ifdef OFFLOAD_DEBUG_DATA
        write (*,*) ga_nodeid(), 'release f_singles', size, 'elements'
        write (*,*) ga_nodeid(), 'release f_doubles', size, 'elements'
#endif
      endif
#endif
#endif
      deallocate( f_doubles, stat=alloc_error)
      if (alloc_error.ne.0) then
        call errquit('ccsd_t doubles',3,MA_ERR)
      endif
      deallocate( f_singles, stat=alloc_error)
      if (alloc_error.ne.0) then
        call errquit('ccsd_t singles',4,MA_ERR)
      endif
#else
      if (.not.MA_POP_STACK(l_doubles))
     1     call errquit('ccsd_t doubles',3,MA_ERR)
      if (.not.MA_POP_STACK(l_singles))
     1     call errquit('ccsd_t singles',4,MA_ERR)
#endif
      next = nxtask0(-nprocs,k_grain,tasks_skip)
c      call ga_sync()
      energy(1) = energy1
      energy(2) = energy2
      call ga_mask_sync(.false.,.true.)
      call ga_dgop(1975,energy,2,'+')
      energy1 = energy(1)
      energy2 = energy(2)
c - T1/X1 LOCALIZATION ------
         if(.not.MA_POP_STACK(l_t1_local))
     &      call errquit('ccsd_t: l_t1_local',4,MA_ERR)
c ---------------------------
      return
      end





!     wrapper to ccsd_t_dot because of offload ugliness
      subroutine ccsd_t_esum(a_singles, a_doubles, restricted,
     &                       h1b,h2b,h3b,p4b,p5b,p6b,
     &                       o_h1,o_h2,o_h3,
     &                       o_p4,o_p5,o_p6,
     &                       r_h1,r_h2,r_h3,
     &                       r_p4,r_p5,r_p6,
     &                       energy1,energy2)
      implicit none
      integer h1b, h2b, h3b, p4b, p5b, p6b
      double precision o_h1(*),o_h2(*),o_h3(*)
      double precision o_p4(*),o_p5(*),o_p6(*)
      integer r_h1,r_h2,r_h3
      integer r_p4,r_p5,r_p6
      double precision a_singles(*)
      double precision a_doubles(*)
      logical restricted
      double precision energy1,energy2
      external offload_enabled
      logical offload_enabled
#ifdef USE_OPENMP
#ifdef USE_OFFLOAD
      if(offload_enabled()) then
         call offl_ccsd_t_dot(a_singles,a_doubles,restricted,
     &        h1b,h2b,h3b,p4b,p5b,p6b,
     &        o_h1,o_h2,o_h3,o_p4,o_p5,o_p6,
     &        r_h1,r_h2,r_h3,r_p4,r_p5,r_p6,
     &        energy1,energy2)
         else
#endif
#endif
            call ccsd_t_dot(a_singles,a_doubles, restricted,
     &        h1b,h2b,h3b,p4b,p5b,p6b,
     &        o_h1,o_h2,o_h3,o_p4,o_p5,o_p6,
     &        r_h1,r_h2,r_h3,r_p4,r_p5,r_p6,
     &        energy1,energy2)
#ifdef USE_OPENMP
#ifdef USE_OFFLOAD
         endif
#endif
#endif

      return
      end



      subroutine ccsd_t_loop(k,energy1,energy2,
     &     k_list,k_range,k_offset,a_singles,a_doubles,
     &     k_t1_local,k_t1_offset,d_t2,d_v2,k_t2_offset,k_v2_offset,
     &     restricted,k_evl_sorted,size)
      implicit none
#include "mafdecls.fh"
#include "ccsd_t_ps.fh"
      integer k
      double precision energy1,energy2
      integer k_list(7,*)
      integer k_range(*),k_offset(*)
      double precision a_singles(*),a_doubles(*)
      logical restricted
      integer k_evl_sorted
      integer k_t1_local,k_t1_offset
      integer d_t2,d_v2
      integer k_t2_offset,k_v2_offset
      integer size
#if USE_OFFLOAD
      external offload_enabled
      logical offload_enabled
#endif
c
      integer t_p4b,t_p5b,t_p6b,t_h1b,t_h2b,t_h3b
      integer range_p4,range_p5,range_p6
      integer range_h1,range_h2,range_h3
      integer offset_p4,offset_p5,offset_p6
      integer offset_h1,offset_h2,offset_h3
c
      t_p4b=k_list(1,k)
      t_p5b=k_list(2,k)
      t_p6b=k_list(3,k)
      t_h1b=k_list(4,k)
      t_h2b=k_list(5,k)
      t_h3b=k_list(6,k)
      range_p4 = k_range(t_p4b)
      range_p5 = k_range(t_p5b)
      range_p6 = k_range(t_p6b)
      range_h1 = k_range(t_h1b)
      range_h2 = k_range(t_h2b)
      range_h3 = k_range(t_h3b)
      offset_p4 = k_evl_sorted+k_offset(t_p4b)-1
      offset_p5 = k_evl_sorted+k_offset(t_p5b)-1
      offset_p6 = k_evl_sorted+k_offset(t_p6b)-1
      offset_h1 = k_evl_sorted+k_offset(t_h1b)-1
      offset_h2 = k_evl_sorted+k_offset(t_h2b)-1
      offset_h3 = k_evl_sorted+k_offset(t_h3b)-1
c
      size = range_p4 * range_p5 * range_p6
     &     * range_h1 * range_h2 * range_h3
c zeroing ---
#ifdef TCE_OPENACC
      call acc_zero(a_singles, size)
      call acc_zero(a_doubles, size)
#endif
#if USE_OPENMP
#if USE_OFFLOAD
      if (offload_enabled()) then
        call offl_zero(a_singles, size)
        call offl_zero(a_doubles, size)
      else
#endif
#endif
        call dcopy(size, 0.0d0, 0, a_singles, 1)
        call dcopy(size, 0.0d0, 0, a_doubles, 1)
#if USE_OPENMP
#if USE_OFFLOAD
      endif
#endif
#endif
c -----------
      if (otceps) call pstat_on(ps_cctsng)
#ifdef USE_OFFLOAD
      call offl_ccsd_t_singles_l(a_singles,
     1                 k_t1_local,d_v2,k_t1_offset,k_v2_offset,
     1                 t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
#else
      call ccsd_t_singles_l(a_singles,
     1                 k_t1_local,d_v2,k_t1_offset,k_v2_offset,
     1                 t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
#endif
      if (otceps) call pstat_off(ps_cctsng)


      if (otceps) call pstat_on(ps_cctdbl)      
      call ccsd_t_doubles_l(a_doubles,
     1                 d_t2,d_v2,k_t2_offset,k_v2_offset,
     1                 t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
      if (otceps) call pstat_off(ps_cctdbl)
      call ccsd_t_esum(a_singles,a_doubles,restricted,
     &                 t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,
     &                 dbl_mb(offset_h1+1),dbl_mb(offset_h2+1),
     &                 dbl_mb(offset_h3+1),dbl_mb(offset_p4+1),
     &                 dbl_mb(offset_p5+1),dbl_mb(offset_p6+1),
     &                 range_h1,range_h2,range_h3,
     &                 range_p4,range_p5,range_p6,
     &                 energy1,energy2)
      return
      end

#ifdef USE_OPENMP
#ifdef USE_OFFLOAD
      subroutine offl_zero(array, size)
      implicit none
      double precision array(size)
      integer size
      integer ii
!$omp target
!$omp TEAMS_DISTRIBUTE parallel do OMP_SCHEDULE
      do ii=1,size
         array(ii)=0d0
      enddo
!$omp end TEAMS_DISTRIBUTE parallel do
!$omp end target
      return
      end
#endif
#endif
