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
      double precision, allocatable :: f_singles(:),f_doubles(:)
#ifdef USE_FASTMEM
      !dec$ attributes fastmem :: f_singles, f_doubles
#endif
      integer alloc_error
      integer size_t3,i
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
      integer range_p4,range_p5,range_p6
      integer range_h1,range_h2,range_h3
      integer offset_p4,offset_p5,offset_p6
      integer offset_h1,offset_h2,offset_h3
      nodezero=(ga_nodeid().eq.0)
      if (.not.MA_PUSH_GET(mt_dbl,size_t1,'t1_local',
     1                     l_t1_local,k_t1_local))
     1    call errquit('t1_local',1,MA_ERR)
      call ma_zero(dbl_mb(k_t1_local),size_t1)
      call ga_get(d_t1,1,size_t1,1,1,dbl_mb(k_t1_local),size_t1)
      nprocs = GA_NNODES()
      !
      range_p4=0
      do t_p4b = noab+1,noab+nvab
         range_p4 = max(range_p4,int_mb(k_range+t_p4b-1))
      enddo
      range_h1=0
      do t_h1b = 1,noab
         range_h1 = max(range_h1,int_mb(k_range+t_h1b-1))
      enddo
      size_t3=(range_p4**3)*(range_h1**3)
      !
      allocate( f_singles(1:size_t3), stat=alloc_error)
      if (alloc_error.ne.0) then
        call errquit('ccsd_t: MA error sgl',0,MA_ERR)
      endif
      allocate( f_doubles(1:size_t3), stat=alloc_error)
      if (alloc_error.ne.0) then
        call errquit('ccsd_t: MA error dbl',0,MA_ERR)
      endif
      !
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
            size_t3 = int_mb(k_range+t_p4b-1)
     1              * int_mb(k_range+t_p5b-1)
     2              * int_mb(k_range+t_p6b-1)
     3              * int_mb(k_range+t_h1b-1)
     4              * int_mb(k_range+t_h2b-1)
     5              * int_mb(k_range+t_h3b-1)
            call dfill(size_t3, 0.0d0, f_singles, 1)
            call dfill(size_t3, 0.0d0, f_doubles, 1)
            call ccsd_t_singles_omp(f_singles,
     1        k_t1_local,d_v2,k_t1_offset,k_v2_offset,
     2        t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
            call ccsd_t_doubles_omp(f_doubles,d_t2,d_v2,
     1        k_t2_offset,k_v2_offset,
     2        t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
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
            call ccsd_t_dot(f_singles,f_doubles,
     &                      restricted,
     &                      t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,
     &                      dbl_mb(offset_h1),dbl_mb(offset_h2),
     &                      dbl_mb(offset_h3),dbl_mb(offset_p4),
     &                      dbl_mb(offset_p5),dbl_mb(offset_p6),
     &                      range_h1,range_h2,range_h3,
     &                      range_p4,range_p5,range_p6,
     &                      energy1,energy2)
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
      deallocate( f_doubles, stat=alloc_error)
      if (alloc_error.ne.0) then
        call errquit('ccsd_t doubles',3,MA_ERR)
      endif
      deallocate( f_singles, stat=alloc_error)
      if (alloc_error.ne.0) then
        call errquit('ccsd_t singles',4,MA_ERR)
      endif
      return
      end
