      SUBROUTINE ccsd_t_6dts(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1            d_v2,k_v2_offset,energy1,energy2,size_t1,xmem)
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
      integer d_v2
      integer k_v2_offset
      integer t_h1b, t_h1
      integer t_h2b, t_h2
      integer t_h3b, t_h3
      integer t_p4b, t_p4
      integer t_p5b, t_p5
      integer t_p6b, t_p6
#ifdef USE_F90_ALLOCATABLE
      double precision, allocatable :: f_singles(:), f_doubles(:)
#ifdef USE_FASTMEM
      !dec$ attributes fastmem :: f_singles, f_doubles
#endif
      integer alloc_error
#else // USE_F90_ALLOCATABLE
      integer k_singles,l_singles
      integer k_doubles,l_doubles
#endif // USE_F90_ALLOCATABLE
      integer size,i
      integer nxtask
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
      integer xmem
      integer dp4,dp5,ii,jj,istart,istop,jstart,jstop
      integer mdp4,mdp5
      integer maxp4,maxp5
      integer slice_dp4,slice_dp5,qp4,qp5
c      double precision xlocal
      integer xlocal
c ---------------------------------------
      double precision energy(2)
      double precision energy1,energy2
      double precision factor,denom
      double precision denom_p4,denom_p5,denom_p6
      double precision denom_h1,denom_h2,denom_h3
      external nxtask
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
      count = 0
      next = nxtask(nprocs,1)
      energy(1)=0.0d0
      energy(2)=0.0d0
      energy1 = 0.0d0
      energy2 = 0.0d0
      do t_p4b = noab+1,noab+nvab
       range_p4 = int_mb(k_range+t_p4b-1)
       offset_p4 = k_evl_sorted+int_mb(k_offset+t_p4b-1)-1
       do t_p5b = t_p4b,noab+nvab
        range_p5 = int_mb(k_range+t_p5b-1)
        offset_p5 = k_evl_sorted+int_mb(k_offset+t_p5b-1)-1
        do t_p6b = t_p5b,noab+nvab
         range_p6 = int_mb(k_range+t_p6b-1)
         offset_p6 = k_evl_sorted+int_mb(k_offset+t_p6b-1)-1
         do t_h1b = 1,noab
          range_h1 = int_mb(k_range+t_h1b-1)
          offset_h1 = k_evl_sorted+int_mb(k_offset+t_h1b-1)-1
          do t_h2b = t_h1b,noab
           range_h2 = int_mb(k_range+t_h2b-1)
           offset_h2 = k_evl_sorted+int_mb(k_offset+t_h2b-1)-1
           do t_h3b = t_h2b,noab
            range_h3 = int_mb(k_range+t_h3b-1)
            offset_h3 = k_evl_sorted+int_mb(k_offset+t_h3b-1)-1
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
c -- 6DTS --
            dp4=1
            dp5=1
  300       continue
cx            if(dp4.eq.int_mb(k_range+t_p4b-1).and.
cx     1         dp5.eq.int_mb(k_range+t_p5b-1))
cx     2      call errquit('6DTS error',3,MA_ERR)
            xlocal=((int_mb(k_range+t_p4b-1)/dp4)
     1           *  (int_mb(k_range+t_p5b-1)/dp5))
     2           * int_mb(k_range+t_p6b-1)
     3           * int_mb(k_range+t_h1b-1)
     4           * int_mb(k_range+t_h2b-1)
     5           * int_mb(k_range+t_h3b-1)
c
            slice_dp4=range_p4/dp4
            qp4=range_p4/slice_dp4
            slice_dp5=range_p5/dp5
            qp5=range_p5/slice_dp5
c
            if((xlocal)/((1024*1024)/8).gt.(xmem)) then
             if(dp4.lt.int_mb(k_range+t_p4b-1)) then
              dp4=dp4+1
             else
              dp5=dp5+1
             end if
             go to 300
            end if
c *** debug ***
c       write(6,1021) dp4,dp5,slice_dp4,slice_dp5,qp4,qp5
c       write(6,1023) range_p4,range_p5
c       call util_flush(6)
c 1021  format('dp4 dp5 slice_dp4 slice_dp5 qp4 qp5 ',6i5)
c 1023  format('range_p4 range_p5 ',2i5)
c *************
c
            if(slice_dp4*qp4.ne.range_p4) then
              mdp4=qp4+1
            else
              mdp4=qp4
            end if
            if(slice_dp5*qp5.ne.range_p5) then
              mdp5=qp5+1
            else
              mdp5=qp5
            end if
c
c *** debug ***
c       write(6,1022) mdp4,mdp5
c       call util_flush(6)
c 1022  format('mdp4 mdp5 ',2i5)
c *************
c
cx           size = (int_mb(k_range+t_p4b-1)/dp4+
cx     6    int_mb(k_range+t_p4b-1)-(int_mb(k_range+t_p4b-1)/dp4)*dp4)
cx     6           *(int_mb(k_range+t_p5b-1)/dp5+
cx     6    int_mb(k_range+t_p5b-1)-(int_mb(k_range+t_p5b-1)/dp5)*dp5)
cx     2           * int_mb(k_range+t_p6b-1)
cx     3           * int_mb(k_range+t_h1b-1)
cx     4           * int_mb(k_range+t_h2b-1)
cx     5           * int_mb(k_range+t_h3b-1)
c
c *** debug ***
c            write(6,*)'dp4 dp5 ',dp4,dp5
c            call util_flush(6)
c *************
           size =  slice_dp4
     6           * slice_dp5
     2           * range_p6
     3           * range_h1
     4           * range_h2
     5           * range_h3

c
cx          maxp4=int_mb(k_range+t_p4b-1)/dp4+
cx     6    int_mb(k_range+t_p4b-1)-(int_mb(k_range+t_p4b-1)/dp4)*dp4
cx          maxp5=int_mb(k_range+t_p5b-1)/dp5+
cx     6    int_mb(k_range+t_p5b-1)-(int_mb(k_range+t_p5b-1)/dp5)*dp5
          maxp4=slice_dp4
          maxp5=slice_dp5
c
c
c
c ----------
c ----------
#if USE_F90_ALLOCATABLE
            allocate(f_singles(1:size), stat=alloc_error)
            if (alloc_error.ne.0) then
              call errquit('ccsd_t_6dts: allocate failed',
     1                     alloc_error, MA_ERR)
            endif
            allocate(f_doubles(1:size), stat=alloc_error)
            if (alloc_error.ne.0) then
              call errquit('ccsd_t_6dts: allocate failed',
     1                     alloc_error, MA_ERR)
            endif
#else  // USE_F90_ALLOCATABLE
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) singles',l_singles,
     1        k_singles)) call errquit('ccsd_t: MA error',1,MA_ERR)
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) doubles',l_doubles,
     1        k_doubles)) call errquit('ccsd_t: MA error',2,MA_ERR)
#endif  // USE_F90_ALLOCATABLE
c  6DTS ----
           do ii=1,mdp4 !do ii ---
           do jj=1,mdp5 !do jj ---
            if(ii.ne.mdp4) then
             istart=(ii-1)*(slice_dp4)+1
             istop=ii*(slice_dp4)
            else
             istart=(ii-1)*(slice_dp4)+1
             istop=int_mb(k_range+t_p4b-1)
            end if
            if(jj.ne.mdp5) then
             jstart=(jj-1)*(slice_dp5)+1
             jstop=jj*(slice_dp5)
            else
             jstart=(jj-1)*(slice_dp5)+1
             jstop=int_mb(k_range+t_p5b-1)
            end if
c zeroing ---
#if USE_F90_ALLOCATABLE
        call dfill(size, 0.0d0, f_singles, 1)
        call dfill(size, 0.0d0, f_doubles, 1)
#else  // USE_F90_ALLOCATABLE
        call dfill(size, 0.0d0, dbl_mb(k_singles), 1)
        call dfill(size, 0.0d0, dbl_mb(k_doubles), 1)
#endif  // USE_F90_ALLOCATABLE
c -----------
cjeff            print*,'#DTS#',istart,istop,jstart,jstop,maxp4,maxp5
#if USE_F90_ALLOCATABLE
            call ccsd_t_singles_l_6dts(f_singles,
#else  // USE_F90_ALLOCATABLE
            call ccsd_t_singles_l_6dts(dbl_mb(k_singles),
#endif  // USE_F90_ALLOCATABLE
     1        k_t1_local,d_v2,k_t1_offset,
     1        k_v2_offset,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2,
     1        istart,istop,jstart,jstop,maxp4,maxp5)
#if USE_F90_ALLOCATABLE
            call ccsd_t_doubles_l_6dts(f_doubles,
#else  // USE_F90_ALLOCATABLE
            call ccsd_t_doubles_l_6dts(dbl_mb(k_doubles),
#endif  // USE_F90_ALLOCATABLE
     1        d_t2,d_v2,k_t2_offset,
     1        k_v2_offset,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2,
     1        istart,istop,jstart,jstop,maxp4,maxp5)
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
c
c factor = [ 1/36, 1/18, 1/12, 1/6, 1/4, 1/3, 1/2, 1, 2]
c
            i = 0
            do t_p4 = istart,istop
             denom_p4 = dbl_mb(offset_p4+t_p4)
             do t_p5 = jstart,jstop
              denom_p5 = dbl_mb(offset_p5+t_p5)
              do t_p6 = 1, range_p6
               denom_p6 = dbl_mb(offset_p6+t_p6)
               do t_h1 = 1, range_h1
                denom_h1 = dbl_mb(offset_h1+t_h1)
                do t_h2 = 1, range_h2
                 denom_h2 = dbl_mb(offset_h2+t_h2)
                 do t_h3 = 1, range_h3
                  denom_h3 = dbl_mb(offset_h3+t_h3)
                  denom = 1.0d0 / ( ( denom_h1 + denom_h2 + denom_h3 )
     1                            - ( denom_p4 + denom_p5 + denom_p6 ) )
#if USE_F90_ALLOCATABLE
                  energy1 = energy1 + factor*denom
     1                    * f_doubles(i+1) * f_doubles(i+1)
                  energy2 = energy2 + factor*denom*f_doubles(i+1)
     1                    * (f_doubles(i+1)+f_singles(i+1))
#else // USE_F90_ALLOCATABLE
                  energy1 = energy1 + factor*denom
     1                    * dbl_mb(k_doubles+i)*dbl_mb(k_doubles+i)
                  energy2 = energy2 + factor*denom*dbl_mb(k_doubles+i)
     1                    * (dbl_mb(k_doubles+i)+dbl_mb(k_singles+i))
#endif // USE_F90_ALLOCATABLE
                  i = i + 1
                 enddo
                enddo
               enddo
              enddo
             enddo
            enddo
c
            enddo !do ii ---
            enddo !do jj ---
c
#if USE_F90_ALLOCATABLE
            deallocate(f_singles, stat=alloc_error)
            if (alloc_error.ne.0) then
              call errquit('ccst_t_6dts: deallocate failed',
     1                     alloc_error, MA_ERR)
            endif
            deallocate(f_doubles, stat=alloc_error)
            if (alloc_error.ne.0)  then
              call errquit('ccsd_t_6dts: deallocate failed',
     1                     alloc_error, MA_ERR)
            endif
#else // USE_F90_ALLOCATABLE
            if (.not.MA_POP_STACK(l_doubles))
     1        call errquit('ccsd_t',3,MA_ERR)
            if (.not.MA_POP_STACK(l_singles))
     1        call errquit('ccsd_t',4,MA_ERR)
#endif // USE_F90_ALLOCATABLE
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
      energy(1) = energy1
      energy(2) = energy2
      call ga_dgop(mt_dbl,energy,2,'+')
      energy1 = energy(1)
      energy2 = energy(2)
c - T1/X1 LOCALIZATION ------
         if(.not.MA_POP_STACK(l_t1_local))
     &      call errquit('ccsd_t: l_t1_local',4,MA_ERR)
c ---------------------------
      return
      end

