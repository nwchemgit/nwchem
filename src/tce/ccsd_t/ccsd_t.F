      SUBROUTINE ccsd_t(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                  d_v2,k_v2_offset,energy1,energy2)
C
C     $Id: ccsd_t.F,v 1.11 2008-10-08 16:41:07 jhammond Exp $
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
      integer nxtask
      integer next
      integer nprocs
      integer count
      integer offset_p4,offset_p5,offset_p6
      integer offset_h1,offset_h2,offset_h3
      integer range_p4,range_p5,range_p6
      integer range_h1,range_h2,range_h3
      double precision energy(2)
      double precision energy1,energy2
      double precision factor,denom
      double precision denom_p4,denom_p5,denom_p6
      double precision denom_h1,denom_h2,denom_h3
      external nxtask
C
      nprocs = GA_NNODES()
      count = 0
      next = nxtask(nprocs,1)
      energy(1)=0.0d0
      energy(2)=0.0d0
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
            size = range_p4 * range_p5 * range_p6
     3           * range_h1 * range_h2 * range_h3
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) singles',l_singles,
     1        k_singles)) call errquit('ccsd_t: MA error',1,MA_ERR)
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) doubles',l_doubles,
     1        k_doubles)) call errquit('ccsd_t: MA error',2,MA_ERR)
            do i = 1, size
              dbl_mb(k_singles+i-1) = 0.0d0
            enddo
            do i = 1, size
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
            do t_p4 = 1, range_p4
             denom_p4 = dbl_mb(offset_p4+t_p4)
             do t_p5 = 1, range_p5
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
                  energy1 = energy1 + factor*denom
     1                    * dbl_mb(k_doubles+i)*dbl_mb(k_doubles+i)
                  energy2 = energy2 + factor*denom*dbl_mb(k_doubles+i)
     1                    * (dbl_mb(k_doubles+i)+dbl_mb(k_singles+i))
                  i = i + 1
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
      energy(1) = energy1
      energy(2) = energy2
      call ga_dgop(mt_dbl,energy,2,'+')
      energy1 = energy(1)
      energy2 = energy(2)
      return
      end




#ifndef DEBUG_CCSD_T
c --------------------------------------------------------------
c
c                 Only does CCSD(T) not CCSD[T]
c
c --------------------------------------------------------------
      SUBROUTINE ccsd_t_vector(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                         d_v2,k_v2_offset,energy)
C
C     $Id: ccsd_t.F,v 1.11 2008-10-08 16:41:07 jhammond Exp $
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
      integer nxtask
      integer next
      integer nprocs
      integer count
      integer offset_p4,offset_p5,offset_p6
      integer offset_h1,offset_h2,offset_h3
      integer range_p4,range_p5,range_p6
      integer range_h1,range_h2,range_h3
      integer range_h3_b,range_h3_r
      double precision energy
      double precision factor,denom,denom_p,denom_h12,denom_ph12
      double precision denom_p4,denom_p5,denom_p6
      double precision denom_h1,denom_h2,denom_h3
      external nxtask
C
      nprocs = GA_NNODES()
      count = 0
      next = nxtask(nprocs,1)
      energy = 0.0d0
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
            range_h3_r = mod(range_h3,4)              ! remainder
            range_h3_b = ( range_h3 - range_h3_r )/4  ! block part
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
            size = range_p4 * range_p5 * range_p6
     3           * range_h1 * range_h2 * range_h3
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) singles',l_singles,
     1        k_singles)) call errquit('ccsd_t: MA error',1,MA_ERR)
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) doubles',l_doubles,
     1        k_doubles)) call errquit('ccsd_t: MA error',2,MA_ERR)
#ifdef IFCV8
!DEC$ loop count min(4096), avg(64000000)
!DEC$ unroll(4)
!DEC$ vector always
#endif
            do i = 1, size
              dbl_mb(k_singles+i-1) = 0.0d0
            enddo
#ifdef IFCV8
!DEC$ loop count min(4096), avg(64000000)
!DEC$ unroll(8)
!DEC$ vector always
#endif
            do i = 1, size
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
#ifdef IFCV8
!DEC$ ivdep
!DEC$ loop count min(4), max(32), avg(16)
#endif
            do t_p4 = 1, range_p4
             denom_p4 = dbl_mb(offset_p4+t_p4)
#ifdef IFCV8
!DEC$ loop count min(4), max(32), avg(16)
#endif
             do t_p5 = 1, range_p5
              denom_p5 = dbl_mb(offset_p5+t_p5)
#ifdef IFCV8
!DEC$ loop count min(4), max(32), avg(16)
#endif
              do t_p6 = 1, range_p6
               denom_p6 = dbl_mb(offset_p6+t_p6)
               denom_p = -1.0d0 * (denom_p6 + denom_p5 + denom_p4)
#ifdef IFCV8
!DEC$ loop count min(4), max(32), avg(16)
#endif
               do t_h1 = 1, range_h1
                denom_h1 = dbl_mb(offset_h1+t_h1)
#ifdef IFCV8
!DEC$ loop count min(4), max(32), avg(16)
#endif
                do t_h2 = 1, range_h2
                 denom_h2 = dbl_mb(offset_h2+t_h2)
                 denom_h12 = denom_h1 + denom_h2
                 denom_ph12 = denom_p + denom_h1 + denom_h2
#ifdef IFCV8
!DEC$ loop count min(4), max(32), avg(16)
!DEC$ unroll(4)
!DEC$ vector always
#endif
                 do t_h3 = 1, range_h3_b,4
                  denom_h3 = dbl_mb(offset_h3+t_h3)
                  denom = 1.0d0 / ( denom_ph12 + denom_h3 )
                  energy = energy + factor*denom*dbl_mb(k_doubles+i)
     2                   * (dbl_mb(k_doubles+i)+dbl_mb(k_singles+i))
                  i = i + 1
                  denom_h3 = dbl_mb(offset_h3+t_h3+1)
                  denom = 1.0d0 / ( denom_ph12 + denom_h3 )
                  energy = energy + factor*denom*dbl_mb(k_doubles+i)
     2                   * (dbl_mb(k_doubles+i)+dbl_mb(k_singles+i))
                  i = i + 1
                  denom_h3 = dbl_mb(offset_h3+t_h3+2)
                  denom = 1.0d0 / ( denom_ph12 + denom_h3 )
                  energy = energy + factor*denom*dbl_mb(k_doubles+i)
     2                   * (dbl_mb(k_doubles+i)+dbl_mb(k_singles+i))
                  i = i + 1
                  denom_h3 = dbl_mb(offset_h3+t_h3+3)
                  denom = 1.0d0 / ( denom_ph12 + denom_h3 )
                  energy = energy + factor*denom*dbl_mb(k_doubles+i)
     2                   * (dbl_mb(k_doubles+i)+dbl_mb(k_singles+i))
                  i = i + 1
                 enddo ! t_h3
#ifdef IFCV8
!DEC$ loop count min(0), max(3), avg(2)
#endif
                 do t_h3 = 1, range_h3_r
                  denom_h3 = dbl_mb(offset_h3+t_h3)
                  denom = 1.0d0 / ( denom_p + denom_h12 + denom_h3 )
                  energy = energy + factor*denom*dbl_mb(k_doubles+i)
     2                   * (dbl_mb(k_doubles+i)+dbl_mb(k_singles+i))
                  i = i + 1
                 enddo ! t_h3
                enddo ! t_h2
               enddo ! t_h1
              enddo ! t_p6
             enddo ! t_p5
            enddo ! t_p4
            if (.not.MA_POP_STACK(l_doubles)) 
     1        call errquit('ccsd_t',3,MA_ERR)
            if (.not.MA_POP_STACK(l_singles)) 
     1        call errquit('ccsd_t',4,MA_ERR)
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
      call ga_dgop(mt_dbl,energy,1,'+')
      return
      end

#endif

#ifdef DEBUG_CCSD_T

c --------------------------------------------------------------
c
c                 Only does CCSD(T) not CCSD[T]
c
c --------------------------------------------------------------
      SUBROUTINE ccsd_t_vector(d_t1,k_t1_offset,d_t2,k_t2_offset,
     1                         d_v2,k_v2_offset,energy)
C
C     $Id: ccsd_t.F,v 1.11 2008-10-08 16:41:07 jhammond Exp $
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
      integer nxtask
      integer next
      integer nprocs
      integer count
      integer offset_p4,offset_p5,offset_p6
      integer offset_h1,offset_h2,offset_h3
      integer range_p4,range_p5,range_p6
      integer range_h1,range_h2,range_h3
      double precision energy,energy0
      double precision factor,denom,denom_p,denom_h
      double precision denom_p4,denom_p5,denom_p6
      double precision denom_h1,denom_h2,denom_h3
      external nxtask
C
      nprocs = GA_NNODES()
      count = 0
      next = nxtask(nprocs,1)
      energy = 0.0d0
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
            size = range_p4 * range_p5 * range_p6
     3           * range_h1 * range_h2 * range_h3
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) singles',l_singles,
     1        k_singles)) call errquit('ccsd_t: MA error',1,MA_ERR)
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) doubles',l_doubles,
     1        k_doubles)) call errquit('ccsd_t: MA error',2,MA_ERR)
            do i = 1, size
              dbl_mb(k_singles+i-1) = 0.0d0
            enddo
            do i = 1, size
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
#ifdef FORTRAN_VERSION
            print*,'======================================'
            print*,'         Fortran Version              '
            print*,'======================================'
c
c            print*,'i,dbl_mb(k_singles+i),dbl_mb(k_doubles+i)'
c            do i = 1,size
c              if ( ( dabs(dbl_mb(k_singles+i)).gt.(1.0d-8) ).or.
c     1             ( dabs(dbl_mb(k_doubles+i)).gt.(1.0d-8) )
c     2           ) then
c                write(6,350) i,dbl_mb(k_singles+i),dbl_mb(k_doubles+i)
c              endif
c            enddo
c
            energy = 0.0d0
c
            i = 0
            do t_p4 = 1, range_p4
             denom_p4 = dbl_mb(offset_p4+t_p4)
             do t_p5 = 1, range_p5
              denom_p5 = dbl_mb(offset_p5+t_p5)
              do t_p6 = 1, range_p6
               denom_p6 = dbl_mb(offset_p6+t_p6)
c
               denom_p = denom_p6 + denom_p5 + denom_p4
c
               do t_h1 = 1, range_h1
                denom_h1 = dbl_mb(offset_h1+t_h1)
                do t_h2 = 1, range_h2
                 denom_h2 = dbl_mb(offset_h2+t_h2)
                 do t_h3 = 1, range_h3
                  denom_h3 = dbl_mb(offset_h3+t_h3)
c
                  denom_h = denom_h1 + denom_h2 + denom_h3
c
                  denom = 1.0d0 / ( denom_h - denom_p )
c
                  energy0 = factor * denom * dbl_mb(k_doubles+i)
     2                    * (dbl_mb(k_doubles+i) + dbl_mb(k_singles+i))
                  energy = energy + energy0
c
c                  if (dabs(energy0).gt.(1.0d-8)) then
c
c                  write(6,600) t_p4,t_p5,t_p6,t_h1,t_h2,t_h3,
c     1              denom_p4,denom_p5,denom_p6,
c     2              denom_h1,denom_h2,denom_h3,
c     3              factor,denom_p,denom_h,denom,
c     4              dbl_mb(k_singles+i),dbl_mb(k_doubles+i),
c     5              energy0
c
c                  endif
c
                  i = i + 1
c
                 enddo ! t_h3
                enddo ! t_h2
               enddo ! t_h1
              enddo ! t_p6
             enddo ! t_p5
            enddo ! t_p4
c
c            write(6,325) 'Energy = ',energy,' after Fortran'
c
#endif
c            print*,'======================================'
c            print*,'         C + SSE2 Version             '
c            print*,'======================================'
c
c            energy = 999999999.0d0
c
c            write(6,325) 'Energy = ',energy,' before C+SSE2'
c
c        write(6,100) 'range_p4,range_p5,range_p6 = ',
c     1                range_p4,range_p5,range_p6
c        write(6,100) 'range_h1,range_h2,range_h3 = ',
c     1                range_h1,range_h2,range_h3
c        write(6,200) 'dbl_mb(o_p4+1),dbl_mb(o_p5+1),dbl_mb(o_p6+1) = ',
c     1       dbl_mb(offset_p4+1),dbl_mb(offset_p5+1),dbl_mb(offset_p6+1)
c        write(6,200) 'dbl_mb(o_h1+1),dbl_mb(o_h2+1),dbl_mb(o_h3+1) = ',
c     1       dbl_mb(offset_h1+1),dbl_mb(offset_h2+1),dbl_mb(offset_h3+1)
c        write(6,300) 'dbl_mb(k_singles),dbl_mb(k_doubles) = ',
c     1                dbl_mb(k_singles),dbl_mb(k_doubles)
            call ccsd_t_sse2(range_p4,range_p5,range_p6,
     1                       range_h1,range_h2,range_h3,
     2      dbl_mb(offset_p4),dbl_mb(offset_p5),dbl_mb(offset_p6),
     3      dbl_mb(offset_h1),dbl_mb(offset_h2),dbl_mb(offset_h3),
     4                       dbl_mb(k_singles),dbl_mb(k_doubles),
     5                       factor,energy0)
c
            energy = energy + energy0
c
c            write(6,325) 'Energy = ',energy,' after C+SSE2'
c
c            print*,'======================================'
            if (.not.MA_POP_STACK(l_doubles)) 
     1        call errquit('ccsd_t',3,MA_ERR)
            if (.not.MA_POP_STACK(l_singles)) 
     1        call errquit('ccsd_t',4,MA_ERR)
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
      call ga_dgop(mt_dbl,energy,1,'+')
      return
100   format(1x,a60,3i18)
200   format(1x,a60,3f18.8)
300   format(1x,a60,2f18.8)
325   format(1x,a9,f18.8,a20)
350   format(1x,i8,2f18.8)
400   format(1x,'t_*',6i8)
500   format(1x,a12,f18.8)
600   format(6i4,13f12.8)
      end

#endif
