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
c - T1/X1 LOCALIZATION -------------------
      integer l_t1_local,k_t1_local
      integer size_t1
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
c     init sx for offload
      triplesx_alloced=.false.
      triplesx1_alloced=.false.
      v2sub_alloced=.false.
      t1sub_alloced=.false.
      t2sub_alloced=.false.
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
c      call util_align64(size)
      triplesx_mxlgth=size
      if (.not.MA_PUSH_GET(mt_dbl,size,'(T) singles',l_singles,
     1     k_singles)) call errquit('ccsd_t: MA error sgl',
     2     size,MA_ERR)
      if (.not.MA_PUSH_GET(mt_dbl,size,'(T) doubles',l_doubles,
     1     k_doubles)) call errquit('ccsd_t: MA error dbl',
     2     size,MA_ERR)
#ifdef USE_OFFLOAD
      call offl_alloc(dbl_mb(k_singles),size)
      triplesx1_alloced=.true.
      call offl_alloc(dbl_mb(k_doubles),size)
      triplesx_alloced=.true.
#endif      
c     stagger start of loop
      call util_mpinap(100)
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
            size = range_p4 * range_p5 * range_p6
     3           * range_h1 * range_h2 * range_h3
c zeroing ---
        call dcopy(size, 0.0d0, 0, dbl_mb(k_singles), 1)
        call offl_zerofill(dbl_mb(k_singles),size)
        call dcopy(size, 0.0d0, 0, dbl_mb(k_doubles), 1)
      call offl_zerofill(dbl_mb(k_doubles),size)
c -----------
#ifdef USE_OFFLOAD
#define USE_OMP_TRIPLES_DOT_PRODUCT 1
#warning offload code enabled
            call offl_ccsd_t_singles_l(dbl_mb(k_singles),
     1        k_t1_local,d_v2,k_t1_offset,
     1        k_v2_offset,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
#else
            call ccsd_t_singles_l(dbl_mb(k_singles),
     1        k_t1_local,d_v2,k_t1_offset,
     1        k_v2_offset,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
#endif
            call ccsd_t_doubles_l(dbl_mb(k_doubles),
     1        d_t2,d_v2,k_t2_offset,
     1        k_v2_offset,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
#ifdef USE_OMP_TRIPLES_DOT_PRODUCT
            call ccsd_t_esum(
     C           dbl_mb(k_singles),dbl_mb(k_doubles),
     R           restricted,
     &                      t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,
     &                      dbl_mb(offset_h1+1),dbl_mb(offset_h2+1),
     D           dbl_mb(offset_h3+1),
     &                      dbl_mb(offset_p4+1),dbl_mb(offset_p5+1),
     D           dbl_mb(offset_p6+1),
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
c
c factor = [ 1/36, 1/18, 1/12, 1/6, 1/4, 1/3, 1/2, 1, 2]
c
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
#endif
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
c cleanup to dealloc
c      triplesx_free=.true.
c      call ccsd_t_doubles_l(dbl_mb(k_doubles),
c     1     d_t2,d_v2,k_t2_offset,
c     1     k_v2_offset,1,1,1,1,1,1,2)
#ifdef USE_OFFLOAD
      call offl_free(dbl_mb(k_doubles),size)
      triplesx_alloced=.false.
      call offl_free(dbl_mb(k_singles),size)
      triplesx1_alloced=.false.
#endif

      if (.not.MA_POP_STACK(l_doubles)) 
     1     call errquit('ccsd_t doubles',3,MA_ERR)
      if (.not.MA_POP_STACK(l_singles)) 
     1     call errquit('ccsd_t singles',4,MA_ERR)
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
c     wrapper to ccsd_t_dot because of offload ugliness
      subroutine ccsd_t_esum(
     A     a_singles, a_doubles,
     C     restricted,
     &                      h1b,h2b,h3b,p4b,p5b,p6b,
     &                      o_h1,o_h2,o_h3,
     &                      o_p4,o_p5,o_p6,
     &                      r_h1,r_h2,r_h3,
     &                      r_p4,r_p5,r_p6,
     &                      energy1,energy2)
      implicit none
#ifdef OFFLOAD_CODE
#ifdef OPENMP_OFFLOAD
!$omp declare target (offl_ccsd_t_dot)
#else
cdir$ ATTRIBUTES OFFLOAD : mic :: offl_ccsd_t_dot
#endif
#endif
      integer h1b, h2b, h3b, p4b, p5b, p6b
      double precision o_h1(*),o_h2(*),o_h3(*)
      double precision o_p4(*),o_p5(*),o_p6(*)
      integer r_h1,r_h2,r_h3
      integer r_p4,r_p5,r_p6
      double precision a_singles(*)
      double precision a_doubles(*)
      logical restricted
      double precision energy1,energy2
#ifdef USE_OFFLOAD
      integer util_micdev,offload_master
      external util_micdev,offload_master
      integer mic_device
cdir$ ATTRIBUTES OFFLOAD : mic :: offl_ccsd_t_dot
      if(offload_master()) then
         mic_device=util_micdev()
CDIR$ OFFLOAD TARGET(mic:mic_device)
     N	IN(a_singles:length(0) REUSE)  
     N	IN(a_doubles:length(0) REUSE)  
     N  IN(restricted)
     I  IN(h1b,h2b,h3b,p4b,p5b,p6b)
     I  IN(r_h1,r_h2,r_h3,r_p4,r_p5,r_p6)
     N	IN(o_h1:length(r_h1))  
     N	IN(o_h2:length(r_h2))  
     N	IN(o_h3:length(r_h3))  
     N	IN(o_p4:length(r_p4))  
     N	IN(o_p5:length(r_p5))  
     N	IN(o_p6:length(r_p6))  
     I  INOUT(energy1,energy2)
         call offl_ccsd_t_dot(
     C        a_singles,a_doubles,
     R        restricted,
     &        h1b,h2b,h3b,p4b,p5b,p6b,
     &        o_h1,o_h2,o_h3,o_p4,o_p5,o_p6,
     &        r_h1,r_h2,r_h3,r_p4,r_p5,r_p6,
     &        energy1,energy2)
         else
#endif
            call ccsd_t_dot(
     C        a_singles,a_doubles,
     R        restricted,
     &        h1b,h2b,h3b,p4b,p5b,p6b,
     &        o_h1,o_h2,o_h3,o_p4,o_p5,o_p6,
     &        r_h1,r_h2,r_h3,r_p4,r_p5,r_p6,
     &        energy1,energy2)
#ifdef USE_OFFLOAD
         endif
#endif

      return
      end
      subroutine offl_free(ttt,l_ttt)
      implicit none
      double precision ttt(*)
      integer l_ttt
c
      integer util_micdev
      external util_micdev
      logical offload_master
      external offload_master
      integer mic_device
c
#ifdef USE_OFFLOAD
      if(offload_master()) then
        mic_device=util_micdev()
cc Free memory on MIC 
!DIR$ OFFLOAD_TRANSFER TARGET(mic:mic_device) 
     O     NOCOPY(ttt:length(0) 
     &   FREE)
      endif
#endif
      return
      end
      subroutine offl_alloc(ttt,l_ttt)
      implicit none
      double precision ttt(*)
      integer l_ttt
c
      integer util_micdev
      external util_micdev
      logical offload_master
      external offload_master
      integer mic_device
c
#ifdef USE_OFFLOAD
      if(offload_master()) then
        mic_device=util_micdev()
!DIR$ OFFLOAD_TRANSFER TARGET(mic:mic_device) 
     N  NOCOPY(ttt:length(l_ttt) 
     &  ALLOC) 
      endif
#endif
      return
      end
      subroutine offl_zerofill(ttt,l_ttt)
      implicit none
      double precision ttt(*)
      integer l_ttt
c
      integer util_micdev
      external util_micdev
      logical offload_master
      external offload_master
      integer mic_device
cdir$ ATTRIBUTES OFFLOAD : mic :: offl_zero
c
#ifdef USE_OFFLOAD
      if(offload_master()) then
        mic_device=util_micdev()
!DIR$ OFFLOAD TARGET(mic:mic_device) 
     N  IN(ttt:length(0) REUSE) 
     I  IN(l_ttt)
      call offl_zero(ttt,l_ttt)
      endif
#endif
      return
      end
