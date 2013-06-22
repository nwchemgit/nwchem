ckbn ccsd_t.F modified to tce_mrcc_ccsdpt_uncoup_pt3.F 

       SUBROUTINE tce_mrcc_ccsdpt(d_t1,k_t1_offset,
     +                           d_t2,k_t2_offset,
     +                           d_v2,k_v2_offset,
     +                           d_f1,
     +                           k_f1_offset,
ckbn +                           energy1,energy2,energy3,
     +                           energy3,
     +                           iref)
c     1                           iref,nref)

      IMPLICIT NONE
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "sf.fh"
#include "inp.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_hetio.fh"
#include "tce_diis.fh"
#include "tce_prop.fh"
#include "tce_restart.fh"
#include "tce_mrcc.fh"

c      integer nref
      integer iref
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
ckbn-2
c      integer k_q3fnt2,l_q3fnt2
      integer size,i,j
      integer nxtask
      integer next
      integer nprocs
      integer count
      integer offset_p4,offset_p5,offset_p6
      integer offset_h1,offset_h2,offset_h3
      integer range_p4,range_p5,range_p6
      integer range_h1,range_h2,range_h3
      double precision energy(3)
ckbn      double precision energy1,energy2
ckbn-2
      double precision energy3
      double precision factor,denom
      double precision denom_p4,denom_p5,denom_p6
      double precision denom_h1,denom_h2,denom_h3
      external nxtask
ckbn-3
      logical nodezero
      integer d_f1,k_f1_offset
      double precision dnorm1,dnorm2,dmaxt1,dmaxt2
      double precision dmint1,dmint2
      integer ioff
      integer orbspin(6),orbindex(6),off(6),blck(6)
      INTEGER NXTASKsub
      EXTERNAL NXTASKsub

ckbn-2
      nodezero=(ga_nodeid().eq.0)

ckbn-2
c      if(nodezero) write(*,*)"Entering tce_mrcc_ccsdpt"

      nprocs = GA_NNODES()
      count = 0
      next = nxtask(nprocs,1)



      energy(1)=0.0d0
      energy(2)=0.0d0
      energy(3)=0.0d0
c      energy1  =0.0d0
c      energy2  =0.0d0
      energy3  =0.0d0
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

      IF (.not.((log_mb(k_active_tmpm(iref)+t_p4b-1).eqv..true.).and.
     1  (log_mb(k_active_tmpm(iref)+t_p5b-1).eqv..true.).and.
     1  (log_mb(k_active_tmpm(iref)+t_p6b-1).eqv..true.).and.
     1  (log_mb(k_active_tmpm(iref)+t_h1b-1).eqv..true.).and.
     1  (log_mb(k_active_tmpm(iref)+t_h2b-1).eqv..true.).and.
     1  (log_mb(k_active_tmpm(iref)+t_h3b-1).eqv..true.))) then

            size = range_p4 * range_p5 * range_p6
     3           * range_h1 * range_h2 * range_h3
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) singles',l_singles,
     1        k_singles)) call errquit('ccsd_t: MA error',1,MA_ERR)
            if (.not.MA_PUSH_GET(mt_dbl,size,'(T) doubles',l_doubles,
     1        k_doubles)) call errquit('ccsd_t: MA error',2,MA_ERR)
            do i = 0, size-1
              dbl_mb(k_singles+i) = 0.0d0
              dbl_mb(k_doubles+i) = 0.0d0
            enddo

           call tce_mrcc_q3vnt2(dbl_mb(k_doubles),d_t2,d_v2,k_t2_offset,
     1        k_v2_offset,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
           call tce_mrcc_q3vnt1(dbl_mb(k_singles),d_t1,d_v2,k_t1_offset,
     1        k_v2_offset,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
           call tce_mrcc_q3fnt2(dbl_mb(k_singles),d_f1,d_t2,k_f1_offset,
     +           k_t2_offset,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)

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


                  if(.not.lusesamefock_nonit) then
                  denom = 1.0d0 / ((denom_h1 + denom_h2 + denom_h3)
     +                           - (denom_p4 + denom_p5 + denom_p6))


                  else
                  orbspin(1) = int_mb(k_spin+t_p4b-1)-1
                  orbspin(2) = int_mb(k_spin+t_p5b-1)-1
                  orbspin(3) = int_mb(k_spin+t_p6b-1)-1
                  orbspin(4) = int_mb(k_spin+t_h1b-1)-1
                  orbspin(5) = int_mb(k_spin+t_h2b-1)-1
                  orbspin(6) = int_mb(k_spin+t_h3b-1)-1

                  orbindex(1) = (1-orbspin(1)+int_mb(k_mo_indexm(iref)+
     1           int_mb(k_offset+t_p4b-1)+t_p4-1))/2
                  orbindex(2) = (1-orbspin(2)+int_mb(k_mo_indexm(iref)+
     1           int_mb(k_offset+t_p5b-1)+t_p5-1))/2
                  orbindex(3) = (1-orbspin(3)+int_mb(k_mo_indexm(iref)+
     1           int_mb(k_offset+t_p6b-1)+t_p6-1))/2
                  orbindex(4) = (1-orbspin(4)+int_mb(k_mo_indexm(iref)+
     1           int_mb(k_offset+t_h1b-1)+t_h1-1))/2
                  orbindex(5) = (1-orbspin(5)+int_mb(k_mo_indexm(iref)+
     1           int_mb(k_offset+t_h2b-1)+t_h2-1))/2
                  orbindex(6) = (1-orbspin(6)+int_mb(k_mo_indexm(iref)+
     1           int_mb(k_offset+t_h3b-1)+t_h3-1))/2

                 do j=1,6
                  blck(j) = orbinblck(orbindex(j),orbspin(j)+1,1)
                  off(j) = offsetinblck(orbindex(j),orbspin(j)+1,1)
                 enddo

                denom = 1.0d0/(
     +                  -dbl_mb(k_evl_sortedm(1)+int_mb(k_offsetm(1)
     +                         +blck(1)-1)+off(1))
     +                  -dbl_mb(k_evl_sortedm(1)+int_mb(k_offsetm(1)
     +                         +blck(2)-1)+off(2))
     +                  -dbl_mb(k_evl_sortedm(1)+int_mb(k_offsetm(1)
     +                         +blck(3)-1)+off(3))
     +                  +dbl_mb(k_evl_sortedm(1)+int_mb(k_offsetm(1)
     +                         +blck(4)-1)+off(4))
     +                  +dbl_mb(k_evl_sortedm(1)+int_mb(k_offsetm(1)
     +                         +blck(5)-1)+off(5))
     +                  +dbl_mb(k_evl_sortedm(1)+int_mb(k_offsetm(1)
     +                         +blck(6)-1)+off(6)))
                endif ! lusesamefock_nonit

                if((abs(denom)).gt.1.0d2) then
c                 if(nodezero) 
                 write(LuOut,*)
     +           'Warning:Denominator is very low. 1/D=',abs(denom)
c                 if(nodezero) 
                 call util_flush(LuOut)
                endif
               

c                  energy1 = energy1 + factor*denom
c     1                    * dbl_mb(k_doubles+i)*dbl_mb(k_doubles+i)
c                  energy2 = energy2 + factor*denom*dbl_mb(k_doubles+i)
c     1                    * (dbl_mb(k_doubles+i)+dbl_mb(k_singles+i))
                  energy3 = energy3 + factor*denom*dbl_mb(k_doubles+i)
     +                    * (dbl_mb(k_doubles+i)+dbl_mb(k_singles+i))


                  i = i + 1
                 enddo
                enddo
               enddo
              enddo
             enddo
            enddo
c            if (.not.MA_POP_STACK(l_q3fnt2)) 
c     1        call errquit('ccsd_t',2,MA_ERR)
            if (.not.MA_POP_STACK(l_doubles)) 
     1        call errquit('ccsd_t',3,MA_ERR)
            if (.not.MA_POP_STACK(l_singles)) 
     1        call errquit('ccsd_t',4,MA_ERR)
            ENDIF
            endif
            endif
            endif
             if(lusesub) then
              next = NXTASKsub(nprocs,1,mypgid)
             else
              next = nxtask(nprocs,1)
             endif
            endif
            count = count + 1
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      if(lusesub) then
      next = NXTASKsub(-nprocs,1,mypgid)
      else
      next = nxtask(-nprocs,1)
      call ga_sync()
      endif
c      energy(1) = energy1
c      energy(2) = energy2
      energy(3) = energy3
      if(lusesub) then
      call ga_pgroup_sync(mypgid)
c      call ga_pgroup_dgop(p_handle, type, buf, lenbuf, op)
      call ga_pgroup_dgop(mypgid, mt_dbl,energy,3,'+')
      else 
      call ga_dgop(mt_dbl,energy,3,'+')
      endif
c      energy1 = energy(1)
c      energy2 = energy(2)
      energy3 = energy(3)

c       if(lusesub)then
c       if(ga_pgroup_nodeid(mypgid).eq.0) then
cc      write(*,'(A,F17.10,A,F17.10,A,F17.10)')"energy1 ",energy1,
cc     +   " energy2 ",energy2," energy3 ",energy3
c      write(*,'(A,F17.10,A,F17.10,A,F17.10)')"energy3 ",energy3
c       endif
c       else
c      if(nodezero) then
cc      write(*,'(A,F17.10,A,F17.10,A,F17.10)')"energy1 ",energy1,
cc     +   " energy2 ",energy2," energy3 ",energy3
c      write(*,'(A,F17.10,A,F17.10,A,F17.10)')"energy3 ",energy3
c      endif
c
c      endif
c
      return
      end

c $Id$
