!
      subroutine tce_uss_offdiagonal_3(d_t2,d_v2,
     2                 k_t2_offset,k_v2_offset,d_c1,d_c2,iref,jref,
     1                 k_r1_offsetm,k_r2_offsetm)
 
! Routine for R selection
      implicit none
#include "tce.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "sym.fh"
#include "tce_mrcc.fh"
#include "global.fh"
#include "tce_main.fh"

      integer rtdb
      logical nodezero
      integer k_r2_offsetm(maxref)
      integer k_r1_offsetm(maxref)
c
      integer d_t2,d_v2
      integer k_t2_offset,k_v2_offset
c
      integer d_r2m(maxref),d_r1m(maxref)
      integer iref,jref
      integer i,j,p1,h2,k,p2,h3,h4,i1,p,k1,k2
      integer size,l,m,n,o
      integer l_r2,k_r2,l_r3,k_r3
      integer orbindex(8),aorbindex(8),orbindexnew(8)
      integer t, p1new,p2new, h1new, h2new,p3new,p4new,h3new,h4new
      integer p1new1,p2new1, h1new1, h2new1,p3new1,p4new1,h3new1,h4new1
      integer orbspin(8),aorbspin(8),actspin(2)
      integer ioccnew(maxorb,2),iocc0(maxorb,2)
      integer ihash
      integer p1off,p2off,h1off,h2off,p3off,p4off,h3off,h4off
      integer p1off1,p2off1,h1off1,h2off1,p3off1,p4off1,h3off1,h4off1
c
      integer t_p4b,t_p5b,t_p6b,t_h1b,t_h2b,t_h3b
c
c  intermediates
c
c     integer d_i1_1,d_i1_2,d_i1_3
c     integer k_i1_offset_1,k_i1_offset_2,
c    1k_i1_offset_3
c     integer l_i1_offset_1,l_i1_offset_2,
c    1l_i1_offset_3
c
c
      integer ispinfrom,ispinto
      integer totaloff,ilength
      integer l_tmp3,k_tmp3
      integer l_tmp1,k_tmp1
      integer l_d_c1,k_d_c1
      integer l_d_c2,k_d_c2
      integer sizenew,ihashold,ioff,ioff1,ioff2,sizenew1
      integer noabn,nvabn,counter,iactive
      integer d_c2,d_c1
      logical isfound
      logical ap1,ah2,ap2,ah1,ah3,ap3
      logical activep,activeh
      EXTERNAL NXTASKsub
      EXTERNAL nxtask
      INTEGER nxtask
      INTEGER nxt
      INTEGER nprocs
      INTEGER count,next
      integer iexfrom,iexto,iexspin
      integer noper,erank,optyp,eoper,signfact,noper2,eoper2,espin
      integer wtp, wth
      INTEGER logtoint
      INTEGER holep, particlep
      double precision dsmult,p3f,p2f,fact
      dimension eoper(4*maxexcit),optyp(4*maxexcit),eoper2(4*maxexcit),
     1 espin(4*maxexcit)
      dimension iexfrom(8),iexto(8),iexspin(8)



c     if(lusesub) then

c     call ga_pgroup_sync(mypgid)
c     nprocs = GA_pgroup_NNODES(mypgid)
c     count = 0
c     next = NXTASKsub(nprocs, 1,mypgid)

c     else
ccccccccccccccc below lines are commented for VT2 approximation ccccc
c
c      call ccsd_r3(dbl_mb(k_r3),d_f1,d_i1_1,
c    1d_i1_2,d_i1_3,
c    1d_t1,d_t2,d_v2,k_f1_offset,
c    1k_i1_offset_1,k_i1_offset_2,
c    1k_i1_offset_3,k_t1_offset,k_t2_offset,
c    1k_v2_offset,l_i1_offset_1,l_i1_offset_2,
c    1l_i1_offset_3,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,1)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call ga_sync()
      nprocs = GA_NNODES()
      count = -1
      next = NXTASK(nprocs,1)

c     endif


      nodezero = (ga_nodeid().eq.0)
!          
c       if(nodezero)write(*,*)'out of 3'

         noa = nblcks(1,iref)
         nob = nblcks(2,iref)
         nva = nblcks(3,iref)
         nvb = nblcks(4,iref)

         noab = noa+nob
         nvab = nva+nvb

         noabn = nblcks(1,jref)+nblcks(2,jref)
         nvabn = nblcks(3,jref)+nblcks(4,jref)
c-------
      DO t_p4b = noab+1,noab+nvab
      DO t_p5b = t_p4b,noab+nvab
      DO t_p6b = t_p5b,noab+nvab
      DO t_h1b = 1,noab
      DO t_h2b = t_h1b,noab
      DO t_h3b = t_h2b,noab
        if(next.eq.count) then
      IF ((.not.restricted).or.(int_mb(k_spinm(iref)+t_p4b-1)+
     1 int_mb(k_spinm(iref)+t_p5b-1
     1)+int_mb(k_spinm(iref)+t_p6b-1)+int_mb(k_spinm(iref)+t_h1b-1)+
     1int_mb(k_spinm(iref)+t_h2b-1)+int_mb(k_spinm(iref)+t_h3b-1).ne.
     112)) THEN
      IF (int_mb(k_spinm(iref)+t_p4b-1)+int_mb(k_spinm(iref)+t_p5b-1)+
     1 int_mb(k_spinm(iref)+t_p6b-1)
     & .eq. int_mb(k_spinm(iref)+t_h1b-1)+int_mb(k_spinm(iref)+t_h2b-1)+
     & int_mb(k_spinm(iref)+t_h3b-1)) THEN
 
      IF (ieor(int_mb(k_symm(iref)+t_p4b-1),ieor(int_mb(k_symm(iref)+
     1t_p5b-1),ieor(int_mb(k_symm(iref)+t_p6b-1),
     1 ieor(int_mb(k_symm(iref)+t_h1b-1),
     & ieor(int_mb(k_symm(iref)+t_h2b-1),
     1 int_mb(k_symm(iref)+t_h3b-1)))))).eq. irrep_t) THEN
c
      particlep= logtoint(log_mb(k_active+t_p4b-1))+
     1  logtoint(log_mb(k_active+t_p5b-1))+
     1  logtoint(log_mb(k_active+t_p6b-1))
       holep=logtoint(log_mb(k_active+t_h1b-1))+
     1  logtoint(log_mb(k_active+t_h2b-1))+
     1  logtoint(log_mb(k_active+t_h3b-1))

      if((holep.gt.0).and.(particlep.gt.0))
     1 then

      size = int_mb(k_rangem(iref)+t_p4b-1) * 
     1 int_mb(k_rangem(iref)+t_p5b-1)
     1 * int_mb(k_rangem(iref)+t_p6b-1)*int_mb(k_rangem(iref)+t_h1b-1)*
     1int_mb(k_rangem(iref)+t_h2b-1)*int_mb(k_rangem(iref)+t_h3b-1)
 
       l_r3=0
c
        if (.not.ma_push_get(mt_dbl,size,'c2',l_r3,k_r3))
     1   call errquit('tce_uss: MA problem',10,MA_ERR)
 
        CALL DFILL(size,0.0d0,dbl_mb(k_r3),1)
 
       call ccsd_r3(dbl_mb(k_r3),d_t2,d_v2,
     1k_t2_offset,k_v2_offset,
     1t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,2)
c

        counter =0 
        ihash =-1
        ihashold=-1
c-------------------------------
        do i=1,int_mb(k_rangem(iref)+t_p4b-1)
        do j=1,int_mb(k_rangem(iref)+t_p5b-1)
        do k=1,int_mb(k_rangem(iref)+t_p6b-1)
        do m=1,int_mb(k_rangem(iref)+t_h1b-1)
        do n=1,int_mb(k_rangem(iref)+t_h2b-1)
        do p=1,int_mb(k_rangem(iref)+t_h3b-1)
 
        counter = counter +1
 
 
      orbspin(1) = int_mb(k_spinm(iref)+t_p4b-1) -1
      orbspin(2) = int_mb(k_spinm(iref)+t_p5b-1) -1
      orbspin(3) = int_mb(k_spinm(iref)+t_p6b-1) -1
      orbspin(4) = int_mb(k_spinm(iref)+t_h1b-1) -1
      orbspin(5) = int_mb(k_spinm(iref)+t_h2b-1) -1
      orbspin(6) = int_mb(k_spinm(iref)+t_h3b-1) -1
 
 
       orbindex(1) = (1 - orbspin(1)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offsetm(iref)+t_p4b-1)+i-1))/2
 
       orbindex(2) = (1-orbspin(2)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offsetm(iref)+t_p5b-1)+j-1))/2
 
       orbindex(3) = (1 - orbspin(3)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offsetm(iref)+t_p6b-1)+k-1))/2
 
       orbindex(4) = (1-orbspin(4)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offsetm(iref)+t_h1b-1)+m-1))/2
 
       orbindex(5) = (1-orbspin(5)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offsetm(iref)+t_h2b-1)+n-1))/2
 
       orbindex(6) = (1-orbspin(6)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offsetm(iref)+t_h3b-1)+p-1))/2
 
c!sb sb sb sb sb sb sb sb
      iexto(1) = moindexes(orbindex(1),orbspin(1)+1,iref)
      iexto(2) = moindexes(orbindex(2),orbspin(2)+1,iref)
      iexto(3) = moindexes(orbindex(3),orbspin(3)+1,iref)
      iexfrom(1) = moindexes(orbindex(4),orbspin(4)+1,iref)
      iexfrom(2) = moindexes(orbindex(5),orbspin(5)+1,iref)
      iexfrom(3) = moindexes(orbindex(6),orbspin(6)+1,iref)
 
      iexspin(1) = orbspin(1)+1
      iexspin(2) = orbspin(2)+1
      iexspin(3) = orbspin(3)+1

       signfact=0
       dsmult=1.0d0
 
       call perfexcit(erank,eoper,iref,jref,
     &   3,iexfrom,iexto,iexspin,
     &   signfact,espin,wtp,wth)


      if(erank.eq.0) goto 111
      if(erank.eq.2) then
       signfact=mod(signfact,2)
       if(signfact.ne.0) then
        dsmult = -1.0d0
       endif
      orbindex(1)=moindexes(eoper(1),espin(1),jref)
      orbindex(2)=moindexes(eoper(2),espin(2),jref)

      p1new = orbinblck(orbindex(1),espin(1),jref)
      h1new = orbinblck(orbindex(2),espin(2),jref)

      p1off = offsetinblck(orbindex(1),espin(1),jref)
      h1off = offsetinblck(orbindex(2),espin(2),jref)

      ihashold = h1new-1+noabn*(p1new-noabn-1)

        ilength = int_mb(k_r1_offsetm(jref))
        totaloff = -1

        do k1 = 1, ilength
          if(int_mb(k_r1_offsetm(jref)+k1).eq.ihashold) then
             totaloff = 1
            goto 112
          endif
        enddo
 112     continue


        if((p1new.le.noabn).or.
     2     (h1new.gt.noabn))then
           totaloff=-1
        endif

      if(totaloff.ne.-1)then

        l_tmp3=0
      ioff = p1off*int_mb(k_rangem(jref)+h1new-1)+h1off

      sizenew = int_mb(k_rangem(jref)+h1new-1)*
     1           int_mb(k_rangem(jref)+p1new-1)
 
      if (.not.ma_push_get(mt_dbl,sizenew,'tmp1',l_tmp3,k_tmp3))
     1       call errquit('tce_uss: MA problem',3,MA_ERR)

 
        CALL DFILL(sizenew,0.0d0,dbl_mb(k_tmp3),1)

      dbl_mb(k_tmp3+ioff)=dbl_mb(k_r3+counter-1)*dsmult
     1/(p3f(t_p4b,t_p5b,t_p6b)*p3f(t_h1b,t_h2b,t_h3b))

      call add_hash_block(d_c1,dbl_mb(k_tmp3),sizenew,
     1   int_mb(k_r1_offsetm(jref)),h1new-1+noabn*(p1new-noabn-1))
 
       if (.not.ma_pop_stack(l_tmp3))
     1    call errquit('tce_uss: MA problem',4,MA_ERR)
         endif !totaloff
 
       endif!erank=2

      if(erank.eq.4) then

       signfact=mod(signfact,2)
       if(signfact.ne.0) then
        dsmult = -1.0d0
       endif
 
       orbindex(1)=moindexes(eoper(1),espin(1),jref)
       orbindex(2)=moindexes(eoper(2),espin(2),jref)
       orbindex(4)=moindexes(eoper(3),espin(3),jref)
       orbindex(3)=moindexes(eoper(4),espin(4),jref)
      p1new = orbinblck(orbindex(1),espin(1),jref)
      p2new = orbinblck(orbindex(2),espin(2),jref)
      p1off = offsetinblck(orbindex(1),espin(1),jref)
      p2off = offsetinblck(orbindex(2),espin(2),jref)

      h1new = orbinblck(orbindex(3),espin(4),jref)
      h2new = orbinblck(orbindex(4),espin(3),jref)
      h1off = offsetinblck(orbindex(3),espin(4),jref)
      h2off = offsetinblck(orbindex(4),espin(3),jref)
 
      if(p1new.gt.p2new) then
      t = p1new
      p1new = p2new
      p2new = t
      t = p1off
      p1off = p2off
      p2off = t
      dsmult=-1.0d0*dsmult
      end if
 
      if(h1new.gt.h2new) then
      t = h1new
      h1new = h2new
      h2new = t
      t = h1off
      h1off = h2off
      h2off = t
      dsmult=-1.0d0*dsmult
      end if

        if((p1new.le.noabn).or.
     1     (p2new.le.noabn).or.
     2     (h1new.gt.noabn).or.
     3     (h2new.gt.noabn)) goto 111

       if(((espin(1).eq.espin(2)).and.(espin(3).eq.espin(4)))
     1  .and.(((wtp.ge.2).and.(p1new.eq.p2new)).or.((wth.ge.2).and.
     1  (h1new.eq.h2new)))) then
        fact=2.d0
       else
       fact=1.d0
       endif

       ihashold = h2new-1+noabn*(h1new-1+noabn *(p2new-noabn-1+nvabn
     1 *(p1new-noabn-1)))
 
        ilength = int_mb(k_r2_offsetm(jref))
        totaloff = -1

        do k1 = 1, ilength
          if(int_mb(k_r2_offsetm(jref)+k1).eq.ihashold) then
             totaloff = 1
            goto 113
          endif
        enddo

 113     continue
 
      if(totaloff.ne.-1)then
         ioff = p1off*int_mb(k_rangem(jref)+h2new-1)*
     1 int_mb(k_rangem(jref)+h1new-1)*
     2 int_mb(k_rangem(jref)+p2new-1)+
     2 p2off*int_mb(k_rangem(jref)+h2new-1)*
     3 int_mb(k_rangem(jref)+h1new-1)+
     4 h1off*int_mb(k_rangem(jref)+h2new-1)+h2off

       sizenew = int_mb(k_rangem(jref)+p1new-1)*
     1 int_mb(k_rangem(jref)+p2new-1)*
     2 int_mb(k_rangem(jref)+h1new-1)*
     3 int_mb(k_rangem(jref)+h2new-1)
       l_tmp3=0
      if (.not.ma_push_get(mt_dbl,sizenew,'tmp1',l_tmp3,k_tmp3))
     1       call errquit('tce_uss: MA problem',3,MA_ERR)

        CALL DFILL(sizenew,0.0d0,dbl_mb(k_tmp3),1)

      dbl_mb(k_tmp3+ioff)=dbl_mb(k_r3+counter-1)*dsmult
     1*p2f(p1new,p2new)*p2f(h1new,h2new)/(p3f(t_p4b,t_p5b,t_p6b)*
     1p3f(t_h1b,t_h2b,t_h3b)*fact)

      call add_hash_block(d_c2,dbl_mb(k_tmp3),sizenew,
     1   int_mb(k_r2_offsetm(jref)),ihashold)
 
      if (.not.ma_pop_stack(l_tmp3))
     1    call errquit('tce_uss: MA problem',4,MA_ERR)
        endif !totaloff

        if((eoper(1).ne.eoper(2)).and.(eoper(3).ne.eoper(4)))then
        if((espin(1).eq.espin(2)).and.(espin(3).eq.espin(4))) then
        if(wtp.ge.2) then

       p1new1=p2new
       p2new1=p1new
       h1new1=h1new
       h2new1=h2new
       p1off1=p2off
       p2off1=p1off
       h1off1=h1off
       h2off1=h2off

      if(p1new1.gt.p2new1) then
      t = p1new1
      p1new1 = p2new1
      p2new1 = t
      t = p1off1
      p1off1 = p2off1
      p2off1 = t
      dsmult=-1.0d0*dsmult
      end if

      if(h1new1.gt.h2new1) then
      t = h1new1
      h1new1 = h2new1
      h2new1 = t
      t = h1off1
      h1off1 = h2off1
      h2off1 = t
      dsmult=-1.0d0*dsmult
      end if

        if((p1new1.le.noabn).or.
     1     (p2new1.le.noabn).or.
     2     (h1new1.gt.noabn).or.
     3     (h2new1.gt.noabn)) goto 111

       ihashold = h2new1-1+noabn*(h1new1-1+noabn *(p2new1-noabn-1+nvabn
     1 *(p1new1-noabn-1)))
 
        ilength = int_mb(k_r2_offsetm(jref))
        totaloff = -1

        do k1 = 1, ilength
          if(int_mb(k_r2_offsetm(jref)+k1).eq.ihashold) then
             totaloff = 1
            goto 114
          endif
        enddo

 114     continue


      if(totaloff.ne.-1)then

         ioff1 = p1off1*int_mb(k_rangem(jref)+h2new1-1)*
     1 int_mb(k_rangem(jref)+h1new1-1)*
     2 int_mb(k_rangem(jref)+p2new1-1)+
     2 p2off1*int_mb(k_rangem(jref)+h2new1-1)*
     3 int_mb(k_rangem(jref)+h1new1-1)+
     4 h1off1*int_mb(k_rangem(jref)+h2new1-1)+h2off1

       sizenew1 = int_mb(k_rangem(jref)+p1new1-1)*
     1 int_mb(k_rangem(jref)+p2new1-1)*
     2 int_mb(k_rangem(jref)+h1new1-1)*
     3 int_mb(k_rangem(jref)+h2new1-1)

       l_tmp3=0
      dsmult=-1.0d0*dsmult

        if(ioff1.ne.ioff) then
      if (.not.ma_push_get(mt_dbl,sizenew1,'tmp1',l_tmp3,k_tmp3))
     1       call errquit('tce_uss: MA problem',3,MA_ERR)

        CALL DFILL(sizenew1,0.0d0,dbl_mb(k_tmp3),1)

       if(dbl_mb(k_tmp3+ioff1).eq.0)then
       dbl_mb(k_tmp3+ioff1)=dbl_mb(k_r3+counter-1)*dsmult
     1*p2f(h1new1,h2new1)*p2f(p1new1,p2new1)/(p3f(t_p4b,t_p5b,t_p6b)*
     1 p3f(t_h1b,t_h2b,t_h3b)*fact)
       endif

      call add_hash_block(d_c2,dbl_mb(k_tmp3),sizenew1,
     1   int_mb(k_r2_offsetm(jref)),h2new1-1+noabn*(h1new1-1+noabn *
     1   (p2new1-noabn-1+nvabn*(p1new1-noabn-1))))

      if (.not.ma_pop_stack(l_tmp3))
     1    call errquit('tce_uss: MA problem',4,MA_ERR)
         endif
        endif !totaloff
        endif
 
       if(wth.ge.2) then
       p1new1=p1new
       p2new1=p2new
       h1new1=h2new
       h2new1=h1new
       p1off1=p1off
       p2off1=p2off
       h1off1=h2off
       h2off1=h1off

      if(p1new1.gt.p2new1) then
      t = p1new1
      p1new1 = p2new1
      p2new1 = t
      t = p1off1
      p1off1 = p2off1
      p2off1 = t
      dsmult=-1.0d0*dsmult
      end if

      if(h1new1.gt.h2new1) then
      t = h1new1
      h1new1 = h2new1
      h2new1 = t
      t = h1off1
      h1off1 = h2off1
      h2off1 = t
      dsmult=-1.0d0*dsmult
      end if
        if((p1new1.le.noabn).or.
     1     (p2new1.le.noabn).or.
     2     (h1new1.gt.noabn).or.
     3     (h2new1.gt.noabn)) goto 111

       ihashold = h2new1-1+noabn*(h1new1-1+noabn *(p2new1-noabn-1+nvabn
     1 *(p1new1-noabn-1)))

        ilength = int_mb(k_r2_offsetm(jref))
        totaloff = -1

        do k1 = 1, ilength
          if(int_mb(k_r2_offsetm(jref)+k1).eq.ihashold) then
             totaloff = 1
            goto 115
          endif
        enddo

 115     continue


      if(totaloff.ne.-1)then

         ioff2 = p1off1*int_mb(k_rangem(jref)+h2new1-1)*
     1 int_mb(k_rangem(jref)+h1new1-1)*
     2 int_mb(k_rangem(jref)+p2new1-1)+
     2 p2off1*int_mb(k_rangem(jref)+h2new1-1)*
     3 int_mb(k_rangem(jref)+h1new1-1)+
     4 h1off1*int_mb(k_rangem(jref)+h2new1-1)+h2off1

       sizenew1 = int_mb(k_rangem(jref)+p1new1-1)*
     1 int_mb(k_rangem(jref)+p2new1-1)*
     2 int_mb(k_rangem(jref)+h1new1-1)*
     3 int_mb(k_rangem(jref)+h2new1-1)

       l_tmp3=0
      dsmult=-1.0d0*dsmult

        if(ioff2.ne.ioff) then
      if (.not.ma_push_get(mt_dbl,sizenew1,'tmp1',l_tmp3,k_tmp3))
     1       call errquit('tce_uss: MA problem',3,MA_ERR)

        CALL DFILL(sizenew1,0.0d0,dbl_mb(k_tmp3),1)

       dbl_mb(k_tmp3+ioff2)=dbl_mb(k_r3+counter-1)*dsmult
     1*p2f(h1new1,h2new1)*p2f(p1new1,p2new1)/(p3f(t_p4b,t_p5b,t_p6b)*
     1 p3f(t_h1b,t_h2b,t_h3b)*fact)

      call add_hash_block(d_c2,dbl_mb(k_tmp3),sizenew1,
     1   int_mb(k_r2_offsetm(jref)),h2new1-1+noabn*(h1new1-1+noabn *
     1   (p2new1-noabn-1+nvabn*(p1new1-noabn-1))))

      if (.not.ma_pop_stack(l_tmp3))
     1    call errquit('tce_uss: MA problem',4,MA_ERR)
         endif
        endif !totaloff
       endif
 
       endif
        endif!all same spin
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcc
       endif!erank=4

 111  continue 

csb sb sb sb sb sb sb sb 
        end do
        end do
        end do
        end do
        end do
        end do

        if (.not.ma_pop_stack(l_r3))
     1  call errquit('tce_uss: MA problem',2,MA_ERR)

      endif
      END IF
      END IF
      END IF
cc    if(lusesub) then
cc     next = NXTASKsub(nprocs,1,mypgid)
cc    else
c     next = NXTASK(nprocs,1)
cc    endif
      next = NXTASK(nprocs,1)
      END IF
      count = count + 1

      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
c      if (nodezero) call util_flush(LuOut)
c     if(lusesub) then
c      next = NXTASKsub(-nprocs,1,mypgid)
c      call GA_pgroup_SYNC(mypgid)
c     else
      next = NXTASK(-nprocs,1)
      call GA_SYNC()
cc    endif

c      call ccsd_r3(dbl_mb(k_r3),d_f1,d_i1_1,
c    1d_i1_2,d_i1_3,
c    1d_t1,d_t2,d_v2,k_f1_offset,
c    1k_i1_offset_1,k_i1_offset_2,
c    1k_i1_offset_3,k_t1_offset,k_t2_offset,
c    1k_v2_offset,l_i1_offset_1,l_i1_offset_2,
c    1l_i1_offset_3,t_h1b,t_h2b,t_h3b,t_p4b,t_p5b,t_p6b,3)


       return
       end
!
      function p3f(i,j,k)
 
      integer i, j, k
      double precision p3f
      if((i.eq.j).and.(i.eq.k)) then
      p3f=6.0
      else
      if((i.eq.j).or.(j.eq.k))then
      p3f=2.0
      else
      p3f=1.0
      endif
      endif
      return
      end
!

      function p2f(i,j)
 
      integer i, j
      double precision p2f
      if(i.eq.j) then
      p2f=2.0
      else
      p2f=1.0
      endif
      return
      end
 
      function h2f(i,j)
      integer i,j
      double precision h2f
      if(i.eq.j) then
      h2f=2.0
      else
      h2f=1.0
      endif
      return
      end
 
c
