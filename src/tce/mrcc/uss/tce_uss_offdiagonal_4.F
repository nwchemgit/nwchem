!
      subroutine tce_uss_offdiagonal_4(d_t2,d_v2,
     2                 k_t2_offset,
     2                 k_v2_offset,d_c1,d_c2,
     2                 iref,jref,k_r1_offsetm,k_r2_offsetm)

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
ccccccccc
      integer d_v2,d_t2
      integer k_t2_offset,k_v2_offset
      integer d_i1_1,d_i1_2
      integer k_i1_offset_1,k_i1_offset_2
      integer l_i1_offset_1,l_i1_offset_2
ccccccccc
      integer iref,jref
      integer i,j,p1,h2,k,p2,h3,h4,i1,p,k1,k2
      integer size,l,m,n,o
      integer l_r4a,k_r4a
c
      integer t_h1b,t_h2b,t_h3b,t_h4b,t_p5b,t_p6b,t_p7b,t_p8b
c
      integer orbindex(8),aorbindex(8),orbindexnew(8)
      integer t, p1new,p2new, h1new, h2new,p3new,p4new,h3new,h4new
      integer p1new1,p2new1, h1new1, h2new1,p3new1,p4new1,h3new1,h4new1
      integer orbspin(8),aorbspin(8)
      integer ioccnew(maxorb,2),iocc0(maxorb,2)
      integer ioffset(6),ihash,oldhash,oldhash2
      integer p1off,p2off,h1off,h2off,p3off,p4off,h3off,h4off
      integer p1off1,p2off1,h1off1,h2off1,p3off1,p4off1,h3off1,h4off1
      integer ispinfrom,ispinto
      integer totaloff,ilength
      integer l_tmp,k_tmp,l_tmp4,k_tmp4
      integer l_d_c2,k_d_c2
      integer sizenew,ihashold,ioff,ioff1,ioff2,sizenew1
      integer noabn,nvabn,counter,iactive
      integer d_c1,d_c2
      EXTERNAL NXTASKsub
      EXTERNAL NXTASK
      INTEGER NXTASKsub
      INTEGER NXTASK
      INTEGER nprocs
      INTEGER count,next
      integer iexfrom,iexto,iexspin,wtp,wth
      integer noper,erank,optyp,eoper,signfact,noper2,eoper2,espin
      INTEGER logtoint
      INTEGER holep,particlep
      double precision dsmult,fact,p4f,p2f
      dimension eoper(4*maxexcit),optyp(4*maxexcit),eoper2(4*maxexcit),
     1 espin(4*maxexcit)
      dimension iexfrom(8),iexto(8),iexspin(8)

c     write(100,*)"REFERENCE",iref,jref
c     if(lusesub) then

c     call ga_pgroup_sync(mypgid)
c     nprocs = GA_pgroup_NNODES(mypgid)
c     count = 0
c     next = NXTASKsub(nprocs, 1,mypgid)

c     else
c
      call ccsd_r4(dbl_mb(k_r4a),d_i1_1,d_i1_2,d_t2,d_v2,
     &k_i1_offset_1,k_i1
     &_offset_2,k_t2_offset,k_v2_offset,l_i1_offset_1,l_i1_offset_2,t_h1
     &b,t_h2b,t_h3b,t_h4b,t_p5b,t_p6b,t_p7b,t_p8b,1)

c
      call ga_sync()
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)

c     endif
      nodezero = (ga_nodeid().eq.0)
c
c       if(nodezero)write(*,*)'out of 4'
c          
         noa = nblcks(1,iref)
         nob = nblcks(2,iref)
         nva = nblcks(3,iref)
         nvb = nblcks(4,iref)

         noab = noa+nob
         nvab = nva+nvb

         noabn = nblcks(1,jref)+nblcks(2,jref)
         nvabn = nblcks(3,jref)+nblcks(4,jref)
 
!
      DO t_p5b = noab+1,noab+nvab
      DO t_p6b = t_p5b,noab+nvab
      DO t_p7b = t_p6b,noab+nvab
      DO t_p8b = t_p7b,noab+nvab
      DO t_h1b = 1,noab
      DO t_h2b = t_h1b,noab
      DO t_h3b = t_h2b,noab
      DO t_h4b = t_h3b,noab


       if(next.eq.count) then
      IF (int_mb(k_spinm(iref)+t_p5b-1)+int_mb(k_spinm(iref)+t_p6b-1)+
     &int_mb(k_spinm(iref)+t_p7b-1)+int_mb(k_spinm(iref)+t_p8b-1).eq.
     &int_mb(k_spinm(iref)+t_h1b-1)+int_mb(k_spinm(iref)+t_h2b-1)+
     &int_mb(k_spinm(iref)+t_h3b-1)+int_mb(k_spinm(iref)+t_h4b-1))THEN
      IF(ieor(int_mb(k_symm(iref)+t_p5b-1),
     &ieor(int_mb(k_symm(iref)+t_p6b-1),ieor(int_mb(k_symm(iref)+t_p7b
     &-1),ieor(int_mb(k_symm(iref)+t_p8b-1),
     &ieor(int_mb(k_symm(iref)+t_h1b-1),
     &ieor(int_mb(k_symm(iref)+t_h2b-1),
     &ieor(int_mb(k_symm(iref)+t_h3b-1),
     &int_mb(k_symm(iref)+t_h4b-1)))))))) .eq. irrep_t) THEN
      IF ((.not.restricted).or.(int_mb(k_spinm(iref)+t_p5b-1)+
     &int_mb(k_spinm(iref)+t_p6b-1
     &)+int_mb(k_spinm(iref)+t_p7b-1)+int_mb(k_spinm(iref)+t_p8b-1)+
     &int_mb(k_spinm(iref)+t_h1b-1)+int_mb(k_spinm(iref)+t_h2b-1)+
     &int_mb(k_spinm(iref)+t_h3b-1)+
     &int_mb(k_spinm(iref)+t_h4b-1).ne.16)) THEN
!
      particlep= logtoint(log_mb(k_active+t_p5b-1))+
     1  logtoint(log_mb(k_active+t_p6b-1))+
     1  logtoint(log_mb(k_active+t_p7b-1))+
     1  logtoint(log_mb(k_active+t_p8b-1))
       holep=logtoint(log_mb(k_active+t_h1b-1))+
     1  logtoint(log_mb(k_active+t_h2b-1))+
     1  logtoint(log_mb(k_active+t_h3b-1))+
     1  logtoint(log_mb(k_active+t_h4b-1))
      if((holep.eq.2).and.(particlep.eq.2))then

      size = int_mb(k_rangem(iref)+t_p5b-1)*int_mb(k_rangem(iref)
     &+t_p6b-1)*int_mb(k_rangem(iref)+t_p7b-1)*int_mb(k_rangem(iref)
     &+t_p8b-1)*int_mb(k_rangem(iref)+t_h1b-1)
     &* int_mb(k_rangem(iref)+t_h2b-1)*int_mb(k_rangem(iref)+t_h3b-1)* 
     &int_mb(k_rangem(iref)+t_h4b-1)
!
        if (.not.ma_push_get(mt_dbl,size,'c4',l_r4a,k_r4a))
     1   call errquit('tce_uss: MA problem',10,MA_ERR)
!
        CALL DFILL(size,0.0d0,dbl_mb(k_r4a),1)
!
       call ccsd_r4(dbl_mb(k_r4a),d_i1_1,d_i1_2,d_t2,d_v2,
     &k_i1_offset_1,k_i1_offset_2,k_t2_offset,k_v2_offset,
     &l_i1_offset_1,l_i1_offset_2,t_h1
     &b,t_h2b,t_h3b,t_h4b,t_p5b,t_p6b,t_p7b,t_p8b,2)
!

        counter =0
        l_tmp4=0
        ihash =-1
        ihashold =-1
!
       do i=1,int_mb(k_rangem(iref)+t_p5b-1)
       do j=1,int_mb(k_rangem(iref)+t_p6b-1)
       do k=1,int_mb(k_rangem(iref)+t_p7b-1)
       do l=1,int_mb(k_rangem(iref)+t_p8b-1)
       do m=1,int_mb(k_rangem(iref)+t_h1b-1)
       do n=1,int_mb(k_rangem(iref)+t_h2b-1)
       do o=1,int_mb(k_rangem(iref)+t_h3b-1)
       do p=1,int_mb(k_rangem(iref)+t_h4b-1)
c
       counter = counter +1
c
!
      orbspin(1) = int_mb(k_spinm(iref)+t_p5b-1) -1
      orbspin(2) = int_mb(k_spinm(iref)+t_p6b-1) -1
      orbspin(3) = int_mb(k_spinm(iref)+t_p7b-1) -1
      orbspin(4) = int_mb(k_spinm(iref)+t_p8b-1) -1
      orbspin(5) = int_mb(k_spinm(iref)+t_h1b-1) -1
      orbspin(6) = int_mb(k_spinm(iref)+t_h2b-1) -1
      orbspin(7) = int_mb(k_spinm(iref)+t_h3b-1) -1
      orbspin(8) = int_mb(k_spinm(iref)+t_h4b-1) -1
       
!
       orbindex(1) = (1 - orbspin(1)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offsetm(iref)+t_p5b-1)+i-1))/2
!
       orbindex(2) = (1-orbspin(2)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offsetm(iref)+t_p6b-1)+j-1))/2
!
       orbindex(3) = (1 - orbspin(3)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offsetm(iref)+t_p7b-1)+k-1))/2
!
       orbindex(4) = (1-orbspin(4)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offsetm(iref)+t_p8b-1)+l-1))/2
!
       orbindex(5) = (1-orbspin(5)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offsetm(iref)+t_h1b-1)+m-1))/2
!
       orbindex(6) = (1-orbspin(6)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offsetm(iref)+t_h2b-1)+n-1))/2
!
       orbindex(7) = (1-orbspin(7)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offsetm(iref)+t_h3b-1)+o-1))/2
!
       orbindex(8) = (1-orbspin(8)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offsetm(iref)+t_h4b-1)+p-1))/2
!
!
!sb sb sb sb sb sb
      iexto(1) = moindexes(orbindex(1),orbspin(1)+1,iref)
      iexto(2) = moindexes(orbindex(2),orbspin(2)+1,iref)
      iexto(3) = moindexes(orbindex(3),orbspin(3)+1,iref)
      iexto(4) = moindexes(orbindex(4),orbspin(4)+1,iref)
      iexfrom(1) = moindexes(orbindex(5),orbspin(5)+1,iref)
      iexfrom(2) = moindexes(orbindex(6),orbspin(6)+1,iref)
      iexfrom(3) = moindexes(orbindex(7),orbspin(7)+1,iref)
      iexfrom(4) = moindexes(orbindex(8),orbspin(8)+1,iref)
c
      iexspin(1) = orbspin(1)+1
      iexspin(2) = orbspin(2)+1
      iexspin(3) = orbspin(3)+1
      iexspin(4) = orbspin(4)+1
!
       signfact=0
       dsmult=1.0d0
c
       call perfexcit(erank,eoper,iref,jref,
     &   4,iexfrom,iexto,iexspin,
     &   signfact,espin,wtp,wth)
c

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


      ihash = h1new-1+noabn*(p1new-noabn-1)

        ilength = int_mb(k_r1_offsetm(jref))
        totaloff = -1

        do k1 = 1, ilength
          if(int_mb(k_r1_offsetm(jref)+k1).eq.ihash) then
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

      ioff = p1off*int_mb(k_rangem(jref)+h1new-1)+h1off

      sizenew = int_mb(k_rangem(jref)+h1new-1)*
     1           int_mb(k_rangem(jref)+p1new-1)
!
      if (.not.ma_push_get(mt_dbl,sizenew,'tmp',l_tmp4,k_tmp4))
     1       call errquit('tce_uss: MA problem',3,MA_ERR)

        CALL DFILL(sizenew,0.0d0,dbl_mb(k_tmp4),1)
!
      dbl_mb(k_tmp4+ioff)=dbl_mb(k_r4a+counter-1)*dsmult/
     1 (p4f(t_p5b,t_p6b,t_p7b,t_p8b)*p4f(t_h1b,t_h2b,t_h3b,t_h4b))

      call add_hash_block(d_c1,dbl_mb(k_tmp4),sizenew,
     1   int_mb(k_r1_offsetm(jref)),ihash)
!
       if (.not.ma_pop_stack(l_tmp4))
     1    call errquit('tce_uss: MA problem',4,MA_ERR)

         endif !totaloff
!
       endif!erank=2


       if(erank.eq.4) then
c
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
     1  .and.(((wtp.ge.3).and.(p1new.eq.p2new)).or.((wth.ge.3).and.
     1  (h1new.eq.h2new)))) then
        fact=2.d0
       else
       fact=1.d0
       endif

      ihash = h2new-1+noabn*(h1new-1+noabn *(p2new-noabn-1+nvabn
     1 *(p1new-noabn-1)))

        ilength = int_mb(k_r2_offsetm(jref))
        totaloff = -1

        do k1 = 1, ilength
          if(int_mb(k_r2_offsetm(jref)+k1).eq.ihash) then
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


      if (.not.ma_push_get(mt_dbl,sizenew,'tmp1',l_tmp4,k_tmp4))
     1       call errquit('tce_uss: MA problem',3,MA_ERR)

        CALL DFILL(sizenew,0.0d0,dbl_mb(k_tmp4),1)

      dbl_mb(k_tmp4+ioff)=dbl_mb(k_r4a+counter-1)*dsmult
     1 *p2f(p1new,p2new)*p2f(h1new,h2new)/
     1 (fact*p4f(t_p5b,t_p6b,t_p7b,t_p8b)*p4f(t_h1b,t_h2b,t_h3b,t_h4b))

      call add_hash_block(d_c2,dbl_mb(k_tmp4),sizenew,
     1   int_mb(k_r2_offsetm(jref)),ihash)

      if (.not.ma_pop_stack(l_tmp4))
     1    call errquit('tce_uss: MA problem',4,MA_ERR)
        endif !totaloff
ccccc
       if((eoper(1).ne.eoper(2)).and.(eoper(3).ne.eoper(4)))then
       if((espin(1).eq.espin(2)).and.(espin(3).eq.espin(4))) then

        dsmult=-1.0d0*dsmult
       if (wtp.ge.3)then
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
      end if

      if(h1new1.gt.h2new1) then
      t = h1new1
      h1new1 = h2new1
      h2new1 = t
      t = h1off1
      h1off1 = h2off1
      h2off1 = t
      end if

        if((p1new1.le.noabn).or.
     1     (p2new1.le.noabn).or.
     2     (h1new1.gt.noabn).or.
     3     (h2new1.gt.noabn)) goto 111
c
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


       ioff1=0
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

        l_tmp4=0

        if(ioff1.ne.ioff) then
      if (.not.ma_push_get(mt_dbl,sizenew1,'tmp1',l_tmp4,k_tmp4))
     1       call errquit('tce_uss: MA problem',3,MA_ERR)

        CALL DFILL(sizenew1,0.0d0,dbl_mb(k_tmp4),1)

      dbl_mb(k_tmp4+ioff1)=dbl_mb(k_r4a+counter-1)*dsmult
     1 *p2f(p1new,p2new)*p2f(h1new,h2new)/
     1 (fact*p4f(t_p5b,t_p6b,t_p7b,t_p8b)*p4f(t_h1b,t_h2b,t_h3b,t_h4b))


      call add_hash_block(d_c2,dbl_mb(k_tmp4),sizenew1,
     1   int_mb(k_r2_offsetm(jref)),h2new1-1+noabn*(h1new1-1+noabn *
     1   (p2new1-noabn-1+nvabn*(p1new1-noabn-1))))

      if (.not.ma_pop_stack(l_tmp4))
     1    call errquit('tce_uss: MA problem',4,MA_ERR)
         endif
        endif !totaloff
       endif
ccc
       if (wth.ge.3)then
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
      end if

      if(h1new1.gt.h2new1) then
      t = h1new1
      h1new1 = h2new1
      h2new1 = t
      t = h1off1
      h1off1 = h2off1
      h2off1 = t
      end if

        if((p1new1.le.noabn).or.
     1     (p2new1.le.noabn).or.
     2     (h1new1.gt.noabn).or.
     3     (h2new1.gt.noabn)) goto 111
c
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

       ioff1=0
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

        l_tmp4=0

        if(ioff2.ne.ioff) then
      if (.not.ma_push_get(mt_dbl,sizenew1,'tmp1',l_tmp4,k_tmp4))
     1       call errquit('tce_uss: MA problem',3,MA_ERR)

        CALL DFILL(sizenew1,0.0d0,dbl_mb(k_tmp4),1)

      dbl_mb(k_tmp4+ioff2)=dbl_mb(k_r4a+counter-1)*dsmult
     1 *p2f(p1new,p2new)*p2f(h1new,h2new)/
     1 (fact*p4f(t_p5b,t_p6b,t_p7b,t_p8b)*p4f(t_h1b,t_h2b,t_h3b,t_h4b))


      call add_hash_block(d_c2,dbl_mb(k_tmp4),sizenew1,
     1   int_mb(k_r2_offsetm(jref)),h2new1-1+noabn*(h1new1-1+noabn *
     1   (p2new1-noabn-1+nvabn*(p1new1-noabn-1))))

      if (.not.ma_pop_stack(l_tmp4))
     1    call errquit('tce_uss: MA problem',4,MA_ERR)
         endif
        endif !totaloff
       endif

       endif
       endif

       endif!erank=4

111    continue
!sb sb sb sb sb sb

       end do
       end do
       end do
       end do
       end do
       end do
       end do
       end do
        if (.not.ma_pop_stack(l_r4a))
     1  call errquit('tce_uss: MA problem',2,MA_ERR)

      endif
      END IF
      END IF
      END IF
c     if(lusesub) then
c      next = NXTASKsub(nprocs,1,mypgid)
c     else
       next = NXTASK(nprocs, 1)
c     endif
      END IF
      count = count + 1

      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO

c     if(lusesub) then
c      next = NXTASKsub(-nprocs,1,mypgid)
c      call GA_pgroup_SYNC(mypgid)
c     else
       next = NXTASK(-nprocs, 1)
       call GA_SYNC()
c     endif

       call ccsd_r4(dbl_mb(k_r4a),d_i1_1,d_i1_2,d_t2,d_v2,
     &k_i1_offset_1,k_i1_offset_2,k_t2_offset,k_v2_offset,
     &l_i1_offset_1,l_i1_offset_2,t_h1
     &b,t_h2b,t_h3b,t_h4b,t_p5b,t_p6b,t_p7b,t_p8b,3)
c


       return
       end
!

       function p4f(i,j,k,l)
       integer i,j,k,l
       double precision p4f
       if((i.eq.j).and.(i.eq.k).and.(i.eq.l)) then
       p4f=24.0
       else
       if(((i.eq.j).and.(i.eq.k)).or.((j.eq.k).and.(j.eq.l)))then
       p4f=6.0
       else
       if((i.eq.j).and.(k.eq.l))then
       p4f=4.0
       else
       if((i.eq.j).or.(j.eq.k).or.(k.eq.l)) then
       p4f=2.0
       else
       p4f=1.0
       endif
       endif
       endif
       endif
       return
       end
!
       function h4f(i,j,k,l)
       integer i,j,k,l
       double precision h4f
       if((i.eq.j).and.(i.eq.k).and.(i.eq.l)) then
       h4f=24.0
       else
       if(((i.eq.j).and.(i.eq.k)).or.((j.eq.k).and.(j.eq.l)))then
       h4f=6.0
       else
       if((i.eq.j).or.(j.eq.k).or.(k.eq.l)) then
       h4f=2.0
       else
       if((i.eq.j).and.(k.eq.l))then
       h4f=4.0
       else
       h4f=1.0
       endif
       endif
       endif
       endif
       return
       end

