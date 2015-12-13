      subroutine tce_uss_offdiagonal_2(d_r1m,d_r2m,k_r1_offsetm,
     1 k_r2_offsetm,iref,jref,d_c2,d_c1)
! Routine for off diagonal correction
!
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
      integer k_r2_offsetm(maxref),k_r2a_offsetm(maxref)
      integer k_r1_offsetm(maxref)
      integer d_r2m(maxref),d_r2a(maxref),d_r1m(maxref)
      integer iref,jref,iexclevel
      integer i,j,p1,h2,k,p2,h3,h4,i1,k1,k2
      integer size,l,m,n,o
      integer l_r2,k_r2,l_r2a,k_r2a
      integer p1b,p2b,h4b,h2b,h1b,p3b,p3,h3b
      integer orbindex(8),aorbindex(8),orbindexnew(8)
      integer ahblock(4),apblock(4)
      integer t, p1new,p2new, h1new, h2new,p3new,p4new,h3new,h4new
      integer p1new1,p2new1, h1new1, h2new1,p3new1,p4new1,h3new1,h4new1
      integer orbspin(8),aorbspin(8),actspin(2)
      integer ioccnew(maxorb,2),iocc0(maxorb,2)
      integer ioffset(6),ihash,oldhash,oldhash2
      integer p1off,p2off,h1off,h2off,p3off,p4off,h3off,h4off
      integer p1off1,p2off1,h1off1,h2off1,p3off1,p4off1,h3off1,h4off1
      integer ispinfrom,ispinto
      integer iu,is,ifrom,ito
      integer totaloff,hs,ilength
      integer l_tmp2,k_tmp2
      integer l_d_c2,k_d_c2
      integer acth1,acth2,actp1,actp2
      integer sizenew,ihashold,ioff,ioff1,ioff2,sizenew1
      integer noabn,nvabn,counter,iactive
      integer d_c2,d_c1
      logical isfound
      logical ap1,ah2,ap2,ah1
      EXTERNAL NXTASKsub
      EXTERNAL NXTASK
      INTEGER NXTASKsub
      INTEGER NXTASK
      INTEGER nxt
      INTEGER nprocs
      INTEGER count,next
      integer iexfrom,iexto,iexspin
      integer noper,erank,optyp,eoper,signfact,noper2,eoper2,espin
      integer wtp, wth
      double precision dsmult
      double precision p2f,h2f,fact
      dimension eoper(4*maxexcit),optyp(4*maxexcit),eoper2(4*maxexcit),
     1 espin(4*maxexcit)
      dimension iexfrom(8),iexto(8),iexspin(8)
c
      if(lusesub) then

      call ga_pgroup_sync(mypgid)
      nprocs = GA_pgroup_NNODES(mypgid)
      count = 0
      next = NXTASKsub(nprocs, 1,mypgid)

      else

      call ga_sync()
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)

      endif



      nodezero = (ga_nodeid().eq.0)
c       if(nodezero)write(*,*)'out of 2'
!po           
c       if(nodezero)write(100,*)"REFEREENCE",iref,jref
         noa = nblcks(1,iref)
         nob = nblcks(2,iref)
         nva = nblcks(3,iref)
         nvb = nblcks(4,iref)

         noab = noa+nob
         nvab = nva+nvb

         noabn = nblcks(1,jref)+nblcks(2,jref)
         nvabn = nblcks(3,jref)+nblcks(4,jref)
 
c-------
c   R(i,a,x,y)
c-------
! 
      DO p1b = noab+1,noab+nvab
      DO p2b = p1b,noab+nvab
      DO h1b = 1,noab
      DO h2b = h1b,noab
        if(next.eq.count) then
!
      IF (int_mb(k_spinm(iref)+p1b-1)+int_mb(k_spinm(iref)+p2b-1)
     1 .eq.int_mb(k_spinm(iref)+h1b-1)+int_mb(k_spinm(iref)+h2b-1))THEN
!
      IF (ieor(int_mb(k_symm(iref)+p1b-1),ieor(int_mb(k_symm(iref)+p2b
     1 -1),ieor(int_mb(k_symm(iref)+h1b-1),int_mb(k_symm(iref)+h2b-1))))
     1  .eq. irrep_t) THEN
!
      IF ((.not.restricted).or.(int_mb(k_spinm(iref)+p1b-1)+
     1 int_mb(k_spinm(iref)+p2b-1)+int_mb(k_spinm(iref)+h1b-1)
     1 +int_mb(k_spinm(iref)+h2b-1).ne.8)) THEN
!
      size = int_mb(k_rangem(iref)+p1b-1) * int_mb(k_rangem(iref)+p2b-1)
     1 * int_mb(k_rangem(iref)+h1b-1) * int_mb(k_rangem(iref)+h2b-1)
!
!
      oldhash = h2b-1+noab*(h1b-1+noab*(p2b-
     &noab-1+nvab*(p1b - noab - 1)))
!
        if (.not.ma_push_get(mt_dbl,size,'c2',l_r2a,k_r2a))
     1   call errquit('tce_uss: MA problem',10,MA_ERR)
!
        CALL DFILL(size,0.0d0,dbl_mb(k_r2a),1)
!
        call get_hash_block(d_r2m(iref),dbl_mb(k_r2a),size,
     1   int_mb(k_r2_offsetm(iref)),oldhash)
!-----------------------------------!
        counter = 0
        l_tmp2 = 0
        ihash =-1
        ihashold = -1
        do i=1,int_mb(k_rangem(iref)+p1b-1)
        do j=1,int_mb(k_rangem(iref)+p2b-1)
        do m=1,int_mb(k_rangem(iref)+h1b-1)
        do n=1,int_mb(k_rangem(iref)+h2b-1)
!
        counter = counter+1
!
      orbspin(1) = int_mb(k_spinm(iref)+p1b-1) -1
      orbspin(2) = int_mb(k_spinm(iref)+p2b-1) -1
      orbspin(3) = int_mb(k_spinm(iref)+h1b-1) -1
      orbspin(4) = int_mb(k_spinm(iref)+h2b-1) -1
!
       orbindex(1) = (1 - orbspin(1)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offsetm(iref)+p1b-1)+i-1))/2
!
       orbindex(2) = (1-orbspin(2)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offsetm(iref)+p2b-1)+j-1))/2
!
       orbindex(3) = (1 - orbspin(3)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offsetm(iref)+h1b-1)+m-1))/2
!
       orbindex(4) = (1-orbspin(4)+
     1 int_mb(k_mo_indexm(iref)+int_mb(k_offsetm(iref)+h2b-1)+n-1))/2
!
      iexto(1) = moindexes(orbindex(1),orbspin(1)+1,iref)
      iexto(2) = moindexes(orbindex(2),orbspin(2)+1,iref)
      iexfrom(1) = moindexes(orbindex(3),orbspin(3)+1,iref)
      iexfrom(2) = moindexes(orbindex(4),orbspin(4)+1,iref)
      iexspin(1) = orbspin(3)+1
      iexspin(2) = orbspin(4)+1
!
       signfact=0
       dsmult=1.0d0

       call perfexcit(erank,eoper,iref,jref,
     &   2,iexfrom,iexto,iexspin,
     &   signfact,espin,wtp,wth)
!
      if(erank.eq.0) goto 111
      if(erank.eq.2) then
c
       signfact=mod(signfact,2)
       if(signfact.ne.0) then
        dsmult = -1.0d0
       endif
c
      orbindex(1)=moindexes(eoper(1),espin(1),jref)
      orbindex(2)=moindexes(eoper(2),espin(2),jref)

      p1new = orbinblck(orbindex(1),espin(1),jref)
      h1new = orbinblck(orbindex(2),espin(2),jref)

      p1off = offsetinblck(orbindex(1),espin(1),jref)
      h1off = offsetinblck(orbindex(2),espin(2),jref)

      ihash = h1new-1+noabn*(p1new-noabn-1)

        ilength = int_mb(k_r1_offsetm(jref))
        totaloff = -1

        do k = 1, ilength
          if(int_mb(k_r1_offsetm(jref)+k).eq.ihash) then
             totaloff = 1
            goto 112
          endif
        enddo
 112     continue

c      call tce_hash_n(int_mb(k_r1_offsetm(jref)),ihash,totaloff)

        if((p1new.le.noabn).or.
     2     (h1new.gt.noabn))then
           totaloff=-1
        endif

      if(totaloff.ne.-1)then

      ioff = p1off*int_mb(k_rangem(jref)+h1new-1)+h1off

      sizenew = int_mb(k_rangem(jref)+h1new-1)*
     1           int_mb(k_rangem(jref)+p1new-1)
!
      if (.not.ma_push_get(mt_dbl,sizenew,'tmp',l_tmp2,k_tmp2))
     1       call errquit('tce_uss: MA problem',3,MA_ERR)

        CALL DFILL(sizenew,0.0d0,dbl_mb(k_tmp2),1)
!
      dbl_mb(k_tmp2+ioff)=dbl_mb(k_r2a+counter-1)*dsmult/
     1(p2f(p1b,p2b)*h2f(h1b,h2b))

      call add_hash_block(d_c1,dbl_mb(k_tmp2),sizenew,
     1   int_mb(k_r1_offsetm(jref)),ihash)
!
       if (.not.ma_pop_stack(l_tmp2))
     1    call errquit('tce_uss: MA problem',4,MA_ERR)

         endif !totaloff
 
       endif!erank=2
!
       if(erank.eq.4) then
c
       orbindex(1)=moindexes(eoper(1),espin(1),jref)
       orbindex(2)=moindexes(eoper(2),espin(2),jref)
       orbindex(4)=moindexes(eoper(3),espin(3),jref)
       orbindex(3)=moindexes(eoper(4),espin(4),jref)
 
       signfact=mod(signfact,2)
       if(signfact.ne.0) then
        dsmult = -1.0d0
       endif

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
c     if(nodezero)write(100,999)orbindex(1),orbindex(2),orbindex(3),
c    1orbindex(4),wtp,wth
c999   format(6I4)

       if(((espin(1).eq.espin(2)).and.(espin(3).eq.espin(4)))
     1  .and.(((wtp.ge.1).and.(p1new.eq.p2new)).or.((wth.ge.1).and.
     1  (h1new.eq.h2new)))) then
        fact=2.d0
       else
       fact=1.d0
       endif

       ihashold = h2new-1+noabn*(h1new-1+noabn *(p2new-noabn-1+nvabn
     1 *(p1new-noabn-1)))
 
        ilength = int_mb(k_r2_offsetm(jref))
        totaloff = -1

        do k = 1, ilength
          if(int_mb(k_r2_offsetm(jref)+k).eq.ihashold) then
             totaloff = 1
            goto 113
          endif
        enddo
 113     continue

c     call tce_hash_n(int_mb(k_r2_offsetm(jref)),ihashold,totaloff)
        if((p1new.le.noabn).or.
     1     (p2new.le.noabn).or.
     2     (h1new.gt.noabn).or.
     3     (h2new.gt.noabn)) then
           totaloff=-1
        endif 



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
 
      if (.not.ma_push_get(mt_dbl,sizenew,'tmp1',l_tmp2,k_tmp2))
     1       call errquit('tce_uss: MA problem',3,MA_ERR)

        CALL DFILL(sizenew,0.0d0,dbl_mb(k_tmp2),1)
 
      dbl_mb(k_tmp2+ioff)=dbl_mb(k_r2a+counter-1)*dsmult*
     1 p2f(p1new,p2new)*p2f(h1new,h2new)/(p2f(p1b,p2b)*p2f(h1b,h2b)*
     1 fact)

      call add_hash_block(d_c2,dbl_mb(k_tmp2),sizenew,
     1   int_mb(k_r2_offsetm(jref)),ihashold)
 
      if (.not.ma_pop_stack(l_tmp2))
     1    call errquit('tce_uss: MA problem',4,MA_ERR)
        endif !totaloff

c      if((eoper(1).ne.eoper(2)).and.
c    1  (eoper(3).ne.eoper(4)))then
       if((espin(1).eq.espin(2)).and.(espin(3).eq.espin(4)))then

         dsmult=-1.0d0*dsmult
       if(wtp.ge.1) then
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
c         dsmult=-1.0d0*dsmult
      end if

      if(h1new1.gt.h2new1) then
      t = h1new1
      h1new1 = h2new1
      h2new1 = t
      t = h1off1
      h1off1 = h2off1
      h2off1 = t
c         dsmult=-1.0d0*dsmult
      end if

        if((p1new1.le.noabn).or.
     1     (p2new1.le.noabn).or.
     2     (h1new1.gt.noabn).or.
     3     (h2new1.gt.noabn)) goto 111
c          totaloff=-1
c       endif

       ihashold = h2new1-1+noabn*(h1new1-1+noabn *(p2new1-noabn-1+nvabn
     1 *(p1new1-noabn-1)))

        ilength = int_mb(k_r2_offsetm(jref))
c       ioff = 0
        totaloff = -1

        do k = 1, ilength
          if(int_mb(k_r2_offsetm(jref)+k).eq.ihashold) then
             totaloff = 1
            goto 114
          endif
        enddo

 114     continue
c       call tce_hash_n(int_mb(k_r2_offsetm(jref)),ihashold,totaloff)


      if(totaloff.ne.-1)then

         ioff1 = p1off1*int_mb(k_rangem(jref)+h2new1-1)*
     1 int_mb(k_rangem(jref)+h1new1-1)*
     2 int_mb(k_rangem(jref)+p2new1-1)+
     2 p2off1*int_mb(k_rangem(jref)+h2new1-1)*
     3 int_mb(k_rangem(jref)+h1new1-1)+
     4 h1off1*int_mb(k_rangem(jref)+h2new1-1)+h2off1

       sizenew1 = int_mb(k_rangem(jref)+p1new-1)*
     1 int_mb(k_rangem(jref)+p2new-1)*
     2 int_mb(k_rangem(jref)+h1new-1)*
     3 int_mb(k_rangem(jref)+h2new-1)
       l_tmp2=0
 
        if(ioff1.ne.ioff) then
      if (.not.ma_push_get(mt_dbl,sizenew1,'tmp1',l_tmp2,k_tmp2))
     1       call errquit('tce_uss: MA problem',3,MA_ERR)

        CALL DFILL(sizenew1,0.0d0,dbl_mb(k_tmp2),1)

      dbl_mb(k_tmp2+ioff1)=dbl_mb(k_r2a+counter-1)*dsmult
     1 *p2f(p1new,p2new)*p2f(h1new,h2new)/(p2f(p1b,p2b)*p2f(h1b,h2b)*
     1 fact)

      call add_hash_block(d_c2,dbl_mb(k_tmp2),sizenew1,
     1   int_mb(k_r2_offsetm(jref)),h2new1-1+noabn*(h1new1-1+noabn *
     1   (p2new1-noabn-1+nvabn*(p1new1-noabn-1))))

      if (.not.ma_pop_stack(l_tmp2))
     1    call errquit('tce_uss: MA problem',4,MA_ERR)
         endif
       endif !totaloff
      endif
 
       if(wth.ge.1) then
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
c         dsmult=-1.0d0*dsmult
      end if

      if(h1new1.gt.h2new1) then
      t = h1new1
      h1new1 = h2new1
      h2new1 = t
      t = h1off1
      h1off1 = h2off1
      h2off1 = t
c         dsmult=-1.0d0*dsmult
      end if

        if((p1new1.le.noabn).or.
     1     (p2new1.le.noabn).or.
     2     (h1new1.gt.noabn).or.
     3     (h2new1.gt.noabn)) goto 111
c          totaloff=-1
c       endif

       ihashold = h2new1-1+noabn*(h1new1-1+noabn *(p2new1-noabn-1+nvabn
     1 *(p1new1-noabn-1)))

        ilength = int_mb(k_r2_offsetm(jref))
        totaloff = -1

        do k = 1, ilength
          if(int_mb(k_r2_offsetm(jref)+k).eq.ihashold) then
             totaloff = 1
            goto 115
          endif
        enddo

 115     continue
c       call tce_hash_n(int_mb(k_r2_offsetm(jref)),ihashold,totaloff)


      if(totaloff.ne.-1)then

         ioff2 = p1off1*int_mb(k_rangem(jref)+h2new1-1)*
     1 int_mb(k_rangem(jref)+h1new1-1)*
     2 int_mb(k_rangem(jref)+p2new1-1)+
     2 p2off1*int_mb(k_rangem(jref)+h2new1-1)*
     3 int_mb(k_rangem(jref)+h1new1-1)+
     4 h1off1*int_mb(k_rangem(jref)+h2new1-1)+h2off1

       sizenew1 = int_mb(k_rangem(jref)+p1new-1)*
     1 int_mb(k_rangem(jref)+p2new-1)*
     2 int_mb(k_rangem(jref)+h1new-1)*
     3 int_mb(k_rangem(jref)+h2new-1)
       l_tmp2=0

        if(ioff2.ne.ioff) then
      if (.not.ma_push_get(mt_dbl,sizenew1,'tmp1',l_tmp2,k_tmp2))
     1       call errquit('tce_uss: MA problem',3,MA_ERR)

        CALL DFILL(sizenew1,0.0d0,dbl_mb(k_tmp2),1)

      dbl_mb(k_tmp2+ioff2)=dbl_mb(k_r2a+counter-1)*dsmult
     1 *p2f(p1new,p2new)*p2f(h1new,h2new)/(p2f(p1b,p2b)*p2f(h1b,h2b)*
     1 fact)

      call add_hash_block(d_c2,dbl_mb(k_tmp2),sizenew1,
     1   int_mb(k_r2_offsetm(jref)),h2new1-1+noabn*(h1new1-1+noabn *
     1   (p2new1-noabn-1+nvabn*(p1new1-noabn-1))))

      if (.not.ma_pop_stack(l_tmp2))
     1    call errquit('tce_uss: MA problem',4,MA_ERR)
         endif
       endif !totaloff

       endif


       endif
c     endif

       endif!erank=4

111    continue
!sb sb sb sb sb sb sb sb
!
       end do
       end do
       end do
       end do
!
!
!
        if (.not.ma_pop_stack(l_r2a))
     1  call errquit('tce_uss: MA problem',2,MA_ERR)

        END IF!symmetry
        END IF
        END IF
      if(lusesub) then
       next = NXTASKsub(nprocs,1,mypgid)
      else
       next = NXTASK(nprocs, 1)
      endif
      END IF
      count = count + 1

      END DO
      END DO
      END DO
      END DO

      if(lusesub) then
       next = NXTASKsub(-nprocs,1,mypgid)
       call GA_pgroup_SYNC(mypgid)
      else
       next = NXTASK(-nprocs, 1)
       call GA_SYNC()
      endif

!
      return
      end
c
!
