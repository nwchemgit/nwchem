c $Id: amshift.f,v 1.1 1995-10-30 20:57:31 d3e129 Exp $
c----------------------------------------------------
C*
C*  THESE ROUTINES SHIFT THE ANGULAR MOMENTUM
C*
C*         FROM POSITION 1 TO POSITION 2
C*
C*                   AND
C*
C*         FROM POSITION 3 TO POSITION 4
c----------------------------------------------------
c    for re-ordered basis set 
c   nqi.ge.nqj  and  nqk.ge.nql
c 
c   other cases are not included here !
c----------------------------------------------------
      subroutine amshift(bl,nbls,l01,l02,npij,npkl,ngcd)
      implicit real*8 (a-h,o-z)
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
c
      common /cpu/ intsize,iacc,icache,memreal
cxx
      common /logic4/ nfu(1)
cxx
c
      COMMON/SHELL/LSHELLT,LSHELIJ,LSHELKL,LHELP,LCAS2(4),LCAS3(4)
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
c
ccccc common /big/ bl(1)
      common /memor4/ iwt0,iwt1,iwt2,ibuf,ibuf2,
     * ibfij1,ibfij2,ibfkl1,ibfkl2,
     * ibf2l1,ibf2l2,ibf2l3,ibf2l4,ibfij3,ibfkl3,
     * ibf3l,issss,
     * ix2l1,ix2l2,ix2l3,ix2l4,ix3l1,ix3l2,ix3l3,ix3l4,
     * ixij,iyij,izij, iwij,ivij,iuij,isij
c
      common /memor4a/ ibf3l1,ibf3l2,ibf3l3,ibf3l4
c
      common /memor5a/ iaa,ibb,icc,idd,icis,icjs,icks,icls,
     * ixab,ixp,ixpn,ixpp,iabnia,iapb,i1apb,ifij,icij,isab,
     * ixcd,ixq,ixqn,ixqq,icdnia,icpd,i1cpd,ifkl,ickl,iscd
      common /memor5b/ irppq,
     * irho,irr1,irys,irhoapb,irhocpd,iconst,ixwp,ixwq,ip1234,
     * idx1,idx2,indx
c
      dimension bl(*)
c
C**********************************
c* dimensions for "shifts"
c
      lsmx=max(lnij,lnkl)
      lsjl=max(nfu(nqj+1),nfu(nql+1))
c
      lqij=nfu(nqij+1)
      lqkl=nfu(nqkl+1)
      lqmx=max(lqij,lqkl)
c
c* memory for array INDXX :
c
      mindxx=lnijkl
      if(intsize.ne.1) mindxx=lnijkl/intsize+1
c
      call getmem(mindxx,indxx)
c
c------------------------------------------------
          nmr =ngcd*nbls*lnijkl 
c-
          nbuf2=ibuf2
          nbuf=ibuf
          m=1
          if(where.eq.'shif') then
            nbuf=ibuf+nmr
            m=6
          endif
c-
          mnbls=m*nbls
c-
c------------------------------------------------
c------------------------------------------------
C************************************************
C**   SPECIAL CASES WHERE THE SHIFTING OF ANGULAR
C**   MOMENTUM IS NOT NEEDED AT ALL / (DS|SS)..(SS|SD),
C**   (XS|YS),(XS|SY),(SX|YS),(SX|SY) /
C
      IF(NQIJ.EQ.NSIJ .AND. NQKL.EQ.NSKL) THEN
c
c         find apropriate matrices with l-shells
c
          ibfijx=ibfij1
          if(lshelij.eq.2) ibfijx=ibfij2
          ibfklx=ibfkl1
          if(lshelkl.eq.2) ibfklx=ibfkl2
          ibf2lx=ibf2l1
          if(lcas2(2).eq.1) ibf2lx=ibf2l2
          if(lcas2(3).eq.1) ibf2lx=ibf2l3
          if(lcas2(4).eq.1) ibf2lx=ibf2l4
c
c-
          incre =mnbls*lnijkl
          incre2=mnbls*l01*l02
c-
            do 100 iqu=1,ngcd
              jbuf =nbuf +(iqu-1)*incre
              jbuf2=nbuf2+(iqu-1)*incre2
              call noshift(l01,l02,mnbls,  bl(jbuf), bl(jbuf2),
     *          bl(ibfijx),           bl(ibfklx),
     *          bl(ibf2lx),
     *          lqij,lqkl,
     *          bl(indxx),lni*lnj,lnk*lnl)
  100       continue
         call retmem(1)
         return
      ENDIF
CC
ccccccccccccccccccccc
c
       call convr3(bl,m,nbls,npij,npkl,bl(idx1),bl(idx2),
     *             bl(ixab),bl(ixcd),ixabn,ixcdn)
c
c****
c
       if (lshellt.eq.0) then
c-
        incre =mnbls*lnijkl
        incre2=mnbls*l01*l02
c-
         do 200 iqu=1,ngcd
          jbuf =nbuf +(iqu-1)*incre
          jbuf2=nbuf2+(iqu-1)*incre2
          call shift0l(bl(jbuf),bl(jbuf2),
     *                 l01,l02,bl(iwij),lsmx,lsjl,
     *                 bl(ixij),nfu(nqi+1),lqij,lnkl,mnbls,
c--> *                 bl(ixij),nfu(nqi+1),nfu(nqij+1),lnkl,nbls,
     *                 bl(ixabn),bl(ixcdn),
     *                 bl(indxx),lni,lnj,lnk,lnl)
  200    continue
          call retmem(3)
          return
       endif
c
c- 1 l-shell
       if (lshellt.eq.1) then
c---
          jbuf  = nbuf 
          jbuf2 =nbuf2 
          jbfijx=ibfij1 
          if(lshelij.eq.2) jbfijx=ibfij2
          jbfklx=ibfkl1 
          if(lshelkl.eq.2) jbfklx=ibfkl2
c---
          call shift1l(bl(jbuf),bl(jbuf2),l01,l02,
     *                 bl(jbfijx),bl(jbfklx),
     *                 lqij,lqkl,
     *                 bl(iwij),lsmx, bl(ivij),lsjl,
     *                 bl(ixij),nfu(nqi+1),lnkl,bl(iyij),nfu(nqj+1),
     *                 mnbls,
     *                 bl(ixabn),bl(ixcdn),
     *                 bl(indxx),lni,lnj,lnk,lnl)
c
          call retmem(3)
          return
       endif
ccc
c- 2 l-shell
       if (lshellt.eq.2) then
c-
          jbuf  = nbuf 
          jbuf2 =nbuf2  
          jbfij1=ibfij1 
          jbfij2=ibfij2 
          jbfkl1=ibfkl1 
          jbfkl2=ibfkl2 
c-
          jbfij3=ibfij3 
          jbfkl3=ibfkl3 
          jbf2l12=ibf2l1 
          if(lcas2(2).eq.1) jbf2l12=ibf2l2 
          jbf2l34=ibf2l3 
          if(lcas2(4).eq.1) jbf2l34=ibf2l4 
c---
          call shift2l(bl(jbuf),bl(jbuf2),l01,l02,
     *                 bl(jbfij1),bl(jbfij2),bl(jbfkl1),bl(jbfkl2),
     *                 lqij,lqkl,
     *                 bl(jbfij3),bl(jbfkl3),
     *                 bl(jbf2l12),bl(jbf2l34),
     *                 bl(iwij),lsmx, bl(ivij),bl(iuij),bl(isij),lsjl,
     *          bl(ixij),nfu(nqi+1),lnkl,bl(iyij),bl(izij),nfu(nqj+1),
     *                 bl(ix2l1),mnbls,
     *                 bl(ixabn),bl(ixcdn),
     *                 bl(indxx),lni,lnj,lnk,lnl)
c-
          call retmem(3)
          return
       endif
ccc
c- 3 l-shell
c
       if (lshellt.eq.3) then
c---
          jbuf  = nbuf  
          jbuf2 =nbuf2  
          jbfij1=ibfij1 
          jbfij2=ibfij2
          jbfkl1=ibfkl1 
          jbfkl2=ibfkl2 
c-
          jbfij3=ibfij3 
          jbfkl3=ibfkl3 
          jbf2l1=ibf2l1 
          jbf2l2=ibf2l2 
          jbf2l3=ibf2l3
          jbf2l4=ibf2l4
c-
          jbf3l12=ibf3l1 
          ix3l12=ix3l1
          if(lcas3(2).eq.1) then
              jbf3l12=ibf3l2 
              ix3l12=ix3l2
          endif
          jbf3l34=ibf3l3 
          ix3l34=ix3l3
          if(lcas3(4).eq.1) then
              jbf3l34=ibf3l4
              ix3l34=ix3l4
          endif
c---
          call shift3l(bl(jbuf),bl(jbuf2),l01,l02,
     *                 bl(jbfij1),bl(jbfij2),bl(jbfkl1),bl(jbfkl2),
     *                 lqij,lqkl,
     *                 bl(jbfij3),bl(jbfkl3),
     *                 bl(jbf2l1),bl(jbf2l2),bl(jbf2l3),bl(jbf2l4),
     *                 bl(jbf3l12),bl(jbf3l34),lqmx,
     *                 bl(iwij),lsmx, bl(ivij),bl(iuij),bl(isij),lsjl,
     *           bl(ixij),nfu(nqi+1),lnkl,bl(iyij),bl(izij),nfu(nqj+1),
     *                 bl(ix2l1),bl(ix2l2),bl(ix2l3),bl(ix2l4),
     *                 bl(ix3l12),bl(ix3l34),mnbls,
     *                 bl(ixabn),bl(ixcdn),
     *                 bl(indxx),lni,lnj,lnk,lnl)
c-
          call retmem(3)
          return
       endif
cc
c- 4 l-shell
       if (lshellt.eq.4) then
c-
          jbuf  = nbuf  
          jbuf2 =nbuf2 
          jbfij1=ibfij1
          jbfij2=ibfij2 
          jbfkl1=ibfkl1 
          jbfkl2=ibfkl2 
c-
          jbfij3=ibfij3 
          jbfkl3=ibfkl3 
          jbf2l1=ibf2l1 
          jbf2l2=ibf2l2 
          jbf2l3=ibf2l3 
          jbf2l4=ibf2l4 
c
          jbf3l1=ibf3l1 
          jbf3l2=ibf3l2 
          jbf3l3=ibf3l3 
          jbf3l4=ibf3l4 
c-
          jssss =issss  
c---
          call shift4l(bl(jbuf),bl(jbuf2),l01,l02,
     *                 bl(jbfij1),bl(jbfij2),bl(jbfkl1),bl(jbfkl2),
     *                 lqij,lqkl,
     *                 bl(jbfij3),bl(jbfkl3),
     *                 bl(jbf2l1),bl(jbf2l2),bl(jbf2l3),bl(jbf2l4),
     *                 bl(jbf3l1),bl(jbf3l2),bl(jbf3l3),bl(jbf3l4),lqmx,
     *                 bl(jssss),
     *                 bl(iwij),lsmx, bl(ivij),bl(iuij),bl(isij),lsjl,
     *           bl(ixij),nfu(nqi+1),lnkl,bl(iyij),bl(izij),nfu(nqj+1),
     *                 bl(ix2l1),bl(ix2l2),bl(ix2l3),bl(ix2l4),
     *                 bl(ix3l1),bl(ix3l2),bl(ix3l3),bl(ix3l4),mnbls,
     *                 bl(ixabn),bl(ixcdn),
     *                 bl(indxx),lni,lnj,lnk,lnl)
c-
          call retmem(3)
          return
       endif
c----------------------
      return
      end
c=====================================================
      subroutine convr3(bl,m,nbls,npij,npkl,idx1,idx2,
     *                   xab,xcd, ixabn,ixcdn)
      implicit real*8 (a-h,o-z)
cccc  common /big/ bl(1)
      dimension bl(*)
      dimension idx1(*),idx2(*)
      dimension xab(npij,3),xcd(npkl,3)
c
      nbls1=nbls
      nbls2=nbls*2
      nbls3=nbls*3
      nbls1=nbls1*m
      nbls2=nbls2*m
      nbls3=nbls3*m
      call getmem(nbls3,ixabn)
      call getmem(nbls3,ixcdn)
c
       ixab1=ixabn-1
       ixcd1=ixcdn-1
c
      ijklnmr=0
      do 100 ijkl=1,nbls
      ijpar=idx1(ijkl)
      klpar=idx2(ijkl)
c
      xab1=xab(ijpar,1)
      xab2=xab(ijpar,2)
      xab3=xab(ijpar,3)
      xcd1=xcd(klpar,1)
      xcd2=xcd(klpar,2)
      xcd3=xcd(klpar,3)
c
        do 100 nmr=1,m
        ijklnmr=ijklnmr+1
        bl(ixab1+ijklnmr)      =xab1
        bl(ixab1+ijklnmr+nbls1)=xab2
        bl(ixab1+ijklnmr+nbls2)=xab3
c
        bl(ixcd1+ijklnmr)      =xcd1
        bl(ixcd1+ijklnmr+nbls1)=xcd2
        bl(ixcd1+ijklnmr+nbls2)=xcd3
c
  100 continue
      return
      end
c=====================================================
      subroutine daxpy3(n,a,z1,z2,z3,y1,y2,y3,x)
c------------------------------------------------
c* performs the vector operations with a stride=1
c* 
c*     Z = Y + A*X 
c*  
c*   for three values of a and three matrices Z and Y and this same X
c*
c------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension z1(n),z2(n),z3(n),y1(n),y2(n),y3(n),x(n),a(n,3)
c
      do 10 i=1,n
      z1(i)=y1(i) + a(i,1)*x(i)
      z2(i)=y2(i) + a(i,2)*x(i)
      z3(i)=y3(i) + a(i,3)*x(i)
   10 continue
c*
      end
c=====================================================
      subroutine noshift(lt1,lt2,mnbls,
     *                   buf,buf2,
     *                   bfijx,      bfklx,
     *                   bf2lx,
     *                   lt3,lt4,
     *                   indxx,ln12,ln34)
c**
c**   special cases where the shifting of angular
c**   momentum is not needed at all / (ds|ss)..(ss|sd),
c**   (xs|ys),(xs|sy),(sx|ys),(sx|sy) /
c**
      implicit real*8 (a-h,o-z)
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
cc
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
      common/shell/lshellt,lshelij,lshelkl,lhelp,lcas2(4),lcas3(4)
cc
      dimension buf2(mnbls,lt1,lt2),
c    * bfij1(mnbls,lt3,lt2),bfij2(mnbls,lt3,lt2),
c    * bfkl1(mnbls,lt1,lt4),bfkl2(mnbls,lt1,lt4),
c    * bf2l1(mnbls,lt3,lt4),bf2l2(mnbls,lt3,lt4),
c    * bf2l3(mnbls,lt3,lt4),bf2l4(mnbls,lt3,lt4)
c---
     * bfijx(mnbls,lt3,lt2),
     * bfklx(mnbls,lt1,lt4),
     * bf2lx(mnbls,lt3,lt4)
       dimension buf(mnbls,*)
       dimension indxx(ln12,ln34)
c-----------------------------------
c
       ijb1=ijbeg-1
       klb1=klbeg-1
c
         ijkl=0
         do 5031 i=1,lni
         ii=(i-1)*lnj
         do 5031 j=1,lnj
         ij=ii+j
         do 5031 k=1,lnk
         kk=(k-1)*lnl
         do 5031 l=1,lnl
         kl=kk+l
         ijkl=ijkl+1
ccccccc  indxx(ij+ijb1,kl+klb1)=ijkl
         indxx(ij     ,kl     )=ijkl
 5031    continue
c
         do 5034 ij=ijbeg,lnij
         do 5034 kl=klbeg,lnkl
         ijkl=indxx(ij-ijb1,kl-klb1)
            do 5034 i=1,mnbls
         buf(i,ijkl)=buf2(i,ij,kl)
 5034    continue
c
           if(lshelkl.eq.1 .or. lshelkl.eq.2) then
              do 5035 ij=ijbeg,lnij
              ijkl=indxx(ij-ijb1,1)
                 do 5035 i=1,mnbls
c--->         buf(i,ijkl)=bfkl1(i,ij,1)
              buf(i,ijkl)=bfklx(i,ij,1)
 5035         continue
           endif
c-------
           if(lshelij.eq.1 .or. lshelij.eq.2) then
              do 6035 kl=klbeg,lnkl
              ijkl=indxx(1,kl-klb1)
                 do 6035 i=1,mnbls
c----->       buf(i,ijkl)=bfij1(i,1,kl)
              buf(i,ijkl)=bfijx(i,1,kl)
 6035         continue
           endif
c---------
           if(lshellt.eq.2) then
                ijkl=indxx(1,1)
                   do 1010 i=1,mnbls
c------>           buf(i,ijkl)=bf2l1(i,1,1)
                   buf(i,ijkl)=bf2lx(i,1,1)
 1010              continue
cc
           endif
      return
      end
c=====================================================
      subroutine shift0l(buf,buf2,lt1,lt2,
     *                   wij,lt3,lsjl,xij,lt4,lt5,lt6,mnbls,xab,xcd,
     *                   indxx,lni1,lnj1,lnk1,lnl1)
c------------------------------------
c  when l-shells are not present
c------------------------------------
      implicit real*8 (a-h,o-z)
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      common/shell/lshellt,lshelij,lshelkl,lhelp,lcas2(4),lcas3(4)
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbex,klbex
      common /logic4/ nfu(1)
c
       dimension buf2(mnbls,lt1,lt2),wij(mnbls,lt3,lsjl)
       dimension xij(mnbls,lt4,lt5,lt6)
       dimension buf(mnbls,*)
       dimension xab(mnbls,3),xcd(mnbls,3)
       dimension indxx(lni1,lnj1,lnk1,lnl1)
c------------------------------------
c
       nibeg=nfu(nqi)+1
       niend=nfu(nqi+1)
       njbeg=nfu(nqj)+1
       njend=nfu(nqj+1)
ccc
       nkbeg=nfu(nqk)+1
       nkend=nfu(nqk+1)
       nlbeg=nfu(nql)+1
       nlend=nfu(nql+1)
c
       nqix=nqi
       nqjx=nqj
c
       nqkx=nqk
       nqlx=nql
       nqklx=nqkl
c
      do 100 nkl=nfu(nqklx)+1,nfu(nskl+1)
c
           do 102 nij=nfu(nqix)+1,nfu(nsij+1)
           call tfer(buf2(1,nij,nkl),wij(1,nij,1),mnbls)
  102      continue
c
       call horiz12(wij,lt3,lsjl,xab,mnbls,nqi,nqj,nsij1)
c
       do 107 nj=njbeg,njend
       do 107 ni=nibeg,niend
          call tfer(wij(1,ni,nj),xij(1,ni,nj,nkl),mnbls)
  107  continue
  100 continue
c
c------------------------------------
c this part shifts angular momentum
c   from position 3 to position 4
c----                         -------
c
      ncount=0
      do 20 i=1,lni
      do 20 j=1,lnj
      do 20 k=1,lnk
      do 20 l=1,lnl
         ncount=ncount+1
      indxx(i,j,k,l)=ncount
   20 continue
c
c------------------------------------
c
      ixyz=0
      do 300 ni=nibeg,niend
      ixyz=ixyz+1
      jxyz=0
      do 300 nj=njbeg,njend
      jxyz=jxyz+1
c
         do 301 nkl=nfu(nqkx)+1,nfu(nskl+1)
         call tfer(xij(1,ni,nj,nkl),wij(1,nkl,1),mnbls)
  301    continue
c
       call horiz12(wij,lt3,lsjl,xcd,mnbls,nqk,nql,nskl1)
c
      kxyz=0
      do 305 nk=nkbeg,nkend
      kxyz=kxyz+1
      lxyz=0
      do 305 nl=nlbeg,nlend
      lxyz=lxyz+1
      indx=indxx(ixyz,jxyz,kxyz,lxyz)
      call tfer(wij(1,nk,nl),buf(1,indx),mnbls)
  305 continue
  300 continue
c
      end
c=============================================================
      subroutine shift1l(buf,
     *                buf2,lt1,lt2,bijx,bklx,lt3,lt4,
     *                wij,lt5,vij,lt6,
     *                xij,lt7,lt8,yij,lt9,mnbls,xab,xcd,
     *                indxx,lni1,lnj1,lnk1,lnl1)
c***********************************************************
c*
c*  when 1 l-shell is present somewhere
c*          in position 1,2 or 3,4
c***********************************************************
      implicit real*8 (a-h,o-z)
cxxx
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
c
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbex,klbex
c
      common/shell/lshellt,lshelij,lshelkl,lhelp,lcas2(4),lcas3(4)
c
cxx
      common /logic4/ nfu(1)
cxx
c
      dimension buf2(mnbls,lt1,lt2),
     * bijx(mnbls,lt3,lt2),
     * bklx(mnbls,lt1,lt4)
      dimension wij(mnbls,lt5,lt6),vij(mnbls,lt5,lt6)
      dimension xij(mnbls,lt7,lt3,lt8),yij(mnbls,lt7,lt9,lt4)
      dimension buf(mnbls,*)
      dimension xab(mnbls,3),xcd(mnbls,3)
c
      dimension indxx(lni1,lnj1,lnk1,lnl1)
c***********************************************************
c
       nibeg=nfu(nqi)+1
       niend=nfu(nqi+1)
       njbeg=nfu(nqj)+1
       njend=nfu(nqj+1)
c
       nkbeg=nfu(nqk)+1
       nkend=nfu(nqk+1)
       nlbeg=nfu(nql)+1
       nlend=nfu(nql+1)
c
       nqix=nqi
       nqjx=nqj
       if(ityp.eq.3) then
          nibeg=1
          if(jtyp.le.3) nqix=1
       endif
       if(jtyp.eq.3) then
          njbeg=1
          if(ityp.le.3) nqjx=1
       endif
c
       nqkx=nqk
       nqlx=nql
       nqklx=nqkl
       if(ktyp.eq.3) then
          nkbeg=1
          if(ltyp.le.3) then
             nqklx=1
             nqkx=1
          endif
       endif
       if(ltyp.eq.3) then
          nlbeg=1
          if(ktyp.le.3) then
             nqklx=1
             nqlx=1
          endif
       endif
c
c*****
c 
      do 10 nkl=nfu(nqkl)+1,nfu(nskl+1)
c
       do 11 ij=nqi,nsij
       ijbeg=nfu(ij)+1
       ijend=nfu(ij+1)
           do 12 nij=ijbeg,ijend
           call tfer(buf2(1,nij,nkl),wij(1,nij,1),mnbls)
   12      continue
c
   11  continue
c
cccccccc
       call horiz12(wij,lt5,lt6,xab,mnbls,nqi,nqj,nsij1)
cccccccc
c
       do 16 nj=nfu(nqj)+1,nfu(nqj+1)
       do 16 ni=nfu(nqi)+1,nfu(nqi+1)
       call tfer(wij(1,ni,nj),xij(1,ni,nj,nkl),mnbls)
   16  continue
c
ccc  here lshelij can be eq. 0, 1 or 2 only  ccc
c
      if(lshelij.gt.0) then
          if(lshelij.eq.1) then
              if(jtyp.eq.1) then
                call tfer(bijx(1,1,nkl),xij(1,1,1,nkl),mnbls)
              else
                call daxpy3(mnbls,xab,
     *          xij(1,1,2,nkl),xij(1,1,3,nkl),xij(1,1,4,nkl),
     *          bijx(1,2,nkl),bijx(1,3,nkl),bijx(1,4,nkl),bijx(1,1,nkl))
              endif
          else
               do 17 nij=nfu(nqi )+1,nfu(nqi +1)
               call tfer(bijx(1,nij,nkl),xij(1,nij,1,nkl),mnbls)
   17          continue
          endif
      endif
c
   10 continue
c******
ccc  here lshelkl can be eq. 0, 1 or 2 only  ccc
      if(lshelkl.gt.0) then
ccc
      do 100 nkl=nfu(nqklx)+1,nfu(nqkl+1)
c
       do 101 ij=nqix,nsij
       ijbeg=nfu(ij)+1
       ijend=nfu(ij+1)
c
             do 103 nij=ijbeg,ijend
             call tfer(bklx(1,nij,nkl),vij(1,nij,1),mnbls)
  103        continue
  101  continue
c
cccccc
        call horiz12(vij,lt5,lt6,xab,mnbls,nqi,nqj,nsij1)
cccccc
            do 1071 nj=njbeg,njend
            do 1071 ni=nibeg,niend
            call tfer(vij(1,ni,nj),yij(1,ni,nj,nkl),mnbls)
 1071       continue
c
  100 continue
c
      endif
c
c***********************************************************
c*    this part    shifts the angular momentum
c*         from position 3 to position 4
c****                                                   ****
c
      ncount=0
      do 20 i=1,lni
      do 20 j=1,lnj
      do 20 k=1,lnk
      do 20 l=1,lnl
         ncount=ncount+1
      indxx(i,j,k,l)=ncount
   20 continue
c
c***********************************************************
c
      ixyz=0
      do 300 ni=nibeg,niend
      ixyz=ixyz+1
      jxyz=0
      do 300 nj=njbeg,njend
      jxyz=jxyz+1
c
         do 301 nkl=nfu(nqkx)+1,nfu(nskl+1)
         call tfer(xij(1,ni,nj,nkl),wij(1,nkl,1),mnbls)
  301    continue
c
        call horiz12(wij,lt5,lt6,xcd,mnbls,nqk,nql,nskl1)
c
      if(lshelkl.gt.0)  then
c
         if(lshelkl.eq.1) then
            if(ltyp.eq.1) then
             call tfer(yij(1,ni,nj,1),wij(1,1,1),mnbls)
            else
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     *    yij(1,ni,nj,2),yij(1,ni,nj,3),yij(1,ni,nj,4),yij(1,ni,nj,1))
            endif
         else
             do 312 nkl=nfu(nqkx)+1,nfu(nqkl+1)
             call tfer(yij(1,ni,nj,nkl),wij(1,nkl,1),mnbls)
  312        continue
         endif
      endif
c
c
ccccccccccccccccccccccccccccccccccc
      kxyz=0
      do 305 nk=nkbeg,nkend
      kxyz=kxyz+1
      lxyz=0
      do 305 nl=nlbeg,nlend
      lxyz=lxyz+1
      indx=indxx(ixyz,jxyz,kxyz,lxyz)
      call tfer(wij(1,nk,nl),buf(1,indx),mnbls)
  305 continue
  300 continue
c****
c
      return
      end
c=============================================================
      subroutine shift2l(buf,
     *                buf2,lt1,lt2,bij1,bij2,bkl1,bkl2,lt3,lt4,
     *                bij3,bkl3,b2l12,b2l34,
     *                wij,lt5,vij,uij,sij,lt6,
     *                xij,lt7,lt8,yij,zij,lt9,
     *                x2l,mnbls,xab,xcd,
     *                indxx,lni1,lnj1,lnk1,lnl1)
c***********************************************************
c*
c*  when 2 l-shells are present somewhere
c*          in position 1,2 or 3,4
c***********************************************************
      implicit real*8 (a-h,o-z)
cxxxx
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
c
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbex,klbex
c
      common/shell/lshellt,lshelij,lshelkl,lhelp,lcas2(4),lcas3(4)
c
cxx
      common /logic4/ nfu(1)
cxx
c
      dimension x2l(mnbls,lt3,lt4)
      dimension buf2(mnbls,lt1,lt2),
     * bij1(mnbls,lt3,lt2),bij2(mnbls,lt3,lt2),
     * bkl1(mnbls,lt1,lt4),bkl2(mnbls,lt1,lt4),
     * bij3(mnbls,lt2),bkl3(mnbls,lt1),
     * b2l12(mnbls,lt3,lt4),
     * b2l34(mnbls,lt3,lt4)
       dimension wij(mnbls,lt5,lt6),
     * vij(mnbls,lt5,lt6),uij(mnbls,lt5,lt6),sij(mnbls,lt5,lt6)
      dimension xij(mnbls,lt7,lt3,lt8),yij(mnbls,lt7,lt9,lt4),
     *                             zij(mnbls,lt7,lt9,lt4)
      dimension buf(mnbls,*)
      dimension xab(mnbls,3),xcd(mnbls,3)     
      dimension indxx(lni1,lnj1,lnk1,lnl1)
c***********************************************************
c
       nibeg=nfu(nqi)+1
       niend=nfu(nqi+1)
       njbeg=nfu(nqj)+1
       njend=nfu(nqj+1)
c
       nkbeg=nfu(nqk)+1
       nkend=nfu(nqk+1)
       nlbeg=nfu(nql)+1
       nlend=nfu(nql+1)
c
       nqix=nqi
       nqjx=nqj
       if(ityp.eq.3) then
          nibeg=1
          if(jtyp.le.3) nqix=1
       endif
       if(jtyp.eq.3) then
          njbeg=1
          if(ityp.le.3) nqjx=1
       endif
c
       nqkx=nqk
       nqlx=nql
       nqklx=nqkl
       if(ktyp.eq.3) then
          nkbeg=1
          if(ltyp.le.3) then
             nqklx=1
             nqkx=1
          endif
       endif
       if(ltyp.eq.3) then
          nlbeg=1
          if(ktyp.le.3) then
             nqklx=1
             nqlx=1
          endif
       endif
c
c*****
      do 100 nkl=nfu(nqklx)+1,nfu(nskl+1)
c
       do 101 ij=nqix,nsij
       ijbeg=nfu(ij)+1
       ijend=nfu(ij+1)
           do 102 nij=ijbeg,ijend
           call tfer(buf2(1,nij,nkl),wij(1,nij,1),mnbls)
  102      continue
c
       if( nkl.le.nfu(nqkl+1)) then
         if(lshelkl.eq.1.or.lshelkl.eq.3) then
           do 103 nij=ijbeg,ijend
           call tfer(bkl1(1,nij,nkl),vij(1,nij,1),mnbls)
  103      continue
         endif
         if(lshelkl.eq.2.or.lshelkl.eq.3) then
           do 104 nij=ijbeg,ijend
           call tfer(bkl2(1,nij,nkl),uij(1,nij,1),mnbls)
  104      continue
         endif
         if(lshelkl.eq.3.and.nkl.eq.1) then
           do 105 nij=ijbeg,ijend
           call tfer(bkl3(1,nij),sij(1,nij,1),mnbls)
  105      continue
         endif
       endif
c
  101  continue
c
ccccccccccccccccccccccc
c
c
       call horiz12(wij,lt5,lt6,xab,mnbls,nqi,nqj,nsij1)
c      
      if( nkl.le.nfu(nqkl+1)) then
          if(lshelkl.eq.1.or.lshelkl.eq.3) then
            call horiz12(vij,lt5,lt6,xab,mnbls,nqi,nqj,nsij1)
          endif
          if(lshelkl.eq.2.or.lshelkl.eq.3) then
            call horiz12(uij,lt5,lt6,xab,mnbls,nqi,nqj,nsij1)
          endif
      endif
      if(nkl.eq.1) then
          if(lshelkl.eq.3) then
            call horiz12(sij,lt5,lt6,xab,mnbls,nqi,nqj,nsij1)
          endif
      endif
cccc
c
       do 107 nj=njbeg,njend
       do 107 ni=nibeg,niend
       call tfer(wij(1,ni,nj),xij(1,ni,nj,nkl),mnbls)
  107  continue
       if( nkl.le.nfu(nqkl+1)) then
           if(lshelkl.eq.1.or.lshelkl.eq.3) then
                do 1071 nj=njbeg,njend
                do 1071 ni=nibeg,niend
                call tfer(vij(1,ni,nj),yij(1,ni,nj,nkl),mnbls)
 1071           continue
           endif
           if(lshelkl.eq.2.or.lshelkl.eq.3) then
                do 1072 nj=njbeg,njend
                do 1072 ni=nibeg,niend
                call tfer(uij(1,ni,nj),zij(1,ni,nj,nkl),mnbls)
 1072           continue
           endif
       endif
c
       if(lshelij.eq.1.or.lshelij.eq.3) then
          if(jtyp.eq.1) then
             call tfer(bij1(1,1,nkl),xij(1,1,1,nkl),mnbls)
          else
      call daxpy3(mnbls,xab,
     *          xij(1,1,2,nkl),xij(1,1,3,nkl),xij(1,1,4,nkl),
     *          bij1(1,2,nkl),bij1(1,3,nkl),bij1(1,4,nkl),bij1(1,1,nkl))
          endif
       endif
       if(lshelij.eq.2.or.lshelij.eq.3) then
           do 108 nij=nfu(nqi )+1,nfu(nqi +1)
           call tfer(bij2(1,nij,nkl),xij(1,nij,1,nkl),mnbls)
  108      continue
       endif
       if(lshelij.eq.3) then
           call tfer(bij3(1,nkl),xij(1,1,1,nkl),mnbls)
       endif
c*****
      if( nkl.le.nfu(nqkl+1)) then
c****  2 l-shells ****
c
           if(lshelij.eq.1) then
              if(jtyp.eq.1) then
                if(lcas2(1).eq.1 .or. lcas2(2).eq.1) then
                   call tfer(b2l12(1,1,nkl),x2l(1,1,nkl),mnbls)
                endif
              else
                if(lcas2(1).eq.1 .or. lcas2(2).eq.1) then
                   call daxpy3(mnbls,xab,
     *             x2l(1,2,nkl),x2l(1,3,nkl),x2l(1,4,nkl),
     *             b2l12(1,2,nkl),b2l12(1,3,nkl),b2l12(1,4,nkl),
     *             b2l12(1,1,nkl))
                endif
              endif
           endif
           if(lshelij.eq.2) then
                if(lcas2(3).eq.1 .or. lcas2(4).eq.1) then
                  do 109 nij=nfu(nqix)+1,nfu(nqij+1)
                  call tfer(b2l34(1,nij,nkl),x2l(1,nij,nkl),mnbls)
  109             continue
                endif
           endif
c
      endif
c
  100 continue
c
c***********************************************************
c*    this part    shifts the angular momentum
c*         from position 3 to position 4
c****                                                   ****
c
      ncount=0
      do 20 i=1,lni
      do 20 j=1,lnj
      do 20 k=1,lnk
      do 20 l=1,lnl
         ncount=ncount+1
      indxx(i,j,k,l)=ncount
   20 continue
c
c***********************************************************
c
      ixyz=0
      do 300 ni=nibeg,niend
      ixyz=ixyz+1
      jxyz=0
      do 300 nj=njbeg,njend
      jxyz=jxyz+1
c
         do 301 nkl=nfu(nqkx)+1,nfu(nskl+1)
         call tfer(xij(1,ni,nj,nkl),wij(1,nkl,1),mnbls)
  301    continue
c
cccccccccc
c
       call horiz12(wij,lt5,lt6,xcd,mnbls,nqk,nql,nskl1)
cccccccccc
c
       if(lshelkl.eq.1.or.lshelkl.eq.3) then
         if(ltyp.eq.1) then
           call tfer(yij(1,ni,nj,1),wij(1,1,1),mnbls)
         else
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     * yij(1,ni,nj,2),yij(1,ni,nj,3),yij(1,ni,nj,4),yij(1,ni,nj,1))
         endif
       endif
       if(lshelkl.eq.2.or.lshelkl.eq.3) then
         do 312 nkl=nfu(nqkx)+1,nfu(nqkl+1)
         call tfer(zij(1,ni,nj,nkl),wij(1,nkl,1),mnbls)
  312    continue
       endif
c
       if(lshelkl.eq.3) then
         call tfer(sij(1,ni,nj),wij(1,1,1),mnbls)
       endif
c
       if( ni.le.nfu(nqi +1).and.nj.le.nfu(nqj +1) ) then
c
c****  2 l-shells ****
         if(lshelkl.eq.1) then
            if(ltyp.eq.1) then
              if(lcas2(1).eq.1) then
               if(ni.eq.1.and.nj.ge.nfu(nqj)+1) then 
                   call tfer(x2l(1,nj,1),wij(1,1,1),mnbls)
               endif
              endif
              if(lcas2(3).eq.1) then
               if(nj.eq.1.and.ni.ge.nfu(nqi)+1) then
                   call tfer(x2l(1,ni,1),wij(1,1,1),mnbls)
               endif
              endif
            else
              if(lcas2(1).eq.1.and.ni.eq.1.and.nj.ge.nfu(nqj)+1) then
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     *            x2l(1,nj,2),x2l(1,nj,3),x2l(1,nj,4),x2l(1,nj,1))
              endif
              if(lcas2(3).eq.1.and.nj.eq.1.and.ni.ge.nfu(nqi)+1) then
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     *            x2l(1,ni,2),x2l(1,ni,3),x2l(1,ni,4),x2l(1,ni,1))
              endif
            endif
         endif
         if(lshelkl.eq.2) then
              do 313 nkl=nfu(nqkx)+1,nfu(nqkl+1)
              if(lcas2(2).eq.1.and.ni.eq.1) then 
                  call tfer(x2l(1,nj,nkl),wij(1,nkl,1),mnbls)
              endif
              if(lcas2(4).eq.1.and.nj.eq.1) then
                  call tfer(x2l(1,ni,nkl),wij(1,nkl,1),mnbls)
              endif
  313         continue
         endif
       endif
c****
c
ccccccccccccccccccccccccccccccccccc
      kxyz=0
      do 305 nk=nkbeg,nkend
      kxyz=kxyz+1
      lxyz=0
      do 305 nl=nlbeg,nlend
      lxyz=lxyz+1
      indx=indxx(ixyz,jxyz,kxyz,lxyz)
      call tfer(wij(1,nk,nl),buf(1,indx),mnbls)
  305 continue
  300 continue
c
      return
      end
c=============================================================
      subroutine shift3l(buf,
     *                buf2,lt1,lt2,bij1,bij2,bkl1,bkl2,lt3,lt4,
     *                bij3,bkl3,b2l1,b2l2,b2l3,b2l4,
     *                b3l12,b3l34,lt5,
     *                wij,lt6,vij,uij,sij,lt7,
     *                xij,lt8,lt9,yij,zij,lt10,
     *                x2l1,x2l2,x2l3,x2l4,x3l12,x3l34,mnbls,
     *                xab,xcd,
     *                indxx,lni1,lnj1,lnk1,lnl1)
c***********************************************************
c*
c*  when 3 l-shells are present somewhere
c*          in positions 1,2 or 3,4
c***********************************************************
      implicit real*8 (a-h,o-z)
cxxx
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
c
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbex,klbex
c
      common/shell/lshellt,lshelij,lshelkl,lhelp,lcas2(4),lcas3(4)
c
cxx
      common /logic4/ nfu(1)
cxx
c
      dimension x2l1(mnbls,lt3,lt4),x2l2(mnbls,lt3,lt4),
     *          x2l3(mnbls,lt3,lt4),x2l4(mnbls,lt3,lt4)
      dimension x3l12(mnbls,lt4),x3l34(mnbls,lt3)
c
      dimension buf2(mnbls,lt1,lt2),
     * bij1(mnbls,lt3,lt2),bij2(mnbls,lt3,lt2),
     * bkl1(mnbls,lt1,lt4),bkl2(mnbls,lt1,lt4),
     * bij3(mnbls,lt2),bkl3(mnbls,lt1),
     * b2l1(mnbls,lt3,lt4),b2l2(mnbls,lt3,lt4),
     * b2l3(mnbls,lt3,lt4),b2l4(mnbls,lt3,lt4),
     * b3l12(mnbls,lt5),b3l34(mnbls,lt5)
      dimension wij(mnbls,lt6,lt7),
     * vij(mnbls,lt6,lt7),uij(mnbls,lt6,lt7),sij(mnbls,lt6,lt7)
      dimension xij(mnbls,lt8,lt3,lt9),yij(mnbls,lt8,lt10,lt4)
     *                            ,zij(mnbls,lt8,lt10,lt4)
      dimension buf(mnbls,*)
      dimension xab(mnbls,3),xcd(mnbls,3)
      dimension indxx(lni1,lnj1,lnk1,lnl1)
c***********************************************************
c
       nibeg=nfu(nqi)+1
       niend=nfu(nqi+1)
       njbeg=nfu(nqj)+1
       njend=nfu(nqj+1)
c
       nkbeg=nfu(nqk)+1
       nkend=nfu(nqk+1)
       nlbeg=nfu(nql)+1
       nlend=nfu(nql+1)
c
       nqix=nqi
       nqjx=nqj
       if(ityp.eq.3) then
          nibeg=1
          if(jtyp.le.3) nqix=1
       endif
       if(jtyp.eq.3) then
          njbeg=1
          if(ityp.le.3) nqjx=1
       endif
c
       nqkx=nqk
       nqlx=nql
       nqklx=nqkl
       if(ktyp.eq.3) then
          nkbeg=1
          if(ltyp.le.3) then
             nqklx=1
             nqkx=1
          endif
       endif
       if(ltyp.eq.3) then
          nlbeg=1
          if(ktyp.le.3) then
             nqklx=1
             nqlx=1
          endif
       endif
c
c*****
c
      do 100 nkl=nfu(nqklx)+1,nfu(nskl+1)
c
       do 101 ij=nqix,nsij
       ijbeg=nfu(ij)+1
       ijend=nfu(ij+1)
           do 102 nij=ijbeg,ijend
           call tfer(buf2(1,nij,nkl),wij(1,nij,1),mnbls)
  102      continue
c
       if( nkl.le.nfu(nqkl+1)) then
         if(lshelkl.eq.1.or.lshelkl.eq.3) then
           do 103 nij=ijbeg,ijend
           call tfer(bkl1(1,nij,nkl),vij(1,nij,1),mnbls)
  103      continue
         endif
         if(lshelkl.eq.2.or.lshelkl.eq.3) then
           do 104 nij=ijbeg,ijend
           call tfer(bkl2(1,nij,nkl),uij(1,nij,1),mnbls)
  104      continue
         endif
         if(lshelkl.eq.3.and.nkl.eq.1) then
           do 105 nij=ijbeg,ijend
           call tfer(bkl3(1,nij),sij(1,nij,1),mnbls)
  105      continue
         endif
       endif
c
  101  continue
c
cccccccccccc
c
       call horiz12(wij,lt6,lt7,xab,mnbls,nqi,nqj,nsij1)
c
      if( nkl.le.nfu(nqkl+1)) then
          if(lshelkl.eq.1.or.lshelkl.eq.3) then
            call horiz12(vij,lt6,lt7,xab,mnbls,nqi,nqj,nsij1)
          endif
          if(lshelkl.eq.2.or.lshelkl.eq.3) then
            call horiz12(uij,lt6,lt7,xab,mnbls,nqi,nqj,nsij1)
          endif
      endif
      if(nkl.eq.1) then
          if(lshelkl.eq.3) then
            call horiz12(sij,lt6,lt7,xab,mnbls,nqi,nqj,nsij1)
          endif
      endif
cccccccccccc
c
       do 107 nj=njbeg,njend
       do 107 ni=nibeg,niend
       call tfer(wij(1,ni,nj),xij(1,ni,nj,nkl),mnbls)
  107  continue
       if( nkl.le.nfu(nqkl+1)) then
           if(lshelkl.eq.1.or.lshelkl.eq.3) then
                do 1071 nj=njbeg,njend
                do 1071 ni=nibeg,niend
                call tfer(vij(1,ni,nj),yij(1,ni,nj,nkl),mnbls)
 1071           continue
           endif
           if(lshelkl.eq.2.or.lshelkl.eq.3) then
                do 1072 nj=njbeg,njend
                do 1072 ni=nibeg,niend
                call tfer(uij(1,ni,nj),zij(1,ni,nj,nkl),mnbls)
 1072           continue
           endif
       endif
c
       if(lshelij.eq.1.or.lshelij.eq.3) then
          if(jtyp.eq.1) then
             call tfer(bij1(1,1,nkl),xij(1,1,1,nkl),mnbls)
          else
      call daxpy3(mnbls,xab,
     *    xij(1,1,2,nkl),xij(1,1,3,nkl),xij(1,1,4,nkl),
     *    bij1(1,2,nkl),bij1(1,3,nkl),bij1(1,4,nkl),bij1(1,1,nkl))
          endif
       endif
       if(lshelij.eq.2.or.lshelij.eq.3) then
           do 108 nij=nfu(nqi )+1,nfu(nqi +1)
           call tfer(bij2(1,nij,nkl),xij(1,nij,1,nkl),mnbls)
  108      continue
       endif
       if(lshelij.eq.3) then
           call tfer(bij3(1,nkl),xij(1,1,1,nkl),mnbls)
       endif
c*****
      if( nkl.le.nfu(nqkl+1)) then
c****  2 or 3 l-shells ****
c
           if(lshelij.eq.1) then
              if(jtyp.eq.1) then
                   call tfer(b2l1(1,1,nkl),x2l1(1,1,nkl),mnbls)
                   call tfer(b2l2(1,1,nkl),x2l2(1,1,nkl),mnbls)
c-- test ?      if(nkl.eq.1) then
                if(nkl.eq.1 .and. lcas3(3).eq.1) then
                   call tfer(b3l34(1,1),x3l34(1,1),mnbls)
                endif
              else
                call daxpy3(mnbls,xab,
     *          x2l1(1,2,nkl),x2l1(1,3,nkl),x2l1(1,4,nkl),
     *          b2l1(1,2,nkl),b2l1(1,3,nkl),b2l1(1,4,nkl),b2l1(1,1,nkl))
cc
                call daxpy3(mnbls,xab,
     *          x2l2(1,2,nkl),x2l2(1,3,nkl),x2l2(1,4,nkl),
     *          b2l2(1,2,nkl),b2l2(1,3,nkl),b2l2(1,4,nkl),b2l2(1,1,nkl))
c-- test ?      if(nkl.eq.1) then
                if(nkl.eq.1 .and. lcas3(3).eq.1) then
                  call daxpy3(mnbls,xab,
     *            x3l34(1,2),x3l34(1,3),x3l34(1,4),
     *            b3l34(1,2),b3l34(1,3),b3l34(1,4),b3l34(1,1))
                endif
              endif
           endif
           if(lshelij.eq.2) then
                do 109 nij=nfu(nqix)+1,nfu(nqij+1)
                   call tfer(b2l3(1,nij,nkl),x2l3(1,nij,nkl),mnbls)
                   call tfer(b2l4(1,nij,nkl),x2l4(1,nij,nkl),mnbls)
  109           continue
c-- test ?      if(nkl.eq.1) then
                if(nkl.eq.1 .and. lcas3(4).eq.1) then
                   do 110 nij=nfu(nqi )+1,nfu(nqi +1)
                   call tfer(b3l34(1,nij),x3l34(1,nij),mnbls)
  110              continue
                endif
           endif
           if(lshelij.eq.3) then
                if(lcas2(1).eq.1) then
                  call daxpy3(mnbls,xab,
     *            x2l1(1,2,nkl),x2l1(1,3,nkl),x2l1(1,4,nkl),
     *            b2l1(1,2,nkl),b2l1(1,3,nkl),b2l1(1,4,nkl),
     *            b2l1(1,1,nkl))
cc
                  call tfer(b2l3(1,2,nkl),x2l3(1,2,nkl),mnbls)
                  call tfer(b2l3(1,3,nkl),x2l3(1,3,nkl),mnbls)
                  call tfer(b2l3(1,4,nkl),x2l3(1,4,nkl),mnbls)
                endif
                if(lcas2(2).eq.1) then
                  call daxpy3(mnbls,xab,
     *            x2l2(1,2,nkl),x2l2(1,3,nkl),x2l2(1,4,nkl),
     *            b2l2(1,2,nkl),b2l2(1,3,nkl),b2l2(1,4,nkl),
     *            b2l2(1,1,nkl))
cc
                  call tfer(b2l4(1,2,nkl),x2l4(1,2,nkl),mnbls)
                  call tfer(b2l4(1,3,nkl),x2l4(1,3,nkl),mnbls)
                  call tfer(b2l4(1,4,nkl),x2l4(1,4,nkl),mnbls)
                endif
                if(lcas3(1).eq.1) then
                  call tfer(b3l12(1,nkl),x3l12(1,nkl),mnbls)
                endif
                if(lcas3(2).eq.1) then
                  call tfer(b3l12(1,nkl),x3l12(1,nkl),mnbls)
                endif
           endif
      endif
c
  100 continue
c
c***********************************************************
c*    this part    shifts the angular momentum
c*         from position 3 to position 4
c****                                                   ****
c
      ncount=0
      do 20 i=1,lni
      do 20 j=1,lnj
      do 20 k=1,lnk
      do 20 l=1,lnl
         ncount=ncount+1
      indxx(i,j,k,l)=ncount
   20 continue
c
c***********************************************************
c
      ixyz=0
      do 300 ni=nibeg,niend
      ixyz=ixyz+1
      jxyz=0
      do 300 nj=njbeg,njend
      jxyz=jxyz+1
c
         do 301 nkl=nfu(nqkx)+1,nfu(nskl+1)
         call tfer(xij(1,ni,nj,nkl),wij(1,nkl,1),mnbls)
  301    continue
c
ccccccc
        call horiz12(wij,lt6,lt7,xcd,mnbls,nqk,nql,nskl1)
ccccccc
cc
       if(lshelkl.eq.1.or.lshelkl.eq.3) then
         if(ltyp.eq.1) then
           call tfer(yij(1,ni,nj,1),wij(1,1,1),mnbls)
         else
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     *   yij(1,ni,nj,2),yij(1,ni,nj,3),yij(1,ni,nj,4),yij(1,ni,nj,1))
         endif
       endif
       if(lshelkl.eq.2.or.lshelkl.eq.3) then
         do 312 nkl=nfu(nqkx)+1,nfu(nqkl+1)
         call tfer(zij(1,ni,nj,nkl),wij(1,nkl,1),mnbls)
  312    continue
       endif
c
       if(lshelkl.eq.3) then
         call tfer(sij(1,ni,nj),wij(1,1,1),mnbls)
       endif
c
       if( ni.le.nfu(nqi +1).and.nj.le.nfu(nqj +1) ) then
c
c****  2,3  l-shells ****
         if(lshelkl.eq.1) then
            if(ltyp.eq.1) then
               if(ni.eq.1.and.nj.ge.nfu(nqj)+1) then
                  call tfer(x2l1(1,nj,1),wij(1,1,1),mnbls)
               endif
               if(nj.eq.1.and.ni.ge.nfu(nqi)+1) then
                  call tfer(x2l3(1,ni,1),wij(1,1,1),mnbls)
               endif
c--test ?      if(ni.eq.1.and.nj.eq.1) then
               if(ni.eq.1.and.nj.eq.1 .and. lcas3(1).eq.1) then
                  call tfer(x3l12(1,1),wij(1,1,1),mnbls)
               endif
            else
              if(ni.eq.1.and.nj.ge.nfu(nqj)+1) then
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     *  x2l1(1,nj,2),x2l1(1,nj,3),x2l1(1,nj,4),x2l1(1,nj,1))
              endif
              if(nj.eq.1.and.ni.ge.nfu(nqi)+1) then
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     *  x2l3(1,ni,2),x2l3(1,ni,3),x2l3(1,ni,4),x2l3(1,ni,1))
              endif
c--test ?     if(ni.eq.1.and.nj.eq.1) then
              if(ni.eq.1.and.nj.eq.1 .and. lcas3(1).eq.1) then
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     *               x3l12(1,2),x3l12(1,3),x3l12(1,4),x3l12(1,1))
              endif
            endif
c
         endif
         if(lshelkl.eq.2) then
              do 313 nkl=nfu(nqkx)+1,nfu(nqkl+1)
              if(ni.eq.1)  then
                  call tfer(x2l2(1,nj,nkl),wij(1,nkl,1),mnbls)
              endif
              if(nj.eq.1)  then
                  call tfer(x2l4(1,ni,nkl),wij(1,nkl,1),mnbls)
              endif
c--test ?     if(ni.eq.1.and.nj.eq.1) then
              if(ni.eq.1.and.nj.eq.1 .and. lcas3(2).eq.1) then
                  call tfer(x3l12(1,nkl),wij(1,nkl,1),mnbls)
              endif
  313         continue
         endif
c
c****  3 l-shells ****
         if(lshelkl.eq.3) then
              if(lcas2(1).eq.1.and.ni.eq.1) then
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     *  x2l1(1,nj,2),x2l1(1,nj,3),x2l1(1,nj,4),x2l1(1,nj,1))
cc
              call tfer(x2l2(1,nj,2),wij(1,2,1),mnbls)
              call tfer(x2l2(1,nj,3),wij(1,3,1),mnbls)
              call tfer(x2l2(1,nj,4),wij(1,4,1),mnbls)
              endif
              if(lcas2(3).eq.1.and.nj.eq.1) then
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     *  x2l3(1,ni,2),x2l3(1,ni,3),x2l3(1,ni,4),x2l3(1,ni,1))
cc
              call tfer(x2l4(1,ni,2),wij(1,2,1),mnbls)
              call tfer(x2l4(1,ni,3),wij(1,3,1),mnbls)
              call tfer(x2l4(1,ni,4),wij(1,4,1),mnbls)
              endif
c
              if(lcas3(3).eq.1.and.ni.eq.1) then
                  call tfer(x3l34(1,nj),wij(1,1,1),mnbls)         
              endif
              if(lcas3(4).eq.1.and.nj.eq.1) then
                  call tfer(x3l34(1,ni),wij(1,1,1),mnbls)
              endif
         endif
       endif
c
ccccccccccccccccccccccccccccccccccc
      kxyz=0
      do 305 nk=nkbeg,nkend
      kxyz=kxyz+1
      lxyz=0
      do 305 nl=nlbeg,nlend
      lxyz=lxyz+1
      indx=indxx(ixyz,jxyz,kxyz,lxyz)
      call tfer(wij(1,nk,nl),buf(1,indx),mnbls)
  305 continue
  300 continue
c
      return
      end
c=============================================================
      subroutine shift4l(buf,
     *                buf2,lt1,lt2,bij1,bij2,bkl1,bkl2,lt3,lt4,
     *                bij3,bkl3,b2l1,b2l2,b2l3,b2l4,
     *                b3l1,b3l2,b3l3,b3l4,lt5,ssss,
     *                wij,lt6,vij,uij,sij,lt7,
     *                xij,lt8,lt9,yij,zij,lt10,
     *                x2l1,x2l2,x2l3,x2l4,x3l1,x3l2,x3l3,x3l4,mnbls,
     *                xab,xcd,
     *                indxx,lni1,lnj1,lnk1,lnl1)
c***********************************************************
c*        when 4 l-shells are present
c***********************************************************
      implicit real*8 (a-h,o-z)
ctest
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
ctest
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbex,klbex
      common/shell/lshellt,lshelij,lshelkl,lhelp,lcas2(4),lcas3(4)
      common /logic4/ nfu(1)
c
      dimension x2l1(mnbls,lt3,lt4),x2l2(mnbls,lt3,lt4),
     *          x2l3(mnbls,lt3,lt4),x2l4(mnbls,lt3,lt4)
      dimension x3l1(mnbls,lt4),x3l2(mnbls,lt4),
     *          x3l3(mnbls,lt3),x3l4(mnbls,lt3)
c
      dimension buf2(mnbls,lt1,lt2),
     * bij1(mnbls,lt3,lt2),bij2(mnbls,lt3,lt2),
     * bkl1(mnbls,lt1,lt4),bkl2(mnbls,lt1,lt4),
     * bij3(mnbls,lt2),bkl3(mnbls,lt1),
     * b2l1(mnbls,lt3,lt4),b2l2(mnbls,lt3,lt4),
     * b2l3(mnbls,lt3,lt4),b2l4(mnbls,lt3,lt4),
cc-> * b3l(mnbls,lt5,4)
     * b3l1(mnbls,lt5),b3l2(mnbls,lt5),b3l3(mnbls,lt5),b3l4(mnbls,lt5)
      dimension wij(mnbls,lt6,lt7),
     * vij(mnbls,lt6,lt7),uij(mnbls,lt6,lt7),sij(mnbls,lt6,lt7)
      dimension xij(mnbls,lt8,lt3,lt9),yij(mnbls,lt8,lt10,lt4)
     *                            ,zij(mnbls,lt8,lt10,lt4)
c???? dimension ssss(nbls)
      dimension ssss(mnbls)
      dimension buf(mnbls,*)
      dimension xab(mnbls,3),xcd(mnbls,3)
      dimension indxx(lni1,lnj1,lnk1,lnl1)
c***********************************************************
c
       niend=nfu(nqi+1)
       njend=nfu(nqj+1)
c
       nkend=nfu(nqk+1)
       nlend=nfu(nql+1)
c
       nqix=1
       nqjx=1
       nibeg=1
       njbeg=1
c
       nqkx=1
       nqlx=1
       nqklx=1
       nkbeg=1
       nlbeg=1
c
c*****
c     do 100 nkl=nfu(nqklx)+1,nfu(nskl+1)
      do 100 nkl=1,10
c
c      do 101 ij=nqix,nsij
       do 101 ij=1,3
       ijbeg=nfu(ij)+1
       ijend=nfu(ij+1)
           do 102 nij=ijbeg,ijend
           call tfer(buf2(1,nij,nkl),wij(1,nij,1),mnbls)
  102      continue
c
       if( nkl.le.nfu(nqkl+1)) then
           do 103 nij=ijbeg,ijend
           call tfer(bkl1(1,nij,nkl),vij(1,nij,1),mnbls)
  103      continue
           do 104 nij=ijbeg,ijend
           call tfer(bkl2(1,nij,nkl),uij(1,nij,1),mnbls)
  104      continue
         if(                 nkl.eq.1) then
           do 105 nij=ijbeg,ijend
           call tfer(bkl3(1,nij),sij(1,nij,1),mnbls)
  105      continue
         endif
       endif
c
  101  continue
c
ccccccccccc
c
        call horiz12(wij,lt6,lt7,xab,mnbls,nqi,nqj,nsij1)        
c
      if( nkl.le.nfu(nqkl+1)) then
            call horiz12(vij,lt6,lt7,xab,mnbls,nqi,nqj,nsij1)
            call horiz12(uij,lt6,lt7,xab,mnbls,nqi,nqj,nsij1)
      endif
      if( nkl.eq.1) then
            call horiz12(sij,lt6,lt7,xab,mnbls,nqi,nqj,nsij1)
      endif
ccccccccccc
c
        do 107 nj=njbeg,njend
        do 107 ni=nibeg,niend
        call tfer(wij(1,ni,nj),xij(1,ni,nj,nkl),mnbls)
  107   continue
c
      if( nkl.le.nfu(nqkl+1)) then
        do 1071 nj=njbeg,njend
        do 1071 ni=nibeg,niend
        call tfer(vij(1,ni,nj),yij(1,ni,nj,nkl),mnbls)
        call tfer(uij(1,ni,nj),zij(1,ni,nj,nkl),mnbls)
 1071   continue
      endif
c
      call daxpy3(mnbls,xab,
     *   xij(1,1,2,nkl),xij(1,1,3,nkl),xij(1,1,4,nkl),
     *   bij1(1,2,nkl),bij1(1,3,nkl),bij1(1,4,nkl),bij1(1,1,nkl))
c
           do 108 nij=nfu(nqi )+1,nfu(nqi +1)
           call tfer(bij2(1,nij,nkl),xij(1,nij,1,nkl),mnbls)
  108      continue
c
           call tfer(bij3(1,nkl),xij(1,1,1,nkl),mnbls)
c
c*****
      if( nkl.le.nfu(nqkl+1)) then
      call daxpy3(mnbls,xab,x2l1(1,2,nkl),x2l1(1,3,nkl),x2l1(1,4,nkl),
     *    b2l1(1,2,nkl),b2l1(1,3,nkl),b2l1(1,4,nkl),b2l1(1,1,nkl))
                call tfer(b2l3(1,2,nkl),x2l3(1,2,nkl),mnbls)
                call tfer(b2l3(1,3,nkl),x2l3(1,3,nkl),mnbls)
                call tfer(b2l3(1,4,nkl),x2l3(1,4,nkl),mnbls)
cc
      call daxpy3(mnbls,xab,x2l2(1,2,nkl),x2l2(1,3,nkl),x2l2(1,4,nkl),
     *  b2l2(1,2,nkl),b2l2(1,3,nkl),b2l2(1,4,nkl),b2l2(1,1,nkl))
cc
                call tfer(b2l4(1,2,nkl),x2l4(1,2,nkl),mnbls)
                call tfer(b2l4(1,3,nkl),x2l4(1,3,nkl),mnbls)
                call tfer(b2l4(1,4,nkl),x2l4(1,4,nkl),mnbls)
c
                call tfer(b3l1(1,nkl),x3l1(1,nkl),mnbls)
                call tfer(b3l2(1,nkl),x3l2(1,nkl),mnbls)
c
      endif
      if(nkl.eq.1) then
c
      call daxpy3(mnbls,xab,x3l3(1,2),x3l3(1,3),x3l3(1,4),
     *  b3l3(1,2),b3l3(1,3),b3l3(1,4),b3l3(1,1))
c
                do 1101 nij=nfu(nqix)+1,nfu(nqij+1)
                call tfer(b3l4(1,nij ),x3l4(1,nij),mnbls)
 1101           continue
      endif
c
  100 continue
c
c***********************************************************
c*    this part    shifts the angular momentum
c*         from position 3 to position 4
c*                         or
c*         from position 4 to position 3
c****                                                   ****
c
      ncount=0
      do 20 i=1,lni
      do 20 j=1,lnj
      do 20 k=1,lnk
      do 20 l=1,lnl
         ncount=ncount+1
      indxx(i,j,k,l)=ncount
   20 continue
c
c***********************************************************
c
      ixyz=0
      do 300 ni=nibeg,niend
      ixyz=ixyz+1
      jxyz=0
      do 300 nj=njbeg,njend
      jxyz=jxyz+1
c
         do 301 nkl=nfu(nqkx)+1,nfu(nskl+1)
         call tfer(xij(1,ni,nj,nkl),wij(1,nkl,1),mnbls)
  301    continue
c
ccccccccc
c
           call horiz12(wij,lt6,lt7,xcd,mnbls,nqk,nql,nskl1)
c
cccccc
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     * yij(1,ni,nj,2),yij(1,ni,nj,3),yij(1,ni,nj,4),yij(1,ni,nj,1))
cc
         do 312 nkl=nfu(nqkx)+1,nfu(nqkl+1)
         call tfer(zij(1,ni,nj,nkl),wij(1,nkl,1),mnbls)
  312    continue
c
         call tfer(sij(1,ni,nj),wij(1,1,1),mnbls)
c
c****
       if( ni.le.nfu(nqi +1).and.nj.le.nfu(nqj +1) ) then
              if(ni.eq.1 .and. nj.ge.2) then
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     *  x2l1(1,nj,2),x2l1(1,nj,3),x2l1(1,nj,4),x2l1(1,nj,1))
cc
              call tfer(x2l2(1,nj,2),wij(1,2,1),mnbls)
              call tfer(x2l2(1,nj,3),wij(1,3,1),mnbls)
              call tfer(x2l2(1,nj,4),wij(1,4,1),mnbls)
cc
              endif
              if(nj.eq.1 .and. ni.ge.2) then
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     * x2l3(1,ni,2),x2l3(1,ni,3),x2l3(1,ni,4),x2l3(1,ni,1))
cc
              call tfer(x2l4(1,ni,2),wij(1,2,1),mnbls)
              call tfer(x2l4(1,ni,3),wij(1,3,1),mnbls)
              call tfer(x2l4(1,ni,4),wij(1,4,1),mnbls)
              endif
c
              if(ni.eq.1) call tfer(x3l3(1,nj),wij(1,1,1),mnbls)
              if(nj.eq.1) call tfer(x3l4(1,ni),wij(1,1,1),mnbls)
       endif
       if( ni.eq.1.and.nj.eq.1 ) then
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     *  x3l1(1,2),x3l1(1,3),x3l1(1,4),x3l1(1,1))
cc
          call tfer(x3l2(1,2),wij(1,2,1),mnbls)
          call tfer(x3l2(1,3),wij(1,3,1),mnbls)
          call tfer(x3l2(1,4),wij(1,4,1),mnbls)
cc
          call tfer(ssss(1),wij(1,1,1),mnbls)
       endif
c
ccccccccccccccccccccccccccccccccccc
      kxyz=0
      do 305 nk=nkbeg,nkend
      kxyz=kxyz+1
      lxyz=0
      do 305 nl=nlbeg,nlend
      lxyz=lxyz+1
      indx=indxx(ixyz,jxyz,kxyz,lxyz)
      call tfer(wij(1,nk,nl),buf(1,indx),mnbls)
  305 continue
  300 continue
c
      return
      end
c=====================================================
      subroutine horiz12(wij,lw1,lw2,xab,mnbls,nqi,nqj,nsij1)
      implicit real*8 (a-h,o-z)
c
      common /logic4/ nfu(1)
      common /logic5/ icoor(1)
      common /logic6/ icool(1)
      common /logic7/ ifrst(1)
      common /logic8/ ilast(1)
      common /logic11/ npxyz(3,1)
c
      dimension wij(mnbls,lw1,lw2),xab(mnbls,3)
c---------------------------------------------------
          do 110 j=2,nqj
          jbeg=nfu(j)+1
          jend=nfu(j+1)
             do 115 i=nsij1-j,nqi,-1
             ibeg=nfu(i)+1
             iend=nfu(i+1)
                do 120 nj=jbeg,jend
                njm=ilast(nj)
                kcr=icool(nj)
                   do 125 ni=ibeg,iend
                   nij=npxyz(kcr,ni)
                      do 130 n=1,mnbls
        wij(n,ni,nj)=wij(n,nij,njm)+ xab(n,kcr) *wij(n,ni,njm)
  130                 continue
  125              continue
  120          continue
  115        continue
  110     continue
c
      end
c=====================================================
