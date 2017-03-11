c $Id$
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
      subroutine amtfer (a,b,n)
      implicit none
      double precision a(*),b(*)
      integer n, i
c
c     local copy of tfer for automatic inlining by compiler
c
      do i = 1, n
         b(i) = a(i)
      enddo
c
      end
      subroutine amshift(bl,nbls,l01,l02,npij,npkl,ngcd)
      implicit real*8 (a-h,o-z)
      character*11 scftype
      character*8 where
#include "texas_lpar.fh"
      common /runtype/ scftype,where
c
      common /cpu/ intsize,iacc,icache,memreal
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
c only for first & second derivatives (for use in amshift):
cccc  common /memor4b/ibuf0
      common /memor4b/ider0,ider1,ider2
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
      call make_indxx(bl(indxx),lni,lnj,lnk,lnl)
c
c------------------------------------------------
c for odinary integrals :
c
          nbuf2=ibuf2
          nbuf=ibuf
          m=1
c
c for giao integral derivatives :
          if(where.eq.'shif') then
            nmr =ngcd*nbls*lnijkl 
            nbuf=ibuf+nmr
            m=6
          endif
c
c for gradient integral derivatives :
          if(where.eq.'forc') then
            m=9
          endif
c
c for hessian  integral derivatives :
          if(where.eq.'hess') then
            m=45
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
ccccc      if(where.ne.'forc') then
           if(where.eq.'buff' .or. where.eq.'shif') then
              call shift0l(bl(jbuf),bl(jbuf2),
     *                     l01,l02,bl(iwij),lsmx,lsjl,
     *                     bl(ixij),nfu(nqi+1),lqij,lnkl,mnbls,
     *                     bl(ixabn),bl(ixcdn),
     *                     bl(indxx),lni,lnj,lnk,lnl)
           endif
           if(where.eq.'forc') then
cccccc        jbuf0=ibuf0+(iqu-1)*nbls*l01*l02
              jbuf0=ider0+(iqu-1)*nbls*l01*l02
              iwij0=iwij +mnbls*lsmx*lsjl 
              ixij0=ixij +mnbls*nfu(nqi+1)*lqij*lnkl
              call shift0l_der1(bl(jbuf),bl(jbuf2),bl(jbuf0),
     *                         l01,l02,bl(iwij),bl(iwij0),lsmx,lsjl,
     *                         bl(ixij),bl(ixij0),nfu(nqi+1),lqij,lnkl,
     *                         mnbls,nbls,
     *                         bl(ixabn),bl(ixcdn),
     *                         bl(indxx),lni,lnj,lnk,lnl)
           endif
           if(where.eq.'hess') then
              jbuf0=ider0+(iqu-1)*nbls*l01*l02
              jbuf1=ider1+(iqu-1)*9*nbls*l01*l02
              iwij1=iwij +mnbls*lsmx*lsjl 
              ixij1=ixij +mnbls*nfu(nqi+1)*lqij*lnkl
              iwij0=iwij1 +    9*nbls*lsmx*lsjl 
              ixij0=ixij1 +    9*nbls*nfu(nqi+1)*lqij*lnkl
              call shift0l_der2(bl(jbuf),bl(jbuf2),bl(jbuf1),bl(jbuf0),
     *                         l01,l02,bl(iwij),bl(iwij1),bl(iwij0),
     *                         lsmx,lsjl,
     *                         bl(ixij),bl(ixij1),bl(ixij0),
     *                         nfu(nqi+1),lqij,lnkl, mnbls,nbls,
     *                         bl(ixabn),bl(ixcdn),
     *                         bl(indxx),lni,lnj,lnk,lnl)
           endif
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
c=======================================================
c moved into the convert.f file 
c     subroutine convr3(bl,m,nbls,npij,npkl,idx1,idx2,
c    *                   xab,xcd, ixabn,ixcdn)
c=======================================================
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
c2002
c        ijkl=0
c        do 5031 i=1,lni
c        ii=(i-1)*lnj
c        do 5031 j=1,lnj
c        ij=ii+j
c        do 5031 k=1,lnk
c        kk=(k-1)*lnl
c        do 5031 l=1,lnl
c        kl=kk+l
c        ijkl=ijkl+1
c        indxx(ij     ,kl     )=ijkl
c5031    continue
c2002
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
#include "texas_lpar.fh"
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
      common/shell/lshellt,lshelij,lshelkl,lhelp,lcas2(4),lcas3(4)
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbex,klbex
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
ccccc  nqjx=nqj
c
       nqkx=nqk
ccccc  nqlx=nql
       nqklx=nqkl
c
      do 100 nkl=nfu(nqklx)+1,nfu(nskl+1)
c
           do 102 nij=nfu(nqix)+1,nfu(nsij+1)
           call amtfer(buf2(1,nij,nkl),wij(1,nij,1),mnbls)
  102      continue
c
       call horiz12(wij,lt3,lsjl,xab,mnbls,nqi,nqj,nsij1)
c
       do 107 nj=njbeg,njend
       do 107 ni=nibeg,niend
          call amtfer(wij(1,ni,nj),xij(1,ni,nj,nkl),mnbls)
  107  continue
  100 continue
c
c------------------------------------
c this part shifts angular momentum
c   from position 3 to position 4
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
         call amtfer(xij(1,ni,nj,nkl),wij(1,nkl,1),mnbls)
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
      call amtfer(wij(1,nk,nl),buf(1,indx),mnbls)
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
#include "texas_lpar.fh"
cxxx
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbex,klbex
c
      common/shell/lshellt,lshelij,lshelkl,lhelp,lcas2(4),lcas3(4)
c
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
ccccc  nqjx=nqj
       if(ityp.eq.3) then
          nibeg=1
          if(jtyp.le.3) nqix=1
       endif
       if(jtyp.eq.3) then
          njbeg=1
ccccc     if(ityp.le.3) nqjx=1
       endif
c
       nqkx=nqk
ccccc  nqlx=nql
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
ccccc        nqlx=1
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
           call amtfer(buf2(1,nij,nkl),wij(1,nij,1),mnbls)
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
       call amtfer(wij(1,ni,nj),xij(1,ni,nj,nkl),mnbls)
   16  continue
c
ccc  here lshelij can be eq. 0, 1 or 2 only  ccc
c
      if(lshelij.gt.0) then
          if(lshelij.eq.1) then
              if(jtyp.eq.1) then
                call amtfer(bijx(1,1,nkl),xij(1,1,1,nkl),mnbls)
              else
                call daxpy3(mnbls,xab,
     *          xij(1,1,2,nkl),xij(1,1,3,nkl),xij(1,1,4,nkl),
     *          bijx(1,2,nkl),bijx(1,3,nkl),bijx(1,4,nkl),bijx(1,1,nkl))
              endif
          else
               do 17 nij=nfu(nqi )+1,nfu(nqi +1)
               call amtfer(bijx(1,nij,nkl),xij(1,nij,1,nkl),mnbls)
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
             call amtfer(bklx(1,nij,nkl),vij(1,nij,1),mnbls)
  103        continue
  101  continue
c
cccccc
        call horiz12(vij,lt5,lt6,xab,mnbls,nqi,nqj,nsij1)
cccccc
            do 1071 nj=njbeg,njend
            do 1071 ni=nibeg,niend
            call amtfer(vij(1,ni,nj),yij(1,ni,nj,nkl),mnbls)
 1071       continue
c
  100 continue
c
      endif
c
c***********************************************************
c*    this part    shifts the angular momentum
c*         from position 3 to position 4
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
         call amtfer(xij(1,ni,nj,nkl),wij(1,nkl,1),mnbls)
  301    continue
c
        call horiz12(wij,lt5,lt6,xcd,mnbls,nqk,nql,nskl1)
c
      if(lshelkl.gt.0)  then
c
         if(lshelkl.eq.1) then
            if(ltyp.eq.1) then
             call amtfer(yij(1,ni,nj,1),wij(1,1,1),mnbls)
            else
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     *    yij(1,ni,nj,2),yij(1,ni,nj,3),yij(1,ni,nj,4),yij(1,ni,nj,1))
            endif
         else
             do 312 nkl=nfu(nqkx)+1,nfu(nqkl+1)
             call amtfer(yij(1,ni,nj,nkl),wij(1,nkl,1),mnbls)
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
      call amtfer(wij(1,nk,nl),buf(1,indx),mnbls)
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
#include "texas_lpar.fh"
cxxxx
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbex,klbex
c
      common/shell/lshellt,lshelij,lshelkl,lhelp,lcas2(4),lcas3(4)
c
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
ccccc  nqjx=nqj
       if(ityp.eq.3) then
          nibeg=1
          if(jtyp.le.3) nqix=1
       endif
       if(jtyp.eq.3) then
          njbeg=1
ccccc     if(ityp.le.3) nqjx=1
       endif
c
       nqkx=nqk
ccccc  nqlx=nql
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
ccccc        nqlx=1
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
           call amtfer(buf2(1,nij,nkl),wij(1,nij,1),mnbls)
  102      continue
c
       if( nkl.le.nfu(nqkl+1)) then
         if(lshelkl.eq.1.or.lshelkl.eq.3) then
           do 103 nij=ijbeg,ijend
           call amtfer(bkl1(1,nij,nkl),vij(1,nij,1),mnbls)
  103      continue
         endif
         if(lshelkl.eq.2.or.lshelkl.eq.3) then
           do 104 nij=ijbeg,ijend
           call amtfer(bkl2(1,nij,nkl),uij(1,nij,1),mnbls)
  104      continue
         endif
         if(lshelkl.eq.3.and.nkl.eq.1) then
           do 105 nij=ijbeg,ijend
           call amtfer(bkl3(1,nij),sij(1,nij,1),mnbls)
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
       call amtfer(wij(1,ni,nj),xij(1,ni,nj,nkl),mnbls)
  107  continue
       if( nkl.le.nfu(nqkl+1)) then
           if(lshelkl.eq.1.or.lshelkl.eq.3) then
                do 1071 nj=njbeg,njend
                do 1071 ni=nibeg,niend
                call amtfer(vij(1,ni,nj),yij(1,ni,nj,nkl),mnbls)
 1071           continue
           endif
           if(lshelkl.eq.2.or.lshelkl.eq.3) then
                do 1072 nj=njbeg,njend
                do 1072 ni=nibeg,niend
                call amtfer(uij(1,ni,nj),zij(1,ni,nj,nkl),mnbls)
 1072           continue
           endif
       endif
c
       if(lshelij.eq.1.or.lshelij.eq.3) then
          if(jtyp.eq.1) then
             call amtfer(bij1(1,1,nkl),xij(1,1,1,nkl),mnbls)
          else
      call daxpy3(mnbls,xab,
     *          xij(1,1,2,nkl),xij(1,1,3,nkl),xij(1,1,4,nkl),
     *          bij1(1,2,nkl),bij1(1,3,nkl),bij1(1,4,nkl),bij1(1,1,nkl))
          endif
       endif
       if(lshelij.eq.2.or.lshelij.eq.3) then
           do 108 nij=nfu(nqi )+1,nfu(nqi +1)
           call amtfer(bij2(1,nij,nkl),xij(1,nij,1,nkl),mnbls)
  108      continue
       endif
       if(lshelij.eq.3) then
           call amtfer(bij3(1,nkl),xij(1,1,1,nkl),mnbls)
       endif
c*****
      if( nkl.le.nfu(nqkl+1)) then
c****  2 l-shells ****
c
           if(lshelij.eq.1) then
              if(jtyp.eq.1) then
                if(lcas2(1).eq.1 .or. lcas2(2).eq.1) then
                   call amtfer(b2l12(1,1,nkl),x2l(1,1,nkl),mnbls)
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
                  call amtfer(b2l34(1,nij,nkl),x2l(1,nij,nkl),mnbls)
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
         call amtfer(xij(1,ni,nj,nkl),wij(1,nkl,1),mnbls)
  301    continue
c
cccccccccc
c
       call horiz12(wij,lt5,lt6,xcd,mnbls,nqk,nql,nskl1)
cccccccccc
c
       if(lshelkl.eq.1.or.lshelkl.eq.3) then
         if(ltyp.eq.1) then
           call amtfer(yij(1,ni,nj,1),wij(1,1,1),mnbls)
         else
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     * yij(1,ni,nj,2),yij(1,ni,nj,3),yij(1,ni,nj,4),yij(1,ni,nj,1))
         endif
       endif
       if(lshelkl.eq.2.or.lshelkl.eq.3) then
         do 312 nkl=nfu(nqkx)+1,nfu(nqkl+1)
         call amtfer(zij(1,ni,nj,nkl),wij(1,nkl,1),mnbls)
  312    continue
       endif
c
       if(lshelkl.eq.3) then
         call amtfer(sij(1,ni,nj),wij(1,1,1),mnbls)
       endif
c
       if( ni.le.nfu(nqi +1).and.nj.le.nfu(nqj +1) ) then
c
c****  2 l-shells ****
         if(lshelkl.eq.1) then
            if(ltyp.eq.1) then
              if(lcas2(1).eq.1) then
               if(ni.eq.1.and.nj.ge.nfu(nqj)+1) then 
                   call amtfer(x2l(1,nj,1),wij(1,1,1),mnbls)
               endif
              endif
              if(lcas2(3).eq.1) then
               if(nj.eq.1.and.ni.ge.nfu(nqi)+1) then
                   call amtfer(x2l(1,ni,1),wij(1,1,1),mnbls)
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
                  call amtfer(x2l(1,nj,nkl),wij(1,nkl,1),mnbls)
              endif
              if(lcas2(4).eq.1.and.nj.eq.1) then
                  call amtfer(x2l(1,ni,nkl),wij(1,nkl,1),mnbls)
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
      call amtfer(wij(1,nk,nl),buf(1,indx),mnbls)
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
#include "texas_lpar.fh"
cxxx
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbex,klbex
c
      common/shell/lshellt,lshelij,lshelkl,lhelp,lcas2(4),lcas3(4)
c
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
ccccc  nqjx=nqj
       if(ityp.eq.3) then
          nibeg=1
          if(jtyp.le.3) nqix=1
       endif
       if(jtyp.eq.3) then
          njbeg=1
ccccc     if(ityp.le.3) nqjx=1
       endif
c
       nqkx=nqk
ccccc  nqlx=nql
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
ccccc        nqlx=1
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
           call amtfer(buf2(1,nij,nkl),wij(1,nij,1),mnbls)
  102      continue
c
       if( nkl.le.nfu(nqkl+1)) then
         if(lshelkl.eq.1.or.lshelkl.eq.3) then
           do 103 nij=ijbeg,ijend
           call amtfer(bkl1(1,nij,nkl),vij(1,nij,1),mnbls)
  103      continue
         endif
         if(lshelkl.eq.2.or.lshelkl.eq.3) then
           do 104 nij=ijbeg,ijend
           call amtfer(bkl2(1,nij,nkl),uij(1,nij,1),mnbls)
  104      continue
         endif
         if(lshelkl.eq.3.and.nkl.eq.1) then
           do 105 nij=ijbeg,ijend
           call amtfer(bkl3(1,nij),sij(1,nij,1),mnbls)
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
       call amtfer(wij(1,ni,nj),xij(1,ni,nj,nkl),mnbls)
  107  continue
       if( nkl.le.nfu(nqkl+1)) then
           if(lshelkl.eq.1.or.lshelkl.eq.3) then
                do 1071 nj=njbeg,njend
                do 1071 ni=nibeg,niend
                call amtfer(vij(1,ni,nj),yij(1,ni,nj,nkl),mnbls)
 1071           continue
           endif
           if(lshelkl.eq.2.or.lshelkl.eq.3) then
                do 1072 nj=njbeg,njend
                do 1072 ni=nibeg,niend
                call amtfer(uij(1,ni,nj),zij(1,ni,nj,nkl),mnbls)
 1072           continue
           endif
       endif
c
       if(lshelij.eq.1.or.lshelij.eq.3) then
          if(jtyp.eq.1) then
             call amtfer(bij1(1,1,nkl),xij(1,1,1,nkl),mnbls)
          else
      call daxpy3(mnbls,xab,
     *    xij(1,1,2,nkl),xij(1,1,3,nkl),xij(1,1,4,nkl),
     *    bij1(1,2,nkl),bij1(1,3,nkl),bij1(1,4,nkl),bij1(1,1,nkl))
          endif
       endif
       if(lshelij.eq.2.or.lshelij.eq.3) then
           do 108 nij=nfu(nqi )+1,nfu(nqi +1)
           call amtfer(bij2(1,nij,nkl),xij(1,nij,1,nkl),mnbls)
  108      continue
       endif
       if(lshelij.eq.3) then
           call amtfer(bij3(1,nkl),xij(1,1,1,nkl),mnbls)
       endif
c*****
      if( nkl.le.nfu(nqkl+1)) then
c****  2 or 3 l-shells ****
c
           if(lshelij.eq.1) then
              if(jtyp.eq.1) then
                   call amtfer(b2l1(1,1,nkl),x2l1(1,1,nkl),mnbls)
                   call amtfer(b2l2(1,1,nkl),x2l2(1,1,nkl),mnbls)
c-- test ?      if(nkl.eq.1) then
                if(nkl.eq.1 .and. lcas3(3).eq.1) then
                   call amtfer(b3l34(1,1),x3l34(1,1),mnbls)
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
                   call amtfer(b2l3(1,nij,nkl),x2l3(1,nij,nkl),mnbls)
                   call amtfer(b2l4(1,nij,nkl),x2l4(1,nij,nkl),mnbls)
  109           continue
c-- test ?      if(nkl.eq.1) then
                if(nkl.eq.1 .and. lcas3(4).eq.1) then
                   do 110 nij=nfu(nqi )+1,nfu(nqi +1)
                   call amtfer(b3l34(1,nij),x3l34(1,nij),mnbls)
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
                  call amtfer(b2l3(1,2,nkl),x2l3(1,2,nkl),mnbls)
                  call amtfer(b2l3(1,3,nkl),x2l3(1,3,nkl),mnbls)
                  call amtfer(b2l3(1,4,nkl),x2l3(1,4,nkl),mnbls)
                endif
                if(lcas2(2).eq.1) then
                  call daxpy3(mnbls,xab,
     *            x2l2(1,2,nkl),x2l2(1,3,nkl),x2l2(1,4,nkl),
     *            b2l2(1,2,nkl),b2l2(1,3,nkl),b2l2(1,4,nkl),
     *            b2l2(1,1,nkl))
cc
                  call amtfer(b2l4(1,2,nkl),x2l4(1,2,nkl),mnbls)
                  call amtfer(b2l4(1,3,nkl),x2l4(1,3,nkl),mnbls)
                  call amtfer(b2l4(1,4,nkl),x2l4(1,4,nkl),mnbls)
                endif
                if(lcas3(1).eq.1) then
                  call amtfer(b3l12(1,nkl),x3l12(1,nkl),mnbls)
                endif
                if(lcas3(2).eq.1) then
                  call amtfer(b3l12(1,nkl),x3l12(1,nkl),mnbls)
                endif
           endif
      endif
c
  100 continue
c
c***********************************************************
c*    this part    shifts the angular momentum
c*         from position 3 to position 4
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
         call amtfer(xij(1,ni,nj,nkl),wij(1,nkl,1),mnbls)
  301    continue
c
ccccccc
        call horiz12(wij,lt6,lt7,xcd,mnbls,nqk,nql,nskl1)
ccccccc
cc
       if(lshelkl.eq.1.or.lshelkl.eq.3) then
         if(ltyp.eq.1) then
           call amtfer(yij(1,ni,nj,1),wij(1,1,1),mnbls)
         else
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     *   yij(1,ni,nj,2),yij(1,ni,nj,3),yij(1,ni,nj,4),yij(1,ni,nj,1))
         endif
       endif
       if(lshelkl.eq.2.or.lshelkl.eq.3) then
         do 312 nkl=nfu(nqkx)+1,nfu(nqkl+1)
         call amtfer(zij(1,ni,nj,nkl),wij(1,nkl,1),mnbls)
  312    continue
       endif
c
       if(lshelkl.eq.3) then
         call amtfer(sij(1,ni,nj),wij(1,1,1),mnbls)
       endif
c
       if( ni.le.nfu(nqi +1).and.nj.le.nfu(nqj +1) ) then
c
c****  2,3  l-shells ****
         if(lshelkl.eq.1) then
            if(ltyp.eq.1) then
               if(ni.eq.1.and.nj.ge.nfu(nqj)+1) then
                  call amtfer(x2l1(1,nj,1),wij(1,1,1),mnbls)
               endif
               if(nj.eq.1.and.ni.ge.nfu(nqi)+1) then
                  call amtfer(x2l3(1,ni,1),wij(1,1,1),mnbls)
               endif
c--test ?      if(ni.eq.1.and.nj.eq.1) then
               if(ni.eq.1.and.nj.eq.1 .and. lcas3(1).eq.1) then
                  call amtfer(x3l12(1,1),wij(1,1,1),mnbls)
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
                  call amtfer(x2l2(1,nj,nkl),wij(1,nkl,1),mnbls)
              endif
              if(nj.eq.1)  then
                  call amtfer(x2l4(1,ni,nkl),wij(1,nkl,1),mnbls)
              endif
c--test ?     if(ni.eq.1.and.nj.eq.1) then
              if(ni.eq.1.and.nj.eq.1 .and. lcas3(2).eq.1) then
                  call amtfer(x3l12(1,nkl),wij(1,nkl,1),mnbls)
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
              call amtfer(x2l2(1,nj,2),wij(1,2,1),mnbls)
              call amtfer(x2l2(1,nj,3),wij(1,3,1),mnbls)
              call amtfer(x2l2(1,nj,4),wij(1,4,1),mnbls)
              endif
              if(lcas2(3).eq.1.and.nj.eq.1) then
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     *  x2l3(1,ni,2),x2l3(1,ni,3),x2l3(1,ni,4),x2l3(1,ni,1))
cc
              call amtfer(x2l4(1,ni,2),wij(1,2,1),mnbls)
              call amtfer(x2l4(1,ni,3),wij(1,3,1),mnbls)
              call amtfer(x2l4(1,ni,4),wij(1,4,1),mnbls)
              endif
c
              if(lcas3(3).eq.1.and.ni.eq.1) then
                  call amtfer(x3l34(1,nj),wij(1,1,1),mnbls)         
              endif
              if(lcas3(4).eq.1.and.nj.eq.1) then
                  call amtfer(x3l34(1,ni),wij(1,1,1),mnbls)
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
      call amtfer(wij(1,nk,nl),buf(1,indx),mnbls)
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
#include "texas_lpar.fh"
ctest
      character*11 scftype
      character*8 where
      common /runtype/ scftype,where
ctest
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbex,klbex
      common/shell/lshellt,lshelij,lshelkl,lhelp,lcas2(4),lcas3(4)
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
c   ? dimension ssss(nbls)
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
ccccc  nqjx=1
       nibeg=1
       njbeg=1
c
       nqkx=1
ccccc  nqlx=1
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
           call amtfer(buf2(1,nij,nkl),wij(1,nij,1),mnbls)
  102      continue
c
       if( nkl.le.nfu(nqkl+1)) then
           do 103 nij=ijbeg,ijend
           call amtfer(bkl1(1,nij,nkl),vij(1,nij,1),mnbls)
  103      continue
           do 104 nij=ijbeg,ijend
           call amtfer(bkl2(1,nij,nkl),uij(1,nij,1),mnbls)
  104      continue
         if(                 nkl.eq.1) then
           do 105 nij=ijbeg,ijend
           call amtfer(bkl3(1,nij),sij(1,nij,1),mnbls)
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
        call amtfer(wij(1,ni,nj),xij(1,ni,nj,nkl),mnbls)
  107   continue
c
      if( nkl.le.nfu(nqkl+1)) then
        do 1071 nj=njbeg,njend
        do 1071 ni=nibeg,niend
        call amtfer(vij(1,ni,nj),yij(1,ni,nj,nkl),mnbls)
        call amtfer(uij(1,ni,nj),zij(1,ni,nj,nkl),mnbls)
 1071   continue
      endif
c
      call daxpy3(mnbls,xab,
     *   xij(1,1,2,nkl),xij(1,1,3,nkl),xij(1,1,4,nkl),
     *   bij1(1,2,nkl),bij1(1,3,nkl),bij1(1,4,nkl),bij1(1,1,nkl))
c
           do 108 nij=nfu(nqi )+1,nfu(nqi +1)
           call amtfer(bij2(1,nij,nkl),xij(1,nij,1,nkl),mnbls)
  108      continue
c
           call amtfer(bij3(1,nkl),xij(1,1,1,nkl),mnbls)
c
c*****
      if( nkl.le.nfu(nqkl+1)) then
      call daxpy3(mnbls,xab,x2l1(1,2,nkl),x2l1(1,3,nkl),x2l1(1,4,nkl),
     *    b2l1(1,2,nkl),b2l1(1,3,nkl),b2l1(1,4,nkl),b2l1(1,1,nkl))
                call amtfer(b2l3(1,2,nkl),x2l3(1,2,nkl),mnbls)
                call amtfer(b2l3(1,3,nkl),x2l3(1,3,nkl),mnbls)
                call amtfer(b2l3(1,4,nkl),x2l3(1,4,nkl),mnbls)
cc
      call daxpy3(mnbls,xab,x2l2(1,2,nkl),x2l2(1,3,nkl),x2l2(1,4,nkl),
     *  b2l2(1,2,nkl),b2l2(1,3,nkl),b2l2(1,4,nkl),b2l2(1,1,nkl))
cc
                call amtfer(b2l4(1,2,nkl),x2l4(1,2,nkl),mnbls)
                call amtfer(b2l4(1,3,nkl),x2l4(1,3,nkl),mnbls)
                call amtfer(b2l4(1,4,nkl),x2l4(1,4,nkl),mnbls)
c
                call amtfer(b3l1(1,nkl),x3l1(1,nkl),mnbls)
                call amtfer(b3l2(1,nkl),x3l2(1,nkl),mnbls)
c
      endif
      if(nkl.eq.1) then
c
      call daxpy3(mnbls,xab,x3l3(1,2),x3l3(1,3),x3l3(1,4),
     *  b3l3(1,2),b3l3(1,3),b3l3(1,4),b3l3(1,1))
c
                do 1101 nij=nfu(nqix)+1,nfu(nqij+1)
                call amtfer(b3l4(1,nij ),x3l4(1,nij),mnbls)
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
         call amtfer(xij(1,ni,nj,nkl),wij(1,nkl,1),mnbls)
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
         call amtfer(zij(1,ni,nj,nkl),wij(1,nkl,1),mnbls)
  312    continue
c
         call amtfer(sij(1,ni,nj),wij(1,1,1),mnbls)
c
c****
       if( ni.le.nfu(nqi +1).and.nj.le.nfu(nqj +1) ) then
              if(ni.eq.1 .and. nj.ge.2) then
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     *  x2l1(1,nj,2),x2l1(1,nj,3),x2l1(1,nj,4),x2l1(1,nj,1))
cc
              call amtfer(x2l2(1,nj,2),wij(1,2,1),mnbls)
              call amtfer(x2l2(1,nj,3),wij(1,3,1),mnbls)
              call amtfer(x2l2(1,nj,4),wij(1,4,1),mnbls)
cc
              endif
              if(nj.eq.1 .and. ni.ge.2) then
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     * x2l3(1,ni,2),x2l3(1,ni,3),x2l3(1,ni,4),x2l3(1,ni,1))
cc
              call amtfer(x2l4(1,ni,2),wij(1,2,1),mnbls)
              call amtfer(x2l4(1,ni,3),wij(1,3,1),mnbls)
              call amtfer(x2l4(1,ni,4),wij(1,4,1),mnbls)
              endif
c
              if(ni.eq.1) call amtfer(x3l3(1,nj),wij(1,1,1),mnbls)
              if(nj.eq.1) call amtfer(x3l4(1,ni),wij(1,1,1),mnbls)
       endif
       if( ni.eq.1.and.nj.eq.1 ) then
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     *  x3l1(1,2),x3l1(1,3),x3l1(1,4),x3l1(1,1))
cc
          call amtfer(x3l2(1,2),wij(1,2,1),mnbls)
          call amtfer(x3l2(1,3),wij(1,3,1),mnbls)
          call amtfer(x3l2(1,4),wij(1,4,1),mnbls)
cc
          call amtfer(ssss(1),wij(1,1,1),mnbls)
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
      call amtfer(wij(1,nk,nl),buf(1,indx),mnbls)
  305 continue
  300 continue
c
      return
      end
c=====================================================
      subroutine horiz12(wij,lw1,lw2,xab,mnbls,nqi,nqj,nsij1)
      implicit real*8 (a-h,o-z)
c
#include "texas_lpar.fh"
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
      subroutine shift0l_der1(buf,buf2,buf0,lt1,lt2,
     *                   wij,wij0,lt3,lsjl,
     *                   xij,xij0,lt4,lt5,lt6,mnbls,nbls,
     *                   xab,xcd,
     *                   indxx,lni1,lnj1,lnk1,lnl1)
c------------------------------------
c  when l-shells are not present
c------------------------------------
      implicit real*8 (a-h,o-z)
#include "texas_lpar.fh"
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
      common/shell/lshellt,lshelij,lshelkl,lhelp,lcas2(4),lcas3(4)
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbex,klbex
c
      dimension buf0(nbls,lt1,lt2),wij0(nbls,lt3,lsjl)
      dimension buf2(mnbls,lt1,lt2),wij(mnbls,lt3,lsjl)
      dimension xij(mnbls,lt4,lt5,lt6), xij0(nbls,lt4,lt5,lt6)
      dimension buf(mnbls,*)
      dimension xab(mnbls,3),xcd(mnbls,3)
      dimension indxx(lni1,lnj1,lnk1,lnl1)
c------------------------------------
c
      nibeg=nfu(nqi)+1
      niend=nfu(nqi+1)
      njbeg=nfu(nqj)+1
      njend=nfu(nqj+1)
cc
      nkbeg=nfu(nqk)+1
      nkend=nfu(nqk+1)
      nlbeg=nfu(nql)+1
      nlend=nfu(nql+1)
c
      nqix=nqi
ccccc nqjx=nqj
c
      nqkx=nqk
ccccc nqlx=nql
      nqklx=nqkl
c
      do 100 nkl=nfu(nqklx)+1,nfu(nskl+1)
c
           do 102 nij=nfu(nqix)+1,nfu(nsij+1)
           call amtfer(buf2(1,nij,nkl),wij(1,nij,1),mnbls)
           call amtfer(buf0(1,nij,nkl),wij0(1,nij,1),nbls)
  102      continue
c
           call horiz12_der1(wij0,wij,lt3,lsjl,xab, nbls,nqi,nqj,nsij1)
c
           do 107 nj=njbeg,njend
           do 107 ni=nibeg,niend
           call amtfer(wij(1,ni,nj),xij(1,ni,nj,nkl),mnbls)
           call amtfer(wij0(1,ni,nj),xij0(1,ni,nj,nkl), nbls)
  107      continue
  100 continue
c
c------------------------------------
c this part shifts angular momentum
c   from position 3 to position 4
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
         call amtfer( xij(1,ni,nj,nkl), wij(1,nkl,1),mnbls)
         call amtfer(xij0(1,ni,nj,nkl),wij0(1,nkl,1), nbls)
  301    continue
c
         call horiz34_der1(wij0,wij,lt3,lsjl,xcd, nbls,nqk,nql,nskl1)
c
      kxyz=0
      do 305 nk=nkbeg,nkend
      kxyz=kxyz+1
      lxyz=0
      do 305 nl=nlbeg,nlend
      lxyz=lxyz+1
      indx=indxx(ixyz,jxyz,kxyz,lxyz)
      call amtfer(wij(1,nk,nl),buf(1,indx),mnbls)
  305 continue
  300 continue
c
      end
c=============================================================
      subroutine horiz12_der1(wij0,wij,lw1,lw2,xab,nbls,nqi,nqj,nsij1)
      implicit real*8 (a-h,o-z)
c
#include "texas_lpar.fh"
c
      dimension wij0(nbls,lw1,lw2)
      dimension wij(9,nbls,lw1,lw2),xab(9,nbls,3)
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
                      do 130 n=1, nbls
c Ax
                      wij(1,n,ni,nj)=wij(1,n,nij,njm)
     *                              + xab(1,n,kcr)*wij(1,n,ni,njm)
c Bx
                      wij(2,n,ni,nj)=wij(2,n,nij,njm)
     *                              + xab(2,n,kcr)*wij(2,n,ni,njm)
c Cx 
                      wij(3,n,ni,nj)=wij(3,n,nij,njm)
     *                              + xab(3,n,kcr)*wij(3,n,ni,njm)
c Ay
                      wij(4,n,ni,nj)=wij(4,n,nij,njm)
     *                              + xab(4,n,kcr)*wij(4,n,ni,njm)
c By
                      wij(5,n,ni,nj)=wij(5,n,nij,njm)
     *                              + xab(5,n,kcr)*wij(5,n,ni,njm)
c Cy 
                      wij(6,n,ni,nj)=wij(6,n,nij,njm)
     *                              + xab(6,n,kcr)*wij(6,n,ni,njm)
c Az
                      wij(7,n,ni,nj)=wij(7,n,nij,njm)
     *                              + xab(7,n,kcr)*wij(7,n,ni,njm)
c Bz
                      wij(8,n,ni,nj)=wij(8,n,nij,njm)
     *                              + xab(8,n,kcr)*wij(8,n,ni,njm)
c Cz
                      wij(9,n,ni,nj)=wij(9,n,nij,njm)
     *                              + xab(9,n,kcr)*wij(9,n,ni,njm)
c
c add an extra term arising from differentiation of the shifting formula
c ONLY for derivatives over center A and B , not C :
c
                      wij0(n,ni,nj)=wij0(n,nij,njm)
     *                              + xab(1,n,kcr)*wij0(n,ni,njm)
                      if(kcr.eq.1) then
                         wij(1,n,ni,nj)=wij(1,n,ni,nj) + wij0(n,ni,njm)
                         wij(2,n,ni,nj)=wij(2,n,ni,nj) - wij0(n,ni,njm)
                      endif
                      if(kcr.eq.2) then
                         wij(4,n,ni,nj)=wij(4,n,ni,nj) + wij0(n,ni,njm)
                         wij(5,n,ni,nj)=wij(5,n,ni,nj) - wij0(n,ni,njm)
                      endif
                      if(kcr.eq.3) then
                         wij(7,n,ni,nj)=wij(7,n,ni,nj) + wij0(n,ni,njm)
                         wij(8,n,ni,nj)=wij(8,n,ni,nj) - wij0(n,ni,njm)
                      endif
  130                 continue
  125              continue
  120          continue
  115        continue
  110     continue
c
      end
c=====================================================
      subroutine horiz34_der1(wij0,wij,lw1,lw2,xab,nbls,nqi,nqj,nsij1)
      implicit real*8 (a-h,o-z)
c
#include "texas_lpar.fh"
c
      dimension wij0(nbls,lw1,lw2)
      dimension wij(9,nbls,lw1,lw2),xab(9,nbls,3)
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
                      do 130 n=1, nbls
c Ax
                      wij(1,n,ni,nj)=wij(1,n,nij,njm)
     *                              + xab(1,n,kcr)*wij(1,n,ni,njm)
c Bx
                      wij(2,n,ni,nj)=wij(2,n,nij,njm)
     *                              + xab(2,n,kcr)*wij(2,n,ni,njm)
c Cx 
                      wij(3,n,ni,nj)=wij(3,n,nij,njm)
     *                              + xab(3,n,kcr)*wij(3,n,ni,njm)
c Ay
                      wij(4,n,ni,nj)=wij(4,n,nij,njm)
     *                              + xab(4,n,kcr)*wij(4,n,ni,njm)
c By
                      wij(5,n,ni,nj)=wij(5,n,nij,njm)
     *                              + xab(5,n,kcr)*wij(5,n,ni,njm)
c Cy 
                      wij(6,n,ni,nj)=wij(6,n,nij,njm)
     *                              + xab(6,n,kcr)*wij(6,n,ni,njm)
c Az
                      wij(7,n,ni,nj)=wij(7,n,nij,njm)
     *                              + xab(7,n,kcr)*wij(7,n,ni,njm)
c Bz
                      wij(8,n,ni,nj)=wij(8,n,nij,njm)
     *                              + xab(8,n,kcr)*wij(8,n,ni,njm)
c Cz 
                      wij(9,n,ni,nj)=wij(9,n,nij,njm)
     *                              + xab(9,n,kcr)*wij(9,n,ni,njm)
c
c add an extra term arising from differentiation of the shifting formula
c ONLY for derivatives over center C, not over A and B:
c
                      wij0(n,ni,nj)=wij0(n,nij,njm)
     *                              + xab(1,n,kcr)*wij0(n,ni,njm)
                      if(kcr.eq.1) then
                         wij(3,n,ni,nj)=wij(3,n,ni,nj) + wij0(n,ni,njm)
                      endif
                      if(kcr.eq.2) then
                         wij(6,n,ni,nj)=wij(6,n,ni,nj) + wij0(n,ni,njm)
                      endif
                      if(kcr.eq.3) then
                         wij(9,n,ni,nj)=wij(9,n,ni,nj) + wij0(n,ni,njm)
                      endif
  130                 continue
  125              continue
  120          continue
  115        continue
  110     continue
c
      end
c=====================================================
      subroutine shift0l_der2(buf,buf2,buf1,buf0,lt1,lt2,
     *                   wij,wij1,wij0,lt3,lsjl,
     *                   xij,xij1,xij0,lt4,lt5,lt6,mnbls,nbls,
     *                   xab,xcd,
     *                   indxx,lni1,lnj1,lnk1,lnl1)
c------------------------------------
c  when l-shells are not present
c------------------------------------
      implicit real*8 (a-h,o-z)
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
      common/shell/lshellt,lshelij,lshelkl,lhelp,lcas2(4),lcas3(4)
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbex,klbex
#include "texas_lpar.fh"
c
      dimension buf0(nbls,lt1,lt2),wij0(nbls,lt3,lsjl)
      dimension buf1(9*nbls,lt1,lt2),wij1(9*nbls,lt3,lsjl)
      dimension buf2(mnbls,lt1,lt2),wij(mnbls,lt3,lsjl)
      dimension xij(mnbls,lt4,lt5,lt6),xij1(9*nbls,lt4,lt5,lt6),
     *                                 xij0(  nbls,lt4,lt5,lt6)
      dimension buf(mnbls,*)
      dimension xab(mnbls,3),xcd(mnbls,3)
      dimension indxx(lni1,lnj1,lnk1,lnl1)
c------------------------------------
      nbls9=9*nbls
c
      nibeg=nfu(nqi)+1
      niend=nfu(nqi+1)
      njbeg=nfu(nqj)+1
      njend=nfu(nqj+1)
cc
      nkbeg=nfu(nqk)+1
      nkend=nfu(nqk+1)
      nlbeg=nfu(nql)+1
      nlend=nfu(nql+1)
c
      nqix=nqi
ccccc nqjx=nqj
c
      nqkx=nqk
ccccc nqlx=nql
      nqklx=nqkl
c
      do 100 nkl=nfu(nqklx)+1,nfu(nskl+1)
c
           do 102 nij=nfu(nqix)+1,nfu(nsij+1)
           call amtfer(buf2(1,nij,nkl),wij(1,nij,1),mnbls)
           call amtfer(buf1(1,nij,nkl),wij1(1,nij,1),nbls9)
           call amtfer(buf0(1,nij,nkl),wij0(1,nij,1),nbls)
  102      continue
c
           call horiz12_der2(wij0,wij1,wij,
     *                       lt3,lsjl,xab, nbls,nqi,nqj,nsij1)
c
           do 107 nj=njbeg,njend
           do 107 ni=nibeg,niend
           call amtfer(wij(1,ni,nj),xij(1,ni,nj,nkl),mnbls)
           call amtfer(wij1(1,ni,nj),xij1(1,ni,nj,nkl),nbls9)
           call amtfer(wij0(1,ni,nj),xij0(1,ni,nj,nkl),nbls)
  107      continue
  100 continue
c
c------------------------------------
c this part shifts angular momentum
c   from position 3 to position 4
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
         call amtfer( xij(1,ni,nj,nkl), wij(1,nkl,1),mnbls)
         call amtfer(xij1(1,ni,nj,nkl),wij1(1,nkl,1),nbls9)
         call amtfer(xij0(1,ni,nj,nkl),wij0(1,nkl,1), nbls)
  301    continue
c
         call horiz34_der2(wij0,wij1,wij,
     *                     lt3,lsjl,xcd, nbls,nqk,nql,nskl1)
c
      kxyz=0
      do 305 nk=nkbeg,nkend
      kxyz=kxyz+1
      lxyz=0
      do 305 nl=nlbeg,nlend
      lxyz=lxyz+1
      indx=indxx(ixyz,jxyz,kxyz,lxyz)
      call amtfer(wij(1,nk,nl),buf(1,indx),mnbls)
  305 continue
  300 continue
c
      end
c=============================================================
      subroutine horiz12_der2(wij0,wij1,wij2,
     *                        lw1,lw2,xab,nbls,nqi,nqj,nsij1)
      implicit real*8 (a-h,o-z)
c
#include "texas_lpar.fh"
c
      dimension wij0(nbls,lw1,lw2)
      dimension wij1(9,nbls,lw1,lw2)
      dimension wij2(45,nbls,lw1,lw2),xab(45,nbls,3)
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
                      do 130 n=1, nbls
      wij0(n,ni,nj)=wij0(n,nij,njm) + xab(1,n,kcr)*wij0(n,ni,njm)
                      do m=1,9
      wij1(m,n,ni,nj)=wij1(m,n,nij,njm) + xab(m,n,kcr)*wij1(m,n,ni,njm)
                      enddo
                      do m=1,45
      wij2(m,n,ni,nj)=wij2(m,n,nij,njm) + xab(m,n,kcr)*wij2(m,n,ni,njm)
                      enddo
c
c add an extra term arising from differentiation of the shifting formula
c ONLY for derivatives over center A and B , not C :
c
c  first derivatives:
c
          if(kcr.eq.1) then
             wij1(1,n,ni,nj)=wij1(1,n,ni,nj) + wij0(n,ni,njm)
             wij1(2,n,ni,nj)=wij1(2,n,ni,nj) - wij0(n,ni,njm)
c
          endif
          if(kcr.eq.2) then
             wij1(4,n,ni,nj)=wij1(4,n,ni,nj) + wij0(n,ni,njm)
             wij1(5,n,ni,nj)=wij1(5,n,ni,nj) - wij0(n,ni,njm)
          endif
          if(kcr.eq.3) then
             wij1(7,n,ni,nj)=wij1(7,n,ni,nj) + wij0(n,ni,njm)
             wij1(8,n,ni,nj)=wij1(8,n,ni,nj) - wij0(n,ni,njm)
          endif
c
c second derivatives:
c
          if(kcr.eq.1) then
c block aa:
             wij2(1,n,ni,nj)=wij2(1,n,ni,nj) + wij1(1,n,ni,njm)*2.d0
             wij2(2,n,ni,nj)=wij2(2,n,ni,nj) + wij1(4,n,ni,njm)
             wij2(3,n,ni,nj)=wij2(3,n,ni,nj) + wij1(7,n,ni,njm)
c block ab:
             wij2(7,n,ni,nj) =wij2(7,n,ni,nj) - wij1(1,n,ni,njm)
     *                                        + wij1(2,n,ni,njm)
             wij2(8,n,ni,nj) =wij2(8,n,ni,nj) + wij1(5,n,ni,njm)
             wij2(9,n,ni,nj) =wij2(9,n,ni,nj) + wij1(8,n,ni,njm)
             wij2(10,n,ni,nj)=wij2(10,n,ni,nj)- wij1(4,n,ni,njm)
             wij2(13,n,ni,nj)=wij2(13,n,ni,nj)- wij1(7,n,ni,njm)
c block ac:
             wij2(16,n,ni,nj)=wij2(16,n,ni,nj)+ wij1(3,n,ni,njm)
             wij2(17,n,ni,nj)=wij2(17,n,ni,nj)+ wij1(6,n,ni,njm)
             wij2(18,n,ni,nj)=wij2(18,n,ni,nj)+ wij1(9,n,ni,njm)
c block bb:
             wij2(25,n,ni,nj)=wij2(25,n,ni,nj)- wij1(2,n,ni,njm)*2.d0
             wij2(26,n,ni,nj)=wij2(26,n,ni,nj)- wij1(5,n,ni,njm)
             wij2(27,n,ni,nj)=wij2(27,n,ni,nj)- wij1(8,n,ni,njm)
c block bc:
             wij2(31,n,ni,nj)=wij2(31,n,ni,nj)- wij1(3,n,ni,njm)
             wij2(32,n,ni,nj)=wij2(32,n,ni,nj)- wij1(6,n,ni,njm)
             wij2(33,n,ni,nj)=wij2(33,n,ni,nj)- wij1(9,n,ni,njm)
c block cc: no contributions of this kind
          endif
          if(kcr.eq.2) then
c block aa:
             wij2(2,n,ni,nj)=wij2(2,n,ni,nj) + wij1(1,n,ni,njm)
             wij2(4,n,ni,nj)=wij2(4,n,ni,nj) + wij1(4,n,ni,njm)*2.d0
             wij2(5,n,ni,nj)=wij2(5,n,ni,nj) + wij1(7,n,ni,njm)
c block ab:
             wij2(8,n,ni,nj) =wij2(8,n,ni,nj) - wij1(1,n,ni,njm)
             wij2(10,n,ni,nj)=wij2(10,n,ni,nj)+ wij1(2,n,ni,njm)
             wij2(11,n,ni,nj)=wij2(11,n,ni,nj)- wij1(4,n,ni,njm)
     *                                        + wij1(5,n,ni,njm)
             wij2(12,n,ni,nj)=wij2(12,n,ni,nj)+ wij1(8,n,ni,njm)
             wij2(14,n,ni,nj)=wij2(14,n,ni,nj)- wij1(7,n,ni,njm)
c block ac:
             wij2(19,n,ni,nj)=wij2(19,n,ni,nj)+ wij1(3,n,ni,njm)
             wij2(20,n,ni,nj)=wij2(20,n,ni,nj)+ wij1(6,n,ni,njm)
             wij2(21,n,ni,nj)=wij2(21,n,ni,nj)+ wij1(9,n,ni,njm)
c block bb:
             wij2(26,n,ni,nj)=wij2(26,n,ni,nj)- wij1(2,n,ni,njm)
             wij2(28,n,ni,nj)=wij2(28,n,ni,nj)- wij1(5,n,ni,njm)*2.d0
             wij2(29,n,ni,nj)=wij2(29,n,ni,nj)- wij1(8,n,ni,njm)
c block bc:
             wij2(34,n,ni,nj)=wij2(34,n,ni,nj)- wij1(3,n,ni,njm)
             wij2(35,n,ni,nj)=wij2(35,n,ni,nj)- wij1(6,n,ni,njm)
             wij2(36,n,ni,nj)=wij2(36,n,ni,nj)- wij1(9,n,ni,njm)
c block cc: no contributions of this kind
          endif
          if(kcr.eq.3) then
c block aa:
             wij2(3,n,ni,nj)=wij2(3,n,ni,nj) + wij1(1,n,ni,njm)
             wij2(5,n,ni,nj)=wij2(5,n,ni,nj) + wij1(4,n,ni,njm)
             wij2(6,n,ni,nj)=wij2(6,n,ni,nj) + wij1(7,n,ni,njm)*2.d0
c block ab:
             wij2(9,n,ni,nj) =wij2(9,n,ni,nj) - wij1(1,n,ni,njm)
             wij2(12,n,ni,nj)=wij2(12,n,ni,nj)- wij1(4,n,ni,njm)
             wij2(13,n,ni,nj)=wij2(13,n,ni,nj)+ wij1(2,n,ni,njm)
             wij2(14,n,ni,nj)=wij2(14,n,ni,nj)+ wij1(5,n,ni,njm)
             wij2(15,n,ni,nj)=wij2(15,n,ni,nj)- wij1(7,n,ni,njm)
     *                                        + wij1(8,n,ni,njm)
c block ac:
             wij2(22,n,ni,nj)=wij2(22,n,ni,nj)+ wij1(3,n,ni,njm)
             wij2(23,n,ni,nj)=wij2(23,n,ni,nj)+ wij1(6,n,ni,njm)
             wij2(24,n,ni,nj)=wij2(24,n,ni,nj)+ wij1(9,n,ni,njm)
c block bb:
             wij2(27,n,ni,nj)=wij2(27,n,ni,nj)- wij1(2,n,ni,njm)
             wij2(29,n,ni,nj)=wij2(29,n,ni,nj)- wij1(5,n,ni,njm)
             wij2(30,n,ni,nj)=wij2(30,n,ni,nj)- wij1(8,n,ni,njm)*2.d0
c block bc:
             wij2(37,n,ni,nj)=wij2(37,n,ni,nj)- wij1(3,n,ni,njm)
             wij2(38,n,ni,nj)=wij2(38,n,ni,nj)- wij1(6,n,ni,njm)
             wij2(39,n,ni,nj)=wij2(39,n,ni,nj)- wij1(9,n,ni,njm)
c block cc: no contributions of this kind
          endif
  130                 continue
  125              continue
  120          continue
  115        continue
  110     continue
c
      end
c=====================================================
      subroutine horiz34_der2(wij0,wij1,wij2,
     *                        lw1,lw2,xab,nbls,nqi,nqj,nsij1)
      implicit real*8 (a-h,o-z)
c
#include "texas_lpar.fh"
c
      dimension wij0(nbls,lw1,lw2)
      dimension wij1(9,nbls,lw1,lw2)
      dimension wij2(45,nbls,lw1,lw2),xab(45,nbls,3)
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
                      do 130 n=1, nbls
      wij0(n,ni,nj)=wij0(n,nij,njm) + xab(1,n,kcr)*wij0(n,ni,njm)
                         do m=1,9
      wij1(m,n,ni,nj)=wij1(m,n,nij,njm) + xab(m,n,kcr)*wij1(m,n,ni,njm)
                         enddo
                         do m=1,45
      wij2(m,n,ni,nj)=wij2(m,n,nij,njm) + xab(m,n,kcr)*wij2(m,n,ni,njm)
                         enddo
c
c add an extra term arising from differentiation of the shifting formula
c ONLY for derivatives over center C, not over A and B:
c
c first derivatives:
c
      if(kcr.eq.1) then
         wij1(3,n,ni,nj)=wij1(3,n,ni,nj) + wij0(n,ni,njm)
      endif
      if(kcr.eq.2) then
         wij1(6,n,ni,nj)=wij1(6,n,ni,nj) + wij0(n,ni,njm)
      endif
      if(kcr.eq.3) then
         wij1(9,n,ni,nj)=wij1(9,n,ni,nj) + wij0(n,ni,njm)
      endif
c
c second derivatives:
c
          if(kcr.eq.1) then
c block ac:
             wij2(16,n,ni,nj)=wij2(16,n,ni,nj)+ wij1(1,n,ni,njm)
             wij2(19,n,ni,nj)=wij2(19,n,ni,nj)+ wij1(4,n,ni,njm)
             wij2(22,n,ni,nj)=wij2(22,n,ni,nj)+ wij1(7,n,ni,njm)
c block bc:
             wij2(31,n,ni,nj)=wij2(31,n,ni,nj)+ wij1(2,n,ni,njm)
             wij2(34,n,ni,nj)=wij2(34,n,ni,nj)+ wij1(5,n,ni,njm)
             wij2(37,n,ni,nj)=wij2(37,n,ni,nj)+ wij1(8,n,ni,njm)
c block cc:
             wij2(40,n,ni,nj)=wij2(40,n,ni,nj)+ wij1(3,n,ni,njm)*2.d0
             wij2(41,n,ni,nj)=wij2(41,n,ni,nj)+ wij1(6,n,ni,njm)
             wij2(42,n,ni,nj)=wij2(42,n,ni,nj)+ wij1(9,n,ni,njm)
          endif
          if(kcr.eq.2) then
c block ac:
             wij2(17,n,ni,nj)=wij2(17,n,ni,nj)+ wij1(1,n,ni,njm)
             wij2(20,n,ni,nj)=wij2(20,n,ni,nj)+ wij1(4,n,ni,njm)
             wij2(23,n,ni,nj)=wij2(23,n,ni,nj)+ wij1(7,n,ni,njm)
c block bc:
             wij2(32,n,ni,nj)=wij2(32,n,ni,nj)+ wij1(2,n,ni,njm)
             wij2(35,n,ni,nj)=wij2(35,n,ni,nj)+ wij1(5,n,ni,njm)
             wij2(38,n,ni,nj)=wij2(38,n,ni,nj)+ wij1(8,n,ni,njm)
c block cc:
             wij2(41,n,ni,nj)=wij2(41,n,ni,nj)+ wij1(3,n,ni,njm)
             wij2(43,n,ni,nj)=wij2(43,n,ni,nj)+ wij1(6,n,ni,njm)*2.d0
             wij2(44,n,ni,nj)=wij2(44,n,ni,nj)+ wij1(9,n,ni,njm)
          endif
          if(kcr.eq.3) then
c block ac:
             wij2(18,n,ni,nj)=wij2(18,n,ni,nj)+ wij1(1,n,ni,njm)
             wij2(21,n,ni,nj)=wij2(21,n,ni,nj)+ wij1(4,n,ni,njm)
             wij2(24,n,ni,nj)=wij2(24,n,ni,nj)+ wij1(7,n,ni,njm)
c block bc:
             wij2(33,n,ni,nj)=wij2(33,n,ni,nj)+ wij1(2,n,ni,njm)
             wij2(36,n,ni,nj)=wij2(36,n,ni,nj)+ wij1(5,n,ni,njm)
             wij2(39,n,ni,nj)=wij2(39,n,ni,nj)+ wij1(8,n,ni,njm)
c block cc:
             wij2(42,n,ni,nj)=wij2(42,n,ni,nj)+ wij1(3,n,ni,njm)
             wij2(44,n,ni,nj)=wij2(44,n,ni,nj)+ wij1(6,n,ni,njm)
             wij2(45,n,ni,nj)=wij2(45,n,ni,nj)+ wij1(9,n,ni,njm)*2.d0
          endif
  130                 continue
  125              continue
  120          continue
  115        continue
  110     continue
c
      end
c=======================================================================
      subroutine make_indxx(indxx,lni,lnj,lnk,lnl)
      dimension indxx(lni,lnj,lnk,lnl)
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
      end
c=======================================================================
