c $Id: zeroint.f 19696 2010-10-29 16:53:42Z d3y133 $
c* This is used to substitute ZERO in Buffers when some
c* of quartets do not appear First time (they are neglected)
c*
c* Changes for gradient derivatives :
c* one more parameter added to show how many arrays have to be zerout:
c* it can be 4 for gradient der. otherwise it is 1
c*
      subroutine zeroint(bl,nbls,nbls1,nbls2,l01,l02,indx,ngcd ,narr)
      implicit real*8 (a-h,o-z)
c
      character*11 scftype
      character*8 where
      common /runtype/ scftype,where
#include "texas_lpar.fh"      
      COMMON/SHELL/LSHELLT,LSHELIJ,LSHELKL,LHELP,LCAS2(4),LCAS3(4)
      common /lcases/ lcase
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
c
      common /memor4/ iwt0,iwt1,iwt2,ibuf,ibuf2,
     * ibfij1,ibfij2,ibfkl1,ibfkl2,
     * ibf2l1,ibf2l2,ibf2l3,ibf2l4,ibfij3,ibfkl3,
     * ibf3l,issss,
     * ix2l1,ix2l2,ix2l3,ix2l4,ix3l1,ix3l2,ix3l3,ix3l4,
     * ixij,iyij,izij, iwij,ivij,iuij,isij
c
      common /memor4a/ ibf3l1,ibf3l2,ibf3l3,ibf3l4
c
c dimensions for assembling :
      common /dimasse/ lqij,lqkl,lqmx,lij3,lkl3,l3l,lsss
c
      dimension bl(*)
C******************************************************
c
      call getmem(nbls,idxnot)
      call getmem(nbls,idpres)
c
      call setnot(nbls,nbls1,bl(indx),bl(idxnot),bl(idpres))
      call retmem(1)
c
      if(mmax.le.2) then
        call zerosp(lnijkl,nbls,bl(ibuf),bl(idxnot),nbls2,ngcd)
        call retmem(1)
        return
      endif
c
C******************************************************
c-------------------------------
c-   --- for buf2  ---
c
         ibut2=ibuf2
         if(where.eq.'forc' .or. where.eq.'hess') ibut2=ibuf
         call zerout(narr,nbls,nbls2,ngcd, bl(idxnot),
     *               bl(ibut2),l01,l02,nfu(nqij)+1,nfu(nqkl)+1 )
      IF(lshellt.eq.0) go to 100
c-------------------------------
c     if(lcase.eq. 2.or.lcase.eq. 6.or.lcase.eq. 8.or.lcase.eq. 9.or.
c    *   lcase.eq.12.or.lcase.eq.13.or.lcase.eq.14.or.lcase.eq.16) then
      if(lshelij.eq.1 .or. lshelij.eq.3) then
c-   --- for bfij1 ---
c
         ibut2=ibfij1
ccc      incre=nbls*lqij*l02
ccc      do 15 iarr=1,narr
            call zerout(narr,nbls,nbls2,ngcd, bl(idxnot),
ccc  *                  bl(ibfij1),lqij,l02,ijbeg,klbeg)
     *                  bl(ibut2 ),lqij,l02,ijbeg,klbeg)
ccc         ibut2=ibut2+incre
ccc15    continue
c
      endif
c----------
c     if(lcase.eq. 3.or.lcase.eq. 6.or.lcase.eq.10.or.lcase.eq.11.or.
c    *   lcase.eq.12.or.lcase.eq.13.or.lcase.eq.15.or.lcase.eq.16) then
      if(lshelij.eq.2 .or. lshelij.eq.3) then
c-   --- for bfij2 ---
c
         ibut2=ibfij2
ccc      incre=nbls*lqij*l02
ccc      do 20 iarr=1,narr
            call zerout(narr,nbls,nbls2,ngcd, bl(idxnot),
ccc  *                  bl(ibfij2),lqij,l02,ijbeg,klbeg)
     *                  bl(ibut2 ),lqij,l02,ijbeg,klbeg)
ccc         ibut2=ibut2+incre
ccc20    continue
      endif
c----------
c     if(lcase.eq. 4.or.lcase.eq. 7.or.lcase.eq. 8.or.lcase.eq.10.or.
c    *   lcase.eq.12.or.lcase.eq.14.or.lcase.eq.15.or.lcase.eq.16) then
      if(lshelkl.eq.1 .or. lshelkl.eq.3) then
c-   --- for bfkl1 ---
c
         ibut2=ibfkl1
ccc      incre=nbls*l01*lqkl
ccc      do 25 iarr=1,narr
            call zerout(narr,nbls,nbls2,ngcd, bl(idxnot),
cccc *                  bl(ibfkl1),l01,lqkl,ijbeg,klbeg)
     *                  bl(ibut2 ),l01,lqkl,ijbeg,klbeg)
ccc         ibut2=ibut2+incre
ccc25    continue
      endif
c----------
c     if(lcase.eq. 5.or.lcase.eq. 7.or.lcase.eq. 9.or.lcase.eq.11.or.
c    *   lcase.eq.13.or.lcase.eq.14.or.lcase.eq.15.or.lcase.eq.16) then
      if(lshelkl.eq.2 .or. lshelkl.eq.3) then
c-   --- for bfkl2 ---
c
         ibut2=ibfkl2
ccc      incre=nbls*l01*lqkl
ccc      do 30 iarr=1,narr
            call zerout(narr,nbls,nbls2,ngcd, bl(idxnot),
ccc  *                  bl(ibfkl2),l01,lqkl,ijbeg,klbeg)
     *                  bl(ibut2 ),l01,lqkl,ijbeg,klbeg)
ccc         ibut2=ibut2+incre
ccc30    continue
      endif
c
      IF(lshellt.eq.1) go to 100
c----------
c     if(lcase.eq. 6.or.lcase.eq.12.or.lcase.eq.13.or.lcase.eq.16) then
      if(lshelij.eq.3) then
c-   --- for bfij3 ---
c
         ibut2=ibfij3
ccc      incre=nbls*lij3*l02
ccc      do 35 iarr=1,narr
            call zerout(narr,nbls,nbls2,ngcd, bl(idxnot),
ccc  *                  bl(ibfij3),lij3,l02,ijbeg,klbeg)
     *                  bl(ibut2 ),lij3,l02,ijbeg,klbeg)
ccc         ibut2=ibut2+incre
ccc35    continue
      endif
c----------
c     if(lcase.eq. 7.or.lcase.eq.14.or.lcase.eq.15.or.lcase.eq.16) then
      if(lshelkl.eq.3) then
c-   --- for bfkl3 ---
c
         ibut2=ibfkl3
ccc      incre=nbls*l01*lkl3
ccc      do 40 iarr=1,narr
            call zerout(narr,nbls,nbls2,ngcd, bl(idxnot),
ccc  *                  bl(ibfkl3),l01,lkl3,ijbeg,klbeg)
     *                  bl(ibut2 ),l01,lkl3,ijbeg,klbeg)
ccc         ibut2=ibut2+incre
ccc40    continue
      endif
c----------
c     if(lcase.eq. 8.or.lcase.eq.12.or.lcase.eq.14.or.lcase.eq.16) then
c-   --- for bf2l1 ---
      if(lcas2(1).eq.1) then
c
         ibut2=ibf2l1
ccc      incre=nbls*lqij*lqkl
ccc      do 45 iarr=1,narr
            call zerout(narr,nbls,nbls2,ngcd, bl(idxnot),
ccc  *                  bl(ibf2l1),lqij,lqkl,ijbeg,klbeg)
     *                  bl(ibut2 ),lqij,lqkl,ijbeg,klbeg)
ccc         ibut2=ibut2+incre
ccc45    continue
      endif
c----------
c     if(lcase.eq. 9.or.lcase.eq.13.or.lcase.eq.14.or.lcase.eq.16) then
      if(lcas2(2).eq.1) then
c-   --- for bf2l2 ---
c
         ibut2=ibf2l2
ccc      incre=nbls*lqij*lqkl
ccc      do 50 iarr=1,narr
            call zerout(narr,nbls,nbls2,ngcd, bl(idxnot),
ccc  *                  bl(ibf2l2),lqij,lqkl,ijbeg,klbeg)
     *                  bl(ibut2 ),lqij,lqkl,ijbeg,klbeg)
ccc         ibut2=ibut2+incre
ccc50    continue
      endif
c----------
c     if(lcase.eq.10.or.lcase.eq.12.or.lcase.eq.15.or.lcase.eq.16) then
      if(lcas2(3).eq.1) then
c-   --- for bf2l3 ---
c
         ibut2=ibf2l3
ccc      incre=nbls*lqij*lqkl
ccc      do 55 iarr=1,narr
            call zerout(narr,nbls,nbls2,ngcd, bl(idxnot),
ccc  *                  bl(ibf2l3),lqij,lqkl,ijbeg,klbeg)
     *                  bl(ibut2 ),lqij,lqkl,ijbeg,klbeg)
ccc         ibut2=ibut2+incre
ccc55    continue
      endif
c----------
c     if(lcase.eq.11.or.lcase.eq.13.or.lcase.eq.15.or.lcase.eq.16) then
      if(lcas2(4).eq.1) then
c-   --- for bf2l4 ---
c
         ibut2=ibf2l4
ccc      incre=nbls*lqij*lqkl
ccc      do 60 iarr=1,narr
            call zerout(narr,nbls,nbls2,ngcd, bl(idxnot),
ccc  *                  bl(ibf2l4),lqij,lqkl,ijbeg,klbeg)
     *                  bl(ibut2 ),lqij,lqkl,ijbeg,klbeg)
ccc         ibut2=ibut2+incre
ccc60    continue
      endif
c
      IF(lshellt.eq.2) go to 100
c----------
c     if(lcase.eq.12.or.lcase.eq.16) then
      if(lcas3(1).eq.1) then
c-   --- for bf3l  ---
c
         ibut2=ibf3l1
ccc      incre=nbls*lqmx*l3l
ccc      do 65 iarr=1,narr
            call zerout(narr,nbls,nbls2,ngcd, bl(idxnot),
ccc  *                  bl(ibf3l1),lqmx,l3l,1   ,1)
     *                  bl(ibut2 ),lqmx,l3l,1   ,1)
ccc         ibut2=ibut2+incre
ccc65    continue
      endif
c     if(lcase.eq.13.or.lcase.eq.16) then
      if(lcas3(2).eq.1) then
c
         ibut2=ibf3l2
ccc      incre=nbls*lqmx*l3l
ccc      do 70 iarr=1,narr
            call zerout(narr,nbls,nbls2,ngcd, bl(idxnot),
ccc  *                  bl(ibf3l2),lqmx,l3l,1   ,1)
     *                  bl(ibut2 ),lqmx,l3l,1   ,1)
ccc         ibut2=ibut2+incre
ccc70    continue
      endif
c     if(lcase.eq.14.or.lcase.eq.16) then
      if(lcas3(3).eq.1) then
c
         ibut2=ibf3l3
ccc      incre=nbls*lqmx*l3l
ccc      do 75 iarr=1,narr
            call zerout(narr,nbls,nbls2,ngcd, bl(idxnot),
ccc  *                  bl(ibf3l3),l3l,lqmx,1   ,1)
     *                  bl(ibut2 ),l3l,lqmx,1   ,1)
ccc         ibut2=ibut2+incre
ccc75    continue
      endif
c     if(lcase.eq.15.or.lcase.eq.16) then
      if(lcas3(4).eq.1) then
c
         ibut2=ibf3l4
ccc      incre=nbls*lqmx*l3l
ccc      do 80 iarr=1,narr
            call zerout(narr,nbls,nbls2,ngcd, bl(idxnot),
ccc  *                  bl(ibf3l4),l3l,lqmx,1   ,1)
     *                  bl(ibut2 ),l3l,lqmx,1   ,1)
ccc         ibut2=ibut2+incre
ccc80    continue
      endif
c
      IF(lshellt.eq.3) go to 100
c----------
      if(lcase.eq.16) then
c-   --- for ssss(nbls)  ---
c
         ibut2=issss
ccc      incre=nbls*lsss*lsss
ccc      do 85 iarr=1,narr
         call zerout(narr,nbls,nbls2,ngcd, bl(idxnot),
ccc  *               bl(issss ),lsss,lsss,1   ,1)
     *               bl(ibut2 ),lsss,lsss,1   ,1)
ccc         ibut2=ibut2+incre
ccc85    continue
      endif
c
c-------------------------------
  100 CONTINUE
      call retmem(1)
c
      end
c**************************************************************
      subroutine setnot(nbls,nbls1,indx,idxnot,idpres)
      dimension indx(*),idxnot(*),idpres(*)
c
      do 10 i=1,nbls
      idpres(i)=0
   10 continue
c
      do 20 i=1,nbls1
      ijkl=indx(i)
      idpres(ijkl)=1
   20 continue
c
      ijkl2=0
      do 30 i=1,nbls
      if(idpres(i).eq.1 ) go to 30
      ijkl2=ijkl2+1
      idxnot(ijkl2)=i
   30 continue
c
      return
      end
c=================================================
      subroutine zerosp(lnijkl,nbls,buf,idxnot,nbls2,ngcd)
      implicit real*8 (a-h,o-z)
      dimension idxnot(*)
      dimension buf(nbls,lnijkl,ngcd)
c
            do 10 iqu=1,ngcd
            do 10 icx=1,lnijkl
            do 10 i=1,nbls2
            ijkl=idxnot(i)
            buf(ijkl,icx,iqu)=0.0d0
  10        continue
c
      end
c=================================================
cccc  subroutine zerout(nbls,nbls2,ngcd, idxnot, azero,l1,l2,i1,i2)
      subroutine zerout(narr,nbls,nbls2,ngcd, idxnot, azero,l1,l2,i1,i2)
      implicit real*8 (a-h,o-z)
      dimension idxnot(*)
c2002 dimension azero(narr,nbls,l1,l2,ngcd)
      dimension azero(nbls,l1,l2,ngcd,narr)
c     
      if (ngcd .gt. 1) then
       do iarr=1,narr
         do iqu=1,ngcd
            do kl=i2,l2
               do ij=i1,l1
                  do i=1,nbls2
                     ijkl=idxnot(i)
c................... do iarr=1,narr
c...................    azero(iarr,ijkl,ij,kl,iqu)=0.0d0
                        azero(ijkl,ij,kl,iqu,iarr)=0.0d0
c................... enddo
                  enddo
               enddo
            enddo
         enddo
       enddo
      else
       do iarr=1,narr
         do kl=i2,l2
            do ij=i1,l1
               do i=1,nbls2
                  ijkl=idxnot(i)
c................ do iarr=1,narr
c................    azero(iarr,ijkl,ij,kl,1)=0.0d0
                     azero(ijkl,ij,kl,1,iarr)=0.0d0
c................ enddo
               enddo
            enddo
         enddo
       enddo
      endif
c     
      end
c=================================================
