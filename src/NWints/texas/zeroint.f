c $Id: zeroint.f,v 1.2 1996-05-14 17:53:52 d3g681 Exp $
c* This is used to substitute ZERO in Buffers when some
c* of quartets do not appear First time (they are neglected)
c*
      subroutine zeroint(bl,nbls,nbls1,nbls2,l01,l02,indx,ngcd,ibux )
      implicit real*8 (a-h,o-z)
c
      common /logic4/ nfu(1)
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
        ibut=ibuf
        call zerosp(lnijkl,nbls,bl(ibut),bl(idxnot),nbls2,ngcd)
        call retmem(1)
        return
      endif
c
C******************************************************
c-------------------------------
c-   --- for buf2  ---
c
        ibut=ibuf2
         call zerout(nbls,nbls2,ngcd, bl(idxnot),
     *               bl(ibut ),l01,l02,nfu(nqij)+1,nfu(nqkl)+1 )
      IF(lshellt.eq.0) go to 100
c-------------------------------
c     if(lcase.eq. 2.or.lcase.eq. 6.or.lcase.eq. 8.or.lcase.eq. 9.or.
c    *   lcase.eq.12.or.lcase.eq.13.or.lcase.eq.14.or.lcase.eq.16) then
      if(lshelij.eq.1 .or. lshelij.eq.3) then
c-   --- for bfij1 ---
         call zerout(nbls,nbls2,ngcd, bl(idxnot),
     *               bl(ibfij1),lqij,l02,ijbeg,klbeg)
c
      endif
c----------
c     if(lcase.eq. 3.or.lcase.eq. 6.or.lcase.eq.10.or.lcase.eq.11.or.
c    *   lcase.eq.12.or.lcase.eq.13.or.lcase.eq.15.or.lcase.eq.16) then
      if(lshelij.eq.2 .or. lshelij.eq.3) then
c-   --- for bfij2 ---
         call zerout(nbls,nbls2,ngcd, bl(idxnot),
     *               bl(ibfij2),lqij,l02,ijbeg,klbeg)
      endif
c----------
c     if(lcase.eq. 4.or.lcase.eq. 7.or.lcase.eq. 8.or.lcase.eq.10.or.
c    *   lcase.eq.12.or.lcase.eq.14.or.lcase.eq.15.or.lcase.eq.16) then
      if(lshelkl.eq.1 .or. lshelkl.eq.3) then
c-   --- for bfkl1 ---
         call zerout(nbls,nbls2,ngcd, bl(idxnot),
     *               bl(ibfkl1),l01,lqkl,ijbeg,klbeg)
      endif
c----------
c     if(lcase.eq. 5.or.lcase.eq. 7.or.lcase.eq. 9.or.lcase.eq.11.or.
c    *   lcase.eq.13.or.lcase.eq.14.or.lcase.eq.15.or.lcase.eq.16) then
      if(lshelkl.eq.2 .or. lshelkl.eq.3) then
c-   --- for bfkl2 ---
         call zerout(nbls,nbls2,ngcd, bl(idxnot),
     *               bl(ibfkl2),l01,lqkl,ijbeg,klbeg)
      endif
c
      IF(lshellt.eq.1) go to 100
c----------
c     if(lcase.eq. 6.or.lcase.eq.12.or.lcase.eq.13.or.lcase.eq.16) then
      if(lshelij.eq.3) then
c-   --- for bfij3 ---
         call zerout(nbls,nbls2,ngcd, bl(idxnot),
     *               bl(ibfij3),lij3,l02,ijbeg,klbeg)
      endif
c----------
c     if(lcase.eq. 7.or.lcase.eq.14.or.lcase.eq.15.or.lcase.eq.16) then
      if(lshelkl.eq.3) then
c-   --- for bfkl3 ---
         call zerout(nbls,nbls2,ngcd, bl(idxnot),
     *               bl(ibfkl3),l01,lkl3,ijbeg,klbeg)
      endif
c----------
c     if(lcase.eq. 8.or.lcase.eq.12.or.lcase.eq.14.or.lcase.eq.16) then
c-   --- for bf2l1 ---
      if(lcas2(1).eq.1) then
         call zerout(nbls,nbls2,ngcd, bl(idxnot),
     *               bl(ibf2l1),lqij,lqkl,ijbeg,klbeg)
      endif
c----------
c     if(lcase.eq. 9.or.lcase.eq.13.or.lcase.eq.14.or.lcase.eq.16) then
      if(lcas2(2).eq.1) then
c-   --- for bf2l2 ---
         call zerout(nbls,nbls2,ngcd, bl(idxnot),
     *               bl(ibf2l2),lqij,lqkl,ijbeg,klbeg)
      endif
c----------
c     if(lcase.eq.10.or.lcase.eq.12.or.lcase.eq.15.or.lcase.eq.16) then
      if(lcas2(3).eq.1) then
c-   --- for bf2l3 ---
         call zerout(nbls,nbls2,ngcd, bl(idxnot),
     *               bl(ibf2l3),lqij,lqkl,ijbeg,klbeg)
      endif
c----------
c     if(lcase.eq.11.or.lcase.eq.13.or.lcase.eq.15.or.lcase.eq.16) then
      if(lcas2(4).eq.1) then
c-   --- for bf2l4 ---
         call zerout(nbls,nbls2,ngcd, bl(idxnot),
     *               bl(ibf2l4),lqij,lqkl,ijbeg,klbeg)
      endif
c
      IF(lshellt.eq.2) go to 100
c----------
c     if(lcase.eq.12.or.lcase.eq.16) then
      if(lcas3(1).eq.1) then
c-   --- for bf3l  ---
         call zerout(nbls,nbls2,ngcd, bl(idxnot),
     *               bl(ibf3l1),lqmx,l3l,1   ,1)
      endif
c     if(lcase.eq.13.or.lcase.eq.16) then
      if(lcas3(2).eq.1) then
         call zerout(nbls,nbls2,ngcd, bl(idxnot),
     *               bl(ibf3l2),lqmx,l3l,1   ,1)
      endif
c     if(lcase.eq.14.or.lcase.eq.16) then
      if(lcas3(3).eq.1) then
         call zerout(nbls,nbls2,ngcd, bl(idxnot),
     *               bl(ibf3l3),l3l,lqmx,1   ,1)
      endif
c     if(lcase.eq.15.or.lcase.eq.16) then
      if(lcas3(4).eq.1) then
         call zerout(nbls,nbls2,ngcd, bl(idxnot),
     *               bl(ibf3l4),l3l,lqmx,1   ,1)
      endif
c
      IF(lshellt.eq.3) go to 100
c----------
      if(lcase.eq.16) then
c-   --- for ssss(nbls)  ---
         call zerout(nbls,nbls2,ngcd, bl(idxnot),
     *               bl(issss ),lsss,lsss,1   ,1)
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
      data zero /0.d0/
c
            do 10 iqu=1,ngcd
            do 10 icx=1,lnijkl
            do 10 i=1,nbls2
            ijkl=idxnot(i)
            buf(ijkl,icx,iqu)=zero
  10        continue
c
      end
c=================================================
      subroutine zerout(nbls,nbls2,ngcd, idxnot, azero,l1,l2,i1,i2)
      implicit real*8 (a-h,o-z)
      dimension idxnot(*)
      dimension azero(nbls,l1,l2,ngcd)
c     
      if (ngcd .gt. 1) then
         do iqu=1,ngcd
            do kl=i2,l2
               do ij=i1,l1
                  do i=1,nbls2
                     ijkl=idxnot(i)
                     azero(ijkl,ij,kl,iqu)=0.0d0
                  enddo
               enddo
            enddo
         enddo
      else
         do kl=i2,l2
            do ij=i1,l1
               do i=1,nbls2
                  ijkl=idxnot(i)
                  azero(ijkl,ij,kl,1)=0.0d0
               enddo
            enddo
         enddo
      endif
c     
      end
c=================================================
