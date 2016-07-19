* $Id: derivat.f 19696 2010-10-29 16:53:42Z d3y133 $
c===================================================================
c This file contains all dintegral derivatives routines :
c for nmr/giao integral derivatives
c for first derivatives (gradient) integrals
c for second derivatives (hessian) integrals
c===================================================================
c
c
c===========NMR/GIAO INTEGRAL DERIVATIVES ROUTINES==================
c
c     It is called from Calcint2 when WHERE='shif'
c
      subroutine shift_der(nbls,lnijr,lnklr,npij,npkl,ngcd,
     *                   idx1,idx2, ixab,ixcd,ixyab,ixycd)
c
      implicit real*8 (a-h,o-z)
c
#include "texas_lpar.fh"      
      common /big/ bl(1)
      COMMON/SHELL/LSHELLT,LSHELIJ,LSHELKL,LHELP,LCAS2(4),LCAS3(4)
      common /lcases/ lcase
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
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
      common /dimasse/ lqijr,lqklr,lqmxr,lij3,lkl3,l3l,lsss
c----------------------------------------------------------
c     lqijr=nfu(nqij1+1)
c     lqklr=nfu(nqkl1+1)
c     lqmxr=lqijr
c       if(lqklr.gt.lqijr) lqmxr=lqklr
c
      lqij=nfu(nqij+1)
      lqkl=nfu(nqkl+1)
      lqmx=lqij
        if(lqkl.gt.lqij) lqmx=lqkl
c----------------------------------------------------------
c
      call getmem(3*nbls,ixabq)
      call getmem(3*nbls,ixcdq)
      call getmem(3*nbls,ixyabq)
      call getmem(3*nbls,ixycdq)
c------
c
      call conv24x(nbls,npij,npkl,bl(idx1),bl(idx2),
     *            bl(ixab),bl(ixcd),bl(ixyab),bl(ixycd),
     *            bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq) )
c------
c
      ibeg =ibuf2
      incr1=ngcd*nbls*lnijr*lnklr
      ider=ibeg+incr1
c
      call giao_der(ngcd,nbls,bl(ibeg),lnijr,lnklr,lnij,lnkl,nqij,nqkl,
     *            bl(ider),bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c
      ibuf2=ibuf2+incr1
c------
      if(lshellt.eq.0) go to 100
c------
c
c     if(lcase.eq. 2.or.lcase.eq. 6.or.lcase.eq. 8.or.lcase.eq. 9.or.
c    *   lcase.eq.12.or.lcase.eq.13.or.lcase.eq.14.or.lcase.eq.16) then
      if(lshelij.eq.1 .or. lshelij.eq.3) then
c-   --- for bfij1 ---
c
         ibeg =ibfij1
         incr1=nbls*lqijr*lnklr
         ider=ibeg+incr1
         call giao_der(ngcd,nbls,bl(ibeg),lqijr,lnklr,lqij,lnkl, 1  ,1 ,
     *      bl(ider),bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         ibfij1=ibfij1+incr1
      endif
c----------
c     if(lcase.eq. 3.or.lcase.eq. 6.or.lcase.eq.10.or.lcase.eq.11.or.
c    *   lcase.eq.12.or.lcase.eq.13.or.lcase.eq.15.or.lcase.eq.16) then
      if(lshelij.eq.2 .or. lshelij.eq.3) then
c-   --- for bfij2 ---
c
         ibeg =ibfij2
         incr1=nbls*lqijr*lnklr
         ider=ibeg+incr1
         call giao_der(ngcd,nbls,bl(ibeg),lqijr,lnklr,lqij,lnkl, 1  ,1 ,
     *      bl(ider), bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         ibfij2=ibfij2+incr1
      endif
c----------
c     if(lcase.eq. 4.or.lcase.eq. 7.or.lcase.eq. 8.or.lcase.eq.10.or.
c    *   lcase.eq.12.or.lcase.eq.14.or.lcase.eq.15.or.lcase.eq.16) then
      if(lshelkl.eq.1 .or. lshelkl.eq.3) then
c-   --- for bfkl1 ---
c
         ibeg =ibfkl1
         incr1=nbls*lnijr*lqklr
         ider=ibeg+incr1
         call giao_der(ngcd,nbls,bl(ibeg),lnijr,lqklr,lnij,lqkl, 1  ,1 ,
     *      bl(ider), bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         ibfkl1=ibfkl1+incr1
      endif
c----------
c     if(lcase.eq. 5.or.lcase.eq. 7.or.lcase.eq. 9.or.lcase.eq.11.or.
c    *   lcase.eq.13.or.lcase.eq.14.or.lcase.eq.15.or.lcase.eq.16) then
      if(lshelkl.eq.2 .or. lshelkl.eq.3) then
c-   --- for bfkl2 ---
c
         ibeg =ibfkl2
         incr1=nbls*lnijr*lqklr
         ider=ibeg+incr1
         call giao_der(ngcd,nbls,bl(ibeg),lnijr,lqklr,lnij,lqkl, 1  ,1 ,
     *      bl(ider), bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         ibfkl2=ibfkl2+incr1
      endif
c----------
      if(lshellt.eq.1) go to 100
c----------
c     if(lcase.eq. 6.or.lcase.eq.12.or.lcase.eq.13.or.lcase.eq.16) then
      if(lshelij.eq.3) then
c-   --- for bfij3 (nbls,4,lnklr) ; 4 is for nmr only
c
         ibeg =ibfij3
         incr1=4*nbls*lnklr
         ider=ibeg+incr1
         call giao_der(ngcd,nbls,bl(ibeg),4    ,lnklr,1   ,lnkl, 1  ,1 ,
     *      bl(ider), bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         ibfij3=ibfij3+incr1
      endif
c----------
c     if(lcase.eq. 7.or.lcase.eq.14.or.lcase.eq.15.or.lcase.eq.16) then
      if(lshelkl.eq.3) then
c-   --- for bfkl3 (nbls,lnijr,4) ; 4 is for nmr only
c
         ibeg =ibfkl3
         incr1=4*nbls*lnijr
         ider=ibeg+incr1
         call giao_der(ngcd,nbls,bl(ibeg),lnijr,4    ,lnij,1   , 1  ,1 ,
     *      bl(ider), bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         ibfkl3=ibfkl3+incr1
      endif
c----------
c     if(lcase.eq. 8.or.lcase.eq.12.or.lcase.eq.14.or.lcase.eq.16) then
      if(lcas2(1).eq.1) then
c-   --- for bf2l1 ---
c
         ibeg =ibf2l1
         incr1=nbls*lqijr*lqklr
         ider=ibeg+incr1
         call giao_der(ngcd,nbls,bl(ibeg),lqijr,lqklr,lqij,lqkl, 1  ,1 ,
     *      bl(ider), bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         ibf2l1=ibf2l1+incr1
      endif
c----------
c     if(lcase.eq. 9.or.lcase.eq.13.or.lcase.eq.14.or.lcase.eq.16) then
      if(lcas2(2).eq.1) then
c-   --- for bf2l2 ---
c
         ibeg =ibf2l2
         incr1=nbls*lqijr*lqklr
         ider=ibeg+incr1
         call giao_der(ngcd,nbls,bl(ibeg),lqijr,lqklr,lqij,lqkl, 1  ,1 ,
     *      bl(ider), bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         ibf2l2=ibf2l2+incr1
      endif
c----------
c     if(lcase.eq.10.or.lcase.eq.12.or.lcase.eq.15.or.lcase.eq.16) then
      if(lcas2(3).eq.1) then
c-   --- for bf2l3 ---
c
         ibeg =ibf2l3
         incr1=nbls*lqijr*lqklr
         ider=ibeg+incr1
         call giao_der(ngcd,nbls,bl(ibeg),lqijr,lqklr,lqij,lqkl, 1  ,1 ,
     *      bl(ider), bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         ibf2l3=ibf2l3+incr1
      endif
c----------
c     if(lcase.eq.11.or.lcase.eq.13.or.lcase.eq.15.or.lcase.eq.16) then
      if(lcas2(4).eq.1) then
c-   --- for bf2l4 ---
c
         ibeg =ibf2l4
         incr1=nbls*lqijr*lqklr
         ider=ibeg+incr1
         call giao_der(ngcd,nbls,bl(ibeg),lqijr,lqklr,lqij,lqkl, 1  ,1 ,
     *      bl(ider), bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         ibf2l4=ibf2l4+incr1
      endif
c----------
      if(lshellt.eq.2) go to 100
c----------
c     if(lcase.eq.12.or.lcase.eq.16) then
      if(lcas3(1).eq.1) then
c-   --- for bf3l(nbls,4,lqmx) -first  
c
         ibeg =ibf3l1
         incr1=4*nbls*lqmxr
         ider=ibeg+incr1
         call giao_der(ngcd,nbls,bl(ibeg),4    ,lqmxr,1   ,lqmx, 1  ,1 ,
     *      bl(ider), bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         ibf3l1=ibf3l1+incr1
      endif
c----------
c     if(lcase.eq.13.or.lcase.eq.16) then
      if(lcas3(2).eq.1) then
c-   --- for bf3l(nbls,4,lqmx) - second
c
         ibeg =ibf3l2
         incr1=4*nbls*lqmxr
         ider=ibeg+incr1
         call giao_der(ngcd,nbls,bl(ibeg),4    ,lqmxr,1   ,lqmx, 1  ,1 ,
     *      bl(ider), bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         ibf3l2=ibf3l2+incr1
      endif
c----------
c     if(lcase.eq.14.or.lcase.eq.16) then
      if(lcas3(3).eq.1) then
c-   --- for bf3l(nbls,lqmx,4) - third 
c
         ibeg =ibf3l3
         incr1=4*nbls*lqmxr
         ider=ibeg+incr1
         call giao_der(ngcd,nbls,bl(ibeg),lqmxr,4    ,lqmx,1   , 1  ,1 ,
     *      bl(ider), bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         ibf3l3=ibf3l3+incr1
      endif
c----------
c     if(lcase.eq.15.or.lcase.eq.16) then
      if(lcas3(4).eq.1) then
c-   --- for bf3l(nbls,lqmx,4) - fourth
c
         ibeg =ibf3l4
         incr1=4*nbls*lqmxr
         ider=ibeg+incr1
         call giao_der(ngcd,nbls,bl(ibeg),lqmxr,4    ,lqmx,1   , 1  ,1 ,
     *      bl(ider), bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         ibf3l4=ibf3l4+incr1
      endif
c----------
      if(lcase.eq.16) then
c-   --- for ssss(nbls)  ---
         ibeg =issss
         incr1=nbls*16
         ider=ibeg+incr1
         call giao_der(ngcd,nbls,bl(ibeg),4    ,4    ,1   ,1   , 1  ,1 ,
     *      bl(ider), bl(ixabq),bl(ixcdq),bl(ixyabq),bl(ixycdq))
c-
         issss=issss+incr1
      endif
c
  100 continue
c...................
      call retmem(4)
      end
c=================================================================
c moved into the convert.f file (1998)
c     subroutine conv24x(nbls,npij,npkl,idx1,idx2 ,
c    *                  xab ,xcd, xyab, xycd ,
c    *                  xabq,xcdq,xyabq,xycdq )
c=================================================================
c moved into the convert.f file (1998)
c     subroutine conv24r(nbls,npij,idx1,xab,xabq)
c=================================================================
      subroutine giao_der(ngcd,nbls,buf2,lnijr,lnklr,lnij,lnkl,
     *                    nqij,nqkl, deriv,  xab,xcd, xyab,xycd)
      implicit real*8 (a-h,o-z)
c
#include "texas_lpar.fh"      
c
      dimension buf2(nbls,lnijr,lnklr,ngcd)
      dimension deriv(6,nbls,lnij,lnkl,ngcd)
c--------------------------------------------
c---> dimension typ1x(nbls,lnij,lnkl,ngcd),
c--->*          typ2x(nbls,lnij,lnkl,ngcd)
c---> dimension typ1y(nbls,lnij,lnkl,ngcd),
c--->*          typ2y(nbls,lnij,lnkl,ngcd)
c---> dimension typ1z(nbls,lnij,lnkl,ngcd),
c--->*          typ2z(nbls,lnij,lnkl,ngcd)
c--------------------------------------------
      dimension xab(nbls,3),xcd(nbls,3),xyab(nbls,3),xycd(nbls,3)
c
c-----
c
      do 200 kl=nfu(nqkl)+1,lnkl
      klpx=npxyz(1,kl)
      klpy=npxyz(2,kl)
      klpz=npxyz(3,kl)
      do 200 ij=nfu(nqij)+1,lnij
      ijpx=npxyz(1,ij)
      ijpy=npxyz(2,ij)
      ijpz=npxyz(3,ij)
c
      do 225 iqu=1,ngcd
c
        do 250 ijkl=1,nbls
c---------------
c--x deriv.
        abzy= -xab(ijkl,3)*buf2(ijkl,ijpy,kl,iqu) 
     *        +xab(ijkl,2)*buf2(ijkl,ijpz,kl,iqu)
     *       +xyab(ijkl,1)*buf2(ijkl,ij,kl,iqu)
        cdzy= -xcd(ijkl,3)*buf2(ijkl,ij,klpy,iqu) 
     *        +xcd(ijkl,2)*buf2(ijkl,ij,klpz,iqu)
     *       +xycd(ijkl,1)*buf2(ijkl,ij,kl,iqu)
c
cccccc  typ1x(ijkl,ij,kl,iqu)= abzy+cdzy
cccccc  typ2x(ijkl,ij,kl,iqu)=-abzy+cdzy
        deriv(1,ijkl,ij,kl,iqu)= abzy+cdzy
        deriv(2,ijkl,ij,kl,iqu)=-abzy+cdzy
c
c
c--y deriv.
        abzx= +xab(ijkl,3)*buf2(ijkl,ijpx,kl,iqu) 
     *        -xab(ijkl,1)*buf2(ijkl,ijpz,kl,iqu)
     *       +xyab(ijkl,2)*buf2(ijkl,ij,kl,iqu)
        cdzx= +xcd(ijkl,3)*buf2(ijkl,ij,klpx,iqu) 
     *        -xcd(ijkl,1)*buf2(ijkl,ij,klpz,iqu)
     *       +xycd(ijkl,2)*buf2(ijkl,ij,kl,iqu)
c
cccccc  typ1y(ijkl,ij,kl,iqu)= abzx+cdzx
cccccc  typ2y(ijkl,ij,kl,iqu)=-abzx+cdzx
        deriv(3,ijkl,ij,kl,iqu)= abzx+cdzx
        deriv(4,ijkl,ij,kl,iqu)=-abzx+cdzx
c
c--z deriv.
        abyx= -xab(ijkl,2)*buf2(ijkl,ijpx,kl,iqu) 
     *        +xab(ijkl,1)*buf2(ijkl,ijpy,kl,iqu)
     *       +xyab(ijkl,3)*buf2(ijkl,ij,kl,iqu)
        cdyx= -xcd(ijkl,2)*buf2(ijkl,ij,klpx,iqu) 
     *        +xcd(ijkl,1)*buf2(ijkl,ij,klpy,iqu)
     *       +xycd(ijkl,3)*buf2(ijkl,ij,kl,iqu)
c
cccccc  typ1z(ijkl,ij,kl,iqu)= abyx+cdyx
cccccc  typ2z(ijkl,ij,kl,iqu)=-abyx+cdyx
        deriv(5,ijkl,ij,kl,iqu)= abyx+cdyx
        deriv(6,ijkl,ij,kl,iqu)=-abyx+cdyx
c
  250   continue
  225 continue
  200 continue
c
      end
c===================================================================
c===========GRADIENT INTEGRAL DERIVATIVES ROUTINES==================
c
c     It is called from Calcint2 when WHERE='forc'
c
      subroutine force_der(bl,nbls,lnijr,lnklr,npij,ngcd,idx1,ixab)
c
      implicit real*8 (a-h,o-z)
c
#include "texas_lpar.fh"      
      COMMON/SHELL/LSHELLT,LSHELIJ,LSHELKL,LHELP,LCAS2(4),LCAS3(4)
      common /lcases/ lcase
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
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
c only for first & second derivatives (for use in amshift):
cccc  common /memor4b/ibuf0
      common /memor4b/ider0,ider1,ider2
c
c dimensions for assembling :
      common /dimasse/ lqijr,lqklr,lqmxr,lij3,lkl3,l3l,lsss
c
      dimension bl(*)
c----------------------------------------------------------
      lqij=nfu(nqij+1)
      lqkl=nfu(nqkl+1)
      lqmx=lqij
      if(lqkl.gt.lqij) lqmx=lqkl
c----------------------------------------------------------
      ndim=4     ! dimension for buf2(ndim,*) used in first_der
c----------------------------------------------------------
      call getmem(3*nbls,ixabq)
      call conv24r(nbls,npij,bl(idx1),bl(ixab),bl(ixabq))
c----------------------------------------------------------
c contracted integrals with higher and lower ang.mom. are in buf now,
c not in buf2. In buf2 derivatives will be placed after call to 
c first_der
c
      ibeg =ibuf
      ider1=ibuf2
      incr9=9*ngcd*nbls*lnij*lnkl
      ider0=ider1+incr9
c
      call first_der(ngcd,nbls,bl(ibeg),ndim,
     *               lnijr,lnklr,lnij,lnkl,nqij,nqkl,
     *               bl(ider1),bl(ider0),bl(ixabq))
c
cccc  ibuf0=ider0
c
c!!!!!!!   poprawic incre and ider0 ponizej !!!!!!!!!
c------
      if(lshellt.eq.0) go to 100
c------
c
c     if(lcase.eq. 2.or.lcase.eq. 6.or.lcase.eq. 8.or.lcase.eq. 9.or.
c    *   lcase.eq.12.or.lcase.eq.13.or.lcase.eq.14.or.lcase.eq.16) then
      if(lshelij.eq.1 .or. lshelij.eq.3) then
c-   --- for bfij1 ---
c
         ibeg =ibfij1
         incr4=4*nbls*lqijr*lnklr
         incr9=9*nbls*lqij*lnkl
         ider1=ibeg+incr4
         ider0=ider1+incr9
         call first_der(ngcd,nbls,bl(ibeg),ndim,
     *                  lqijr,lnklr,lqij,lnkl, 1 ,1 ,
     *                  bl(ider1),bl(ider0),bl(ixabq))
c-
         ibfij1=ider1
         ibfij1_0=ider0
      endif
c----------
c     if(lcase.eq. 3.or.lcase.eq. 6.or.lcase.eq.10.or.lcase.eq.11.or.
c    *   lcase.eq.12.or.lcase.eq.13.or.lcase.eq.15.or.lcase.eq.16) then
      if(lshelij.eq.2 .or. lshelij.eq.3) then
c-   --- for bfij2 ---
c
         ibeg =ibfij2
         incr4=4*nbls*lqijr*lnklr
         incr9=9*nbls*lqij*lnkl
         ider1=ibeg+incr4
         ider0=ider1+incr9
         call first_der(ngcd,nbls,bl(ibeg),ndim,
     *                  lqijr,lnklr,lqij,lnkl, 1 ,1 ,
     *                  bl(ider1),bl(ider0),bl(ixabq))
c-
         ibfij2=ider1
         ibfij2_0=ider0
      endif
c----------
c     if(lcase.eq. 4.or.lcase.eq. 7.or.lcase.eq. 8.or.lcase.eq.10.or.
c    *   lcase.eq.12.or.lcase.eq.14.or.lcase.eq.15.or.lcase.eq.16) then
      if(lshelkl.eq.1 .or. lshelkl.eq.3) then
c-   --- for bfkl1 ---
c
         ibeg =ibfkl1
         incr4=4*nbls*lnijr*lqklr
         incr9=9*nbls*lnij*lqkl
         ider1=ibeg+incr4
         ider0=ider1+incr9
         call first_der(ngcd,nbls,bl(ibeg),ndim,
     *                  lnijr,lqklr,lnij,lqkl, 1 ,1 ,
     *                  bl(ider1),bl(ider0),bl(ixabq))
c-
         ibfkl1=ider1
         ibfkl1_0=ider0
      endif
c----------
c     if(lcase.eq. 5.or.lcase.eq. 7.or.lcase.eq. 9.or.lcase.eq.11.or.
c    *   lcase.eq.13.or.lcase.eq.14.or.lcase.eq.15.or.lcase.eq.16) then
      if(lshelkl.eq.2 .or. lshelkl.eq.3) then
c-   --- for bfkl2 ---
c
         ibeg =ibfkl2
         incr4=4*nbls*lnijr*lqklr
         incr9=9*nbls*lnij*lqkl
         ider1=ibeg+incr4
         ider0=ider1+incr9
         call first_der(ngcd,nbls,bl(ibeg),ndim,
     *                  lnijr,lqklr,lnij,lqkl, 1 ,1 ,
     *                  bl(ider1),bl(ider0),bl(ixabq))
c-
         ibfkl2=ider1
         ibfkl2_0=ider0
      endif
c----------
      if(lshellt.eq.1) go to 100
c----------
c     if(lcase.eq. 6.or.lcase.eq.12.or.lcase.eq.13.or.lcase.eq.16) then
      if(lshelij.eq.3) then
c-   --- for bfij3 (nbls,4,lnklr) ; 4 is for nmr only
c
         ibeg =ibfij3
         incr4=4*4*nbls*lnklr
         incr9=9*4*nbls*lnkl 
         ider1=ibeg+incr4
         ider0=ider1+incr9
         call first_der(ngcd,nbls,bl(ibeg),ndim,
     *                  4    ,lnklr,1   ,lnkl, 1 ,1 ,
     *                  bl(ider1),bl(ider0),bl(ixabq))
c-
         ibfij3=ider1
         ibfij3_0=ider0
      endif
c----------
c     if(lcase.eq. 7.or.lcase.eq.14.or.lcase.eq.15.or.lcase.eq.16) then
      if(lshelkl.eq.3) then
c-   --- for bfkl3 (nbls,lnijr,4) ; 4 is for nmr only
c
         ibeg =ibfkl3
         incr4=4*4*nbls*lnijr
         incr9=9*4*nbls*lnij
         ider1=ibeg+incr4
         ider0=ider1+incr9
         call first_der(ngcd,nbls,bl(ibeg),ndim,
     *                  lnijr,4    ,lnij,1   , 1 ,1 ,
     *                  bl(ider1),bl(ider0),bl(ixabq))
c-
         ibfkl3=ider1
         ibfkl3_0=ider0
      endif
c----------
c     if(lcase.eq. 8.or.lcase.eq.12.or.lcase.eq.14.or.lcase.eq.16) then
      if(lcas2(1).eq.1) then
c-   --- for bf2l1 ---
c
         ibeg =ibf2l1
         incr4=4*nbls*lqijr*lqklr
         incr9=9*nbls*lqij *lqkl 
         ider1=ibeg+incr4
         ider0=ider1+incr9
         call first_der(ngcd,nbls,bl(ibeg),ndim,
     *                  lqijr,lqklr,lqij,lqkl, 1 ,1 ,
     *                  bl(ider1),bl(ider0),bl(ixabq))
c-
         ibf2l1=ider1
         ibf2l1_0=ider0
      endif
c----------
c     if(lcase.eq. 9.or.lcase.eq.13.or.lcase.eq.14.or.lcase.eq.16) then
      if(lcas2(2).eq.1) then
c-   --- for bf2l2 ---
c
         ibeg =ibf2l2
         incr4=4*nbls*lqijr*lqklr
         incr9=9*nbls*lqij *lqkl 
         ider1=ibeg+incr4
         ider0=ider1+incr9
         call first_der(ngcd,nbls,bl(ibeg),ndim,
     *                  lqijr,lqklr,lqij,lqkl, 1 ,1 ,
     *                  bl(ider1),bl(ider0),bl(ixabq))
c-
         ibf2l2=ider1
         ibf2l2_0=ider0
      endif
c----------
c     if(lcase.eq.10.or.lcase.eq.12.or.lcase.eq.15.or.lcase.eq.16) then
      if(lcas2(3).eq.1) then
c-   --- for bf2l3 ---
c
         ibeg =ibf2l3
         incr4=4*nbls*lqijr*lqklr
         incr9=9*nbls*lqij *lqkl 
         ider1=ibeg+incr4
         ider0=ider1+incr9
         call first_der(ngcd,nbls,bl(ibeg),ndim,
     *                  lqijr,lqklr,lqij,lqkl, 1 ,1 ,
     *                  bl(ider1),bl(ider0),bl(ixabq))
c-
         ibf2l3=ider1
         ibf2l3_0=ider0
      endif
c----------
c     if(lcase.eq.11.or.lcase.eq.13.or.lcase.eq.15.or.lcase.eq.16) then
      if(lcas2(4).eq.1) then
c-   --- for bf2l4 ---
c
         ibeg =ibf2l4
         incr4=4*nbls*lqijr*lqklr
         incr9=9*nbls*lqij *lqkl 
         ider1=ibeg+incr4
         ider0=ider1+incr9
         call first_der(ngcd,nbls,bl(ibeg),ndim,
     *                  lqijr,lqklr,lqij,lqkl, 1 ,1 ,
     *                  bl(ider1),bl(ider0),bl(ixabq))
c-
         ibf2l4=ider1
         ibf2l4_0=ider0
      endif
c----------
      if(lshellt.eq.2) go to 100
c----------
c     if(lcase.eq.12.or.lcase.eq.16) then
      if(lcas3(1).eq.1) then
c-   --- for bf3l(nbls,4,lqmx) -first  
c
         ibeg =ibf3l1
         incr4=4*4*nbls*lqmxr
         incr9=9*4*nbls*lqmx 
         ider1=ibeg+incr4
         ider0=ider1+incr9
         call first_der(ngcd,nbls,bl(ibeg),ndim,
     *                  4    ,lqmxr,1   ,lqmx, 1 ,1 ,
     *                  bl(ider1),bl(ider0),bl(ixabq))
c-
         ibf3l1=ider1
         ibf3l1_0=ider0
      endif
c----------
c     if(lcase.eq.13.or.lcase.eq.16) then
      if(lcas3(2).eq.1) then
c-   --- for bf3l(nbls,4,lqmx) - second
c
         ibeg =ibf3l2
         incr4=4*4*nbls*lqmxr
         incr9=9*4*nbls*lqmx
         ider1=ibeg+incr4
         ider0=ider1+incr9
         call first_der(ngcd,nbls,bl(ibeg),ndim,
     *                  4    ,lqmxr,1   ,lqmx, 1 ,1 ,
     *                  bl(ider1),bl(ider0),bl(ixabq))
c-
         ibf3l2=ider1
         ibf3l2_0=ider0
      endif
c----------
c     if(lcase.eq.14.or.lcase.eq.16) then
      if(lcas3(3).eq.1) then
c-   --- for bf3l(nbls,lqmx,4) - third 
c
         ibeg =ibf3l3
         incr4=4*4*nbls*lqmxr
         incr9=9*4*nbls*lqmx
         ider1=ibeg+incr4
         ider0=ider1+incr9
         call first_der(ngcd,nbls,bl(ibeg),ndim,
     *                  lqmxr,4    ,lqmx,1   , 1 ,1 ,
     *                  bl(ider1),bl(ider0),bl(ixabq))
c-
         ibf3l3=ider1
         ibf3l3_0=ider0
      endif
c----------
c     if(lcase.eq.15.or.lcase.eq.16) then
      if(lcas3(4).eq.1) then
c-   --- for bf3l(nbls,lqmx,4) - fourth
c
         ibeg =ibf3l4
         incr4=4*4*nbls*lqmxr
         incr9=9*4*nbls*lqmx
         ider1=ibeg+incr4
         ider0=ider1+incr9
         call first_der(ngcd,nbls,bl(ibeg),ndim,
     *                  lqmxr,4    ,lqmx,1   , 1 ,1 ,
     *                  bl(ider1),bl(ider0),bl(ixabq))
c-
         ibf3l4=ider1
         ibf3l4_0=ider0
      endif
c----------
      if(lcase.eq.16) then
c-   --- for ssss(nbls)  ---
         ibeg =issss
         incr4=4*nbls*16
         incr4=9*nbls*16
         ider1=ibeg+incr4
         ider0=ider1+incr9
         call first_der(ngcd,nbls,bl(ibeg),ndim,
     *                  4    ,4    ,1   ,1   , 1 ,1 ,
     *                  bl(ider1),bl(ider0),bl(ixabq))
c-
         issss=ider1
         issss_0=ider0
      endif
c
  100 continue
c...................
      call retmem(1)
      end
c=================================================================
      subroutine first_der(ngcd,nbls,buf2,ndim,lnijr,lnklr,lnij,lnkl,
     *                     nqij,nqkl,deriv,der00,xab)
      implicit real*8 (a-h,o-z)
c
c this subroutine is also called when second-derivatives are calculated
c because first-derivatives are needed in the shifting procedure for
c second-der. That's why dimension for buf2(ndim,*,*,*,*) has ndim=4
c for first- and ndim=10 for second-derivatives.
c
#include "texas_lpar.fh"      
c
cccc  dimension buf2(4,nbls,lnijr,lnklr,ngcd) OR buf2(10,etc.)
c2002 dimension buf2(ndim,nbls,lnijr,lnklr,ngcd)
      dimension buf2(nbls,lnijr,lnklr,ngcd,ndim)
      dimension deriv(9,nbls,lnij,lnkl,ngcd),der00(nbls,lnij,lnkl,ngcd)
      dimension xab(nbls,3)
c-------------------------------------------------------------------
c This routine calculates UNFINISHED gradient derivatives 
c of the type 
c
c d/dAx (a+b,s|c+d,s)=2a*(a+b+1x,s|c+d,s)-n_ab_x*(a+b-1x,s|c+d,s)
c
c d/dBx (a+b,s|c+d,s)=2b*(a+b,s+1x|c+d,s)-  0*(a+b,s-1x|c+d,s)
c                    =2b[(a+b+1x,s|c+d,s)+(Ax-Bx)*(a+b,s|c+d,s)].
c
c d/dCx (a+b,s|c+d,s)=2c*(a+b,s|c+d+1x,s)-n_cd_x*(a+b,s|c+d-1x,s)
c
c-------------------------------------------------------------------
c INPUT buf2(1,nbls,lnijr,lnklr,ngcd) - ordinary 2-el.integ.
c INPUT buf2(2,nbls,lnijr,lnklr,ngcd) - 2-el.integ. resc by 2*exp_a
c INPUT buf2(3,nbls,lnijr,lnklr,ngcd) - 2-el.integ. resc by 2*exp_b
c INPUT buf2(4,nbls,lnijr,lnklr,ngcd) - 2-el.integ. resc by 2*exp_c
c
c OUTPUT : 
c
c--->           derAX=deriv(1,nbls,lnij,lnkl,ngcd),
c--->           derBX=deriv(2,nbls,lnij,lnkl,ngcd)
c--->           derCX=deriv(3,nbls,lnij,lnkl,ngcd),
c
c--->           derAY=deriv(4,nbls,lnij,lnkl,ngcd),
c--->           derBY=deriv(5,nbls,lnij,lnkl,ngcd),
c--->           derCY=deriv(6,nbls,lnij,lnkl,ngcd),
c
c--->           derAZ=deriv(7,nbls,lnij,lnkl,ngcd),
c--->           derBZ=deriv(8,nbls,lnij,lnkl,ngcd),
c--->           derCZ=deriv(9,nbls,lnij,lnkl,ngcd),
c
c--->           der00(nbls,lnij,lnkl,ngcd) - ordinary integr.
c               needed in shifting for d/dA,d/dB (1->2) and d/dC (3->4)
c-------------------------------------------------------------------
c
      do 200 kl=nfu(nqkl)+1,lnkl
      n_cd_x=nia(1,kl)
      n_cd_y=nia(2,kl)
      n_cd_z=nia(3,kl)
      klpx=npxyz(1,kl)
      klpy=npxyz(2,kl)
      klpz=npxyz(3,kl)
c
      klmx=0
      klmy=0
      klmz=0
      if(n_cd_x.gt.0) klmx=nmxyz(1,kl)
      if(n_cd_y.gt.0) klmy=nmxyz(2,kl)
      if(n_cd_z.gt.0) klmz=nmxyz(3,kl)
c
      do 200 ij=nfu(nqij)+1,lnij
      n_ab_x=nia(1,ij)
      n_ab_y=nia(2,ij)
      n_ab_z=nia(3,ij)
      ijpx=npxyz(1,ij)
      ijpy=npxyz(2,ij)
      ijpz=npxyz(3,ij)
      ijmx=0
      ijmy=0
      ijmz=0
      if(n_ab_x.gt.0) ijmx=nmxyz(1,ij)
      if(n_ab_y.gt.0) ijmy=nmxyz(2,ij)
      if(n_ab_z.gt.0) ijmz=nmxyz(3,ij)
c
      do 225 iqu=1,ngcd
        do 250 ijkl=1,nbls
c---------------
        two_a_x=buf2(ijkl,ijpx,kl,iqu,2) 
        two_a_y=buf2(ijkl,ijpy,kl,iqu,2) 
        two_a_z=buf2(ijkl,ijpz,kl,iqu,2) 
c
        two_b_0=buf2(ijkl,ij,kl,iqu,3) 
        two_b_x=buf2(ijkl,ijpx,kl,iqu,3) 
        two_b_y=buf2(ijkl,ijpy,kl,iqu,3) 
        two_b_z=buf2(ijkl,ijpz,kl,iqu,3)
c
        two_c_x=buf2(ijkl,ij,klpx,iqu,4) 
        two_c_y=buf2(ijkl,ij,klpy,iqu,4) 
        two_c_z=buf2(ijkl,ij,klpz,iqu,4) 
c
c---------------
        der00(ijkl,ij,kl,iqu)=buf2(ijkl,ij,kl,iqu,1)
c---------------
c--X deriv.
c
        if(n_ab_x.gt.0) then
           x_n_ab=n_ab_x*buf2(ijkl,ijmx,kl,iqu,1)
           deriv(1,ijkl,ij,kl,iqu)=two_a_x - x_n_ab
        else
           deriv(1,ijkl,ij,kl,iqu)=two_a_x
        endif
        deriv(2,ijkl,ij,kl,iqu)=two_b_x + xab(ijkl,1)*two_b_0
        if(n_cd_x.gt.0) then
           x_n_cd=n_cd_x*buf2(ijkl,ij,klmx,iqu,1)
           deriv(3,ijkl,ij,kl,iqu)=two_c_x - x_n_cd
        else
           deriv(3,ijkl,ij,kl,iqu)=two_c_x
        endif
c---------------
c--Y deriv.
c
        if(n_ab_y.gt.0) then
           y_n_ab=n_ab_y*buf2(ijkl,ijmy,kl,iqu,1)
           deriv(4,ijkl,ij,kl,iqu)=two_a_y - y_n_ab
        else
           deriv(4,ijkl,ij,kl,iqu)=two_a_y
        endif
        deriv(5,ijkl,ij,kl,iqu)=two_b_y + xab(ijkl,2)*two_b_0
        if(n_cd_y.gt.0) then
           y_n_cd=n_cd_y*buf2(ijkl,ij,klmy,iqu,1)
           deriv(6,ijkl,ij,kl,iqu)=two_c_y - y_n_cd
        else
           deriv(6,ijkl,ij,kl,iqu)=two_c_y
        endif
c---------------
c--Z deriv.
c
        if(n_ab_z.gt.0) then
           z_n_ab=n_ab_z*buf2(ijkl,ijmz,kl,iqu,1)
           deriv(7,ijkl,ij,kl,iqu)=two_a_z - z_n_ab
        else
           deriv(7,ijkl,ij,kl,iqu)=two_a_z
        endif
        deriv(8,ijkl,ij,kl,iqu)=two_b_z + xab(ijkl,3)*two_b_0
        if(n_cd_z.gt.0) then
           z_n_cd=n_cd_z*buf2(ijkl,ij,klmz,iqu,1)
           deriv(9,ijkl,ij,kl,iqu)=two_c_z - z_n_cd
        else
           deriv(9,ijkl,ij,kl,iqu)=two_c_z
        endif
c---------------
c
  250   continue
  225 continue
  200 continue
c
      end
c=================================================================
c
c===========HESSIAN  INTEGRAL DERIVATIVES ROUTINES==================
c
c     It is called from Calcint2 when WHERE='hess'
c
      subroutine hessian_der(bl,nbls,lnijr,lnklr,npij,ngcd,idx1,ixab)
c
      implicit real*8 (a-h,o-z)
c
#include "texas_lpar.fh"      
      COMMON/SHELL/LSHELLT,LSHELIJ,LSHELKL,LHELP,LCAS2(4),LCAS3(4)
      common /lcases/ lcase
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
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
c only for first & second derivatives (for use in amshift):
cccc  common /memor4b/ibuf0
      common /memor4b/ider0,ider1,ider2
c
c dimensions for assembling :
      common /dimasse/ lqijr,lqklr,lqmxr,lij3,lkl3,l3l,lsss
c
      dimension bl(*)
c----------------------------------------------------------
      lqij=nfu(nqij+1)
      lqkl=nfu(nqkl+1)
      lqmx=lqij
        if(lqkl.gt.lqij) lqmx=lqkl
c----------------------------------------------------------
      ndim=10   ! dimension for buf2(ndim,*) used in first_der
c----------------------------------------------------------
      call getmem(3*nbls,ixabq)
      call conv24r(nbls,npij,bl(idx1),bl(ixab),bl(ixabq))
c----------------------------------------------------------
      ibeg =ibuf
      ider2=ibuf2
      incr45=45*ngcd*nbls*lnij*lnkl
      ider1=ider2+incr45
      incr9 = 9*ngcd*nbls*lnij*lnkl
      ider0=ider1+incr9
c
      call first_der(ngcd,nbls,bl(ibeg),ndim,
     *               lnijr,lnklr,lnij,lnkl,nqij,nqkl,
     *               bl(ider1),bl(ider0),bl(ixabq))
c
      call secnd_der(ngcd,nbls,bl(ibeg),lnijr,lnklr,lnij,lnkl,nqij,nqkl,
     *               bl(ider2),bl(ixabq))
c
c
c!!!!!!!   NOTHING IS READY FOR L-shells    !!!!!!!!!
c------
      if(lshellt.eq.0) go to 100
c------
c
c     if(lcase.eq. 2.or.lcase.eq. 6.or.lcase.eq. 8.or.lcase.eq. 9.or.
c    *   lcase.eq.12.or.lcase.eq.13.or.lcase.eq.14.or.lcase.eq.16) then
      if(lshelij.eq.1 .or. lshelij.eq.3) then
c-   --- for bfij1 ---
c
         ibeg =ibfij1
         incre=nbls*lqijr*lnklr
         incr4=4*incre
         incr9=9*incre
         ider=ibeg+incr4
         ider0=ider+incr9
         call secnd_der(ngcd,nbls,bl(ibeg),lqijr,lnklr,lqij,lnkl, 1 ,1 ,
     *                  bl(ider),bl(ixabq))
c-
         ibfij1=ibfij1+incr4
      endif
c----------
c     if(lcase.eq. 3.or.lcase.eq. 6.or.lcase.eq.10.or.lcase.eq.11.or.
c    *   lcase.eq.12.or.lcase.eq.13.or.lcase.eq.15.or.lcase.eq.16) then
      if(lshelij.eq.2 .or. lshelij.eq.3) then
c-   --- for bfij2 ---
c
         ibeg =ibfij2
         incr1=4*nbls*lqijr*lnklr
         ider=ibeg+incr1
         call secnd_der(ngcd,nbls,bl(ibeg),lqijr,lnklr,lqij,lnkl, 1 ,1 ,
     *                  bl(ider),bl(ixabq))
c-
         ibfij2=ibfij2+incr1
      endif
c----------
c     if(lcase.eq. 4.or.lcase.eq. 7.or.lcase.eq. 8.or.lcase.eq.10.or.
c    *   lcase.eq.12.or.lcase.eq.14.or.lcase.eq.15.or.lcase.eq.16) then
      if(lshelkl.eq.1 .or. lshelkl.eq.3) then
c-   --- for bfkl1 ---
c
         ibeg =ibfkl1
         incr1=4*nbls*lnijr*lqklr
         ider=ibeg+incr1
         call secnd_der(ngcd,nbls,bl(ibeg),lnijr,lqklr,lnij,lqkl, 1 ,1 ,
     *                  bl(ider),bl(ixabq))
c-
         ibfkl1=ibfkl1+incr1
      endif
c----------
c     if(lcase.eq. 5.or.lcase.eq. 7.or.lcase.eq. 9.or.lcase.eq.11.or.
c    *   lcase.eq.13.or.lcase.eq.14.or.lcase.eq.15.or.lcase.eq.16) then
      if(lshelkl.eq.2 .or. lshelkl.eq.3) then
c-   --- for bfkl2 ---
c
         ibeg =ibfkl2
         incr1=4*nbls*lnijr*lqklr
         ider=ibeg+incr1
         call secnd_der(ngcd,nbls,bl(ibeg),lnijr,lqklr,lnij,lqkl, 1 ,1 ,
     *                  bl(ider),bl(ixabq))
c-
         ibfkl2=ibfkl2+incr1
      endif
c----------
      if(lshellt.eq.1) go to 100
c----------
c     if(lcase.eq. 6.or.lcase.eq.12.or.lcase.eq.13.or.lcase.eq.16) then
      if(lshelij.eq.3) then
c-   --- for bfij3 (nbls,4,lnklr) ; 4 is for nmr only
c
         ibeg =ibfij3
         incr1=4*4*nbls*lnklr
         ider=ibeg+incr1
         call secnd_der(ngcd,nbls,bl(ibeg),4    ,lnklr,1   ,lnkl, 1 ,1 ,
     *                  bl(ider),bl(ixabq))
c-
         ibfij3=ibfij3+incr1
      endif
c----------
c     if(lcase.eq. 7.or.lcase.eq.14.or.lcase.eq.15.or.lcase.eq.16) then
      if(lshelkl.eq.3) then
c-   --- for bfkl3 (nbls,lnijr,4) ; 4 is for nmr only
c
         ibeg =ibfkl3
         incr1=4*4*nbls*lnijr
         ider=ibeg+incr1
         call secnd_der(ngcd,nbls,bl(ibeg),lnijr,4    ,lnij,1   , 1 ,1 ,
     *                  bl(ider),bl(ixabq))
c-
         ibfkl3=ibfkl3+incr1
      endif
c----------
c     if(lcase.eq. 8.or.lcase.eq.12.or.lcase.eq.14.or.lcase.eq.16) then
      if(lcas2(1).eq.1) then
c-   --- for bf2l1 ---
c
         ibeg =ibf2l1
         incr1=4*nbls*lqijr*lqklr
         ider=ibeg+incr1
         call secnd_der(ngcd,nbls,bl(ibeg),lqijr,lqklr,lqij,lqkl, 1 ,1 ,
     *                  bl(ider),bl(ixabq))
c-
         ibf2l1=ibf2l1+incr1
      endif
c----------
c     if(lcase.eq. 9.or.lcase.eq.13.or.lcase.eq.14.or.lcase.eq.16) then
      if(lcas2(2).eq.1) then
c-   --- for bf2l2 ---
c
         ibeg =ibf2l2
         incr1=4*nbls*lqijr*lqklr
         ider=ibeg+incr1
         call secnd_der(ngcd,nbls,bl(ibeg),lqijr,lqklr,lqij,lqkl, 1 ,1 ,
     *                  bl(ider),bl(ixabq))
c-
         ibf2l2=ibf2l2+incr1
      endif
c----------
c     if(lcase.eq.10.or.lcase.eq.12.or.lcase.eq.15.or.lcase.eq.16) then
      if(lcas2(3).eq.1) then
c-   --- for bf2l3 ---
c
         ibeg =ibf2l3
         incr1=4*nbls*lqijr*lqklr
         ider=ibeg+incr1
         call secnd_der(ngcd,nbls,bl(ibeg),lqijr,lqklr,lqij,lqkl, 1 ,1 ,
     *                  bl(ider),bl(ixabq))
c-
         ibf2l3=ibf2l3+incr1
      endif
c----------
c     if(lcase.eq.11.or.lcase.eq.13.or.lcase.eq.15.or.lcase.eq.16) then
      if(lcas2(4).eq.1) then
c-   --- for bf2l4 ---
c
         ibeg =ibf2l4
         incr1=4*nbls*lqijr*lqklr
         ider=ibeg+incr1
         call secnd_der(ngcd,nbls,bl(ibeg),lqijr,lqklr,lqij,lqkl, 1 ,1 ,
     *                  bl(ider),bl(ixabq))
c-
         ibf2l4=ibf2l4+incr1
      endif
c----------
      if(lshellt.eq.2) go to 100
c----------
c     if(lcase.eq.12.or.lcase.eq.16) then
      if(lcas3(1).eq.1) then
c-   --- for bf3l(nbls,4,lqmx) -first  
c
         ibeg =ibf3l1
         incr1=4*4*nbls*lqmxr
         ider=ibeg+incr1
         call secnd_der(ngcd,nbls,bl(ibeg),4    ,lqmxr,1   ,lqmx, 1 ,1 ,
     *                  bl(ider),bl(ixabq))
c-
         ibf3l1=ibf3l1+incr1
      endif
c----------
c     if(lcase.eq.13.or.lcase.eq.16) then
      if(lcas3(2).eq.1) then
c-   --- for bf3l(nbls,4,lqmx) - second
c
         ibeg =ibf3l2
         incr1=4*4*nbls*lqmxr
         ider=ibeg+incr1
         call secnd_der(ngcd,nbls,bl(ibeg),4    ,lqmxr,1   ,lqmx, 1 ,1 ,
     *                  bl(ider),bl(ixabq))
c-
         ibf3l2=ibf3l2+incr1
      endif
c----------
c     if(lcase.eq.14.or.lcase.eq.16) then
      if(lcas3(3).eq.1) then
c-   --- for bf3l(nbls,lqmx,4) - third 
c
         ibeg =ibf3l3
         incr1=4*4*nbls*lqmxr
         ider=ibeg+incr1
         call secnd_der(ngcd,nbls,bl(ibeg),lqmxr,4    ,lqmx,1   , 1 ,1 ,
     *                  bl(ider),bl(ixabq))
c-
         ibf3l3=ibf3l3+incr1
      endif
c----------
c     if(lcase.eq.15.or.lcase.eq.16) then
      if(lcas3(4).eq.1) then
c-   --- for bf3l(nbls,lqmx,4) - fourth
c
         ibeg =ibf3l4
         incr1=4*4*nbls*lqmxr
         ider=ibeg+incr1
         call secnd_der(ngcd,nbls,bl(ibeg),lqmxr,4    ,lqmx,1   , 1 ,1 ,
     *                  bl(ider),bl(ixabq))
c-
         ibf3l4=ibf3l4+incr1
      endif
c----------
      if(lcase.eq.16) then
c-   --- for ssss(nbls)  ---
         ibeg =issss
         incr1=4*nbls*16
         ider=ibeg+incr1
         call secnd_der(ngcd,nbls,bl(ibeg),4    ,4    ,1   ,1   , 1 ,1 ,
     *                  bl(ider),bl(ixabq))
c-
         issss=issss+incr1
      endif
c
  100 continue
c...................
      call retmem(1)
      end
c=================================================================
      subroutine secnd_der(ngcd,nbls,buf2,lnijr,lnklr,lnij,lnkl,
     *                     nqij,nqkl,der2,xab)
      implicit real*8 (a-h,o-z)
c
#include "texas_lpar.fh"      
c
c2002 dimension buf2(10,nbls,lnijr,lnklr,ngcd)
      dimension buf2(nbls,lnijr,lnklr,ngcd,10)
cccc  dimension der0(nbls,lnij,lnkl,ngcd)        these two are constracted
cccc  dimension der1(9,nbls,lnij,lnkl,ngcd)      when first_der is called
      dimension der2(45,nbls,lnij,lnkl,ngcd)
      dimension xab(nbls,3)
c----------------------------------------------------------------------
c This routine calculates UNFINISHED hessian derivatives (second order)
c of the type 
c----------------------------------------------------------------------
c d2/dAxdAy (a+b,s|c+d,s)=    (2a)*(2a)*(a+b+1x+1y,s|c+d,s)
c                          -n_ab_x*(2a)*(a+b-1x+1y,s|c+d,s)
c                          -n_ab_y*(2a)*(a+b+1x-1y,s|c+d,s)
c                        +n_ab_x*n_ab_y*(a+b-1x-1y,s|c+d,s)
c
c d2/dAxdBy (a+b,s|c+d,s)=
c          (2a)*(2b)*[ (a+b+1x+1y,s|c+d,s) + (Ay-By)*(a+b+1x,s|c+d,s) ]
c            -n_ab_x*[ (a+b-1x+1y,s|c+d,s) + (Ay-By)*(a+b-1x,s|c+d,s) ]
c
c d2/dAxdCy (a+b,s|c+d,s)=    (2a)*(2c)*(a+b+1x,s|c+d+1y,s)
c                          -n_ab_x*(2c)*(a+b-1x,s|c+d+1y,s)
c                          -n_cd_y*(2a)*(a+b+1x,s|c+d-1y,s)
c                        +n_ab_x*n_cd_y*(a+b-1x,s|c+d-1y,s)
c
c d2/dAxdDy (a+b,s|c+d,s)= - [ d2/dAxdAy (a+b,s|c+d,s)
c                             +d2/dAxdBy (a+b,s|c+d,s)
c                             +d2/dAxdCy (a+b,s|c+d,s) ]
c
c----------------------------------------------------------------------
c d2/dBxdAy (a+b,s|c+d,s)= d2/dAydBx (a+b,s|c+d,s) see above
c
c d2/dBxdBy (a+b,s|c+d,s)=          (2b)*(2b)*(a+b+1x+1y,s|c+d,s)
c                          +(Ay-By)*(2b)*(2b)*(a+b+1x,s|c+d,s)
c                          +(Ax-Bx)*(2b)*     (a+b+1y,s|c+d,s)
c                          +(Ax-Bx)*(Ay-By)*  (a+b,s|c+d,s)
c                          +d(Ax-Bx)/dBy*(2b)*(a+b,s|c+d,s)
c
c d2/dBxdCy (a+b,s|c+d,s)=          (2b)*(2c)*(a+b+1x,s|c+d+1y,s)
c                           -n_cd_y*(2b)*     (a+b+1x,s|c+d-1y,s)
c                          +(Ax-Bx)*(2b)*(2c)*(a+b,s|c+d+1y,s)
c                        -(Ax-Bx)*n_cd_y*(2b)*(a+b,s|c+d-1y,s)
c
c d2/dBxdDy (a+b,s|c+d,s)= - [ d2/dBxdAy (a+b,s|c+d,s)
c                             +d2/dBxdBy (a+b,s|c+d,s)
c                             +d2/dBxdCy (a+b,s|c+d,s) ]
c
c----------------------------------------------------------------------
c d2/dCxdAy (a+b,s|c+d,s)= d2/dAydCx (a+b,s|c+d,s)  see above
c d2/dCxdBy (a+b,s|c+d,s)= d2/dBydCx (a+b,s|c+d,s)  see above
c
c d2/dCxdCy (a+b,s|c+d,s)=    (2c)*(2c)*(a+b,s|c+d+1x+1y,s)
c                          -n_cd_x*(2c)*(a+b,s|c+d-1x+1y,s)
c                          -n_cd_y*(2c)*(a+b,s|c+d+1x-1y,s)
c                        +n_cd_x*n_cd_y*(a+b,s|c+d-1x-1y,s)
c
c d2/dCxdDy (a+b,s|c+d,s)= - [ d2/dCxdAy (a+b,s|c+d,s)
c                             +d2/dCxdBy (a+b,s|c+d,s)
c                             +d2/dCxdCy (a+b,s|c+d,s) ]
c
c----------------------------------------------------------------------
c d2/dDxdAy (a+b,s|c+d,s)= d2/dAydDx (a+b,s|c+d,s)  see above
c d2/dDxdBy (a+b,s|c+d,s)= d2/dBydDx (a+b,s|c+d,s)  see above
c d2/dDxdCy (a+b,s|c+d,s)= d2/dCydDx (a+b,s|c+d,s)  see above
c
c d2/dDxdDy (a+b,s|c+d,s)= - [ d2/dDxdAy (a+b,s|c+d,s)
c                             +d2/dDxdBy (a+b,s|c+d,s)
c                             +d2/dDxdCy (a+b,s|c+d,s) ]
c
c----------------------------------------------------------------------
c Order for derivatives:       AA AB AC AD
c                                 BB BC BD
c                                    CC CD
c                                       DD
c
c       Block AA :         Block AB :          Block AC :     Block AD :
c  1 -  d2/dAxdAx     7 -  d2/dAxdBx     16 -  d2/dAxdCx      trans.inv.
c  2 -  d2/dAxdAy     8 -  d2/dAxdBy     17 -  d2/dAxdCy
c  3 -  d2/dAxdAz     9 -  d2/dAxdBz     18 -  d2/dAxdCz
c  4 -  d2/dAydAy    10 -  d2/dAydBx     19 -  d2/dAydCx
c  5 -  d2/dAydAz    11 -  d2/dAydBy     20 -  d2/dAydCy
c  6 -  d2/dAzdAz    12 -  d2/dAydBz     21 -  d2/dAydCz
c                    13 -  d2/dAzdBx     22 -  d2/dAzdCx
c                    14 -  d2/dAzdBy     23 -  d2/dAzdCy
c                    15 -  d2/dAzdBz     24 -  d2/dAzdCz
c
c       Block BA :         Block BB :          Block BC :     Block BD :
c                    25 -  d2/dBxdBx     31 -  d2/dBxdCx      trans.inv.
c                    26 -  d2/dBxdBy     32 -  d2/dBxdCy
c                    27 -  d2/dBxdBz     33 -  d2/dBxdCz
c                    28 -  d2/dBydBy     34 -  d2/dBydCx
c                    29 -  d2/dBydBz     35 -  d2/dBydCy
c                    30 -  d2/dBzdBz     36 -  d2/dBydCz
c                                        37 -  d2/dBzdCx
c                                        38 -  d2/dBzdCy
c                                        39 -  d2/dBzdCz
c
c       Block CA :         Block CB :          Block CC :     Block CD :
c                                        40 -  d2/dCxdCx      trans.inv.
c                                        41 -  d2/dCxdCy
c                                        42 -  d2/dCxdCz
c                                        43 -  d2/dCydCy
c                                        44 -  d2/dCydCz
c                                        45 -  d2/dCzdCz
c----------------------------------------------------------------------
c INPUT buf2(1,nbls,lnijr,lnklr,ngcd) - ordinary 2-el.integ.
c
c  buf2(2,nbls,lnijr,lnklr,ngcd) - 2-el.integ. resc by 2*exp_a
c  buf2(3,nbls,lnijr,lnklr,ngcd) - 2-el.integ. resc by 2*exp_b
c  buf2(4,nbls,lnijr,lnklr,ngcd) - 2-el.integ. resc by 2*exp_c
c
c  buf2(5,nbls,lnijr,lnklr,ngcd) - 2-el.integ. resc by 4*exp_a*exp_b
c  buf2(6,nbls,lnijr,lnklr,ngcd) - 2-el.integ. resc by 4*exp_a*exp_c
c  buf2(7,nbls,lnijr,lnklr,ngcd) - 2-el.integ. resc by 4*exp_b*exp_c
c
c  buf2(8,nbls,lnijr,lnklr,ngcd) - 2-el.integ. resc by 4*exp_a*exp_a
c  buf2(9,nbls,lnijr,lnklr,ngcd) - 2-el.integ. resc by 4*exp_b*exp_b
c  buf2(10nbls,lnijr,lnklr,ngcd) - 2-el.integ. resc by 4*exp_c*exp_c
c
c OUTPUT : 
c    der000=der0(nbls,lnij,lnkl,ngcd) - ordinary integr.
c    needed in shifting for d/dA,d/dB (1->2) and d/dC (3->4)
c
c    der1AX=der1(1,nbls,lnij,lnkl,ngcd),
c    der1BX=der1(2,nbls,lnij,lnkl,ngcd)
c    der1CX=der1(3,nbls,lnij,lnkl,ngcd),
c
c    der1AY=der1(4,nbls,lnij,lnkl,ngcd),
c    der1BY=der1(5,nbls,lnij,lnkl,ngcd),
c    der1CY=der1(6,nbls,lnij,lnkl,ngcd),
c
c    der1AZ=der1(7,nbls,lnij,lnkl,ngcd),
c    der1BZ=der1(8,nbls,lnij,lnkl,ngcd),
c    der1CZ=der1(9,nbls,lnij,lnkl,ngcd),
c
c    der2_AxAx=der2(1,nbls,lnij,lnkl,ngcd),
c    der2_AxAy=der2(2,nbls,lnij,lnkl,ngcd),
c    der2_AxAz=der2(3,nbls,lnij,lnkl,ngcd),
c    der2_AyAy=der2(4,nbls,lnij,lnkl,ngcd),
c    der2_AyAz=der2(5,nbls,lnij,lnkl,ngcd),
c    der2_AzAz=der2(6,nbls,lnij,lnkl,ngcd),
c
c    der2_AxBx= der2(7,nbls,lnij,lnkl,ngcd),
c    der2_AxBy= der2(8,nbls,lnij,lnkl,ngcd),
c    der2_AxBz= der2(9,nbls,lnij,lnkl,ngcd),
c    der2_AyBx=der2(10,nbls,lnij,lnkl,ngcd),
c    der2_AyBy=der2(11,nbls,lnij,lnkl,ngcd),
c    der2_AyBz=der2(12,nbls,lnij,lnkl,ngcd),
c    der2_AzBx=der2(13,nbls,lnij,lnkl,ngcd),
c    der2_AzBy=der2(14,nbls,lnij,lnkl,ngcd),
c    der2_AzBz=der2(15,nbls,lnij,lnkl,ngcd),
c
c    der2_AxCx=der2(16,nbls,lnij,lnkl,ngcd),
c    der2_AxCy=der2(17,nbls,lnij,lnkl,ngcd),
c    der2_AxCz=der2(18,nbls,lnij,lnkl,ngcd),
c    der2_AyCx=der2(19,nbls,lnij,lnkl,ngcd),
c    der2_AyCy=der2(20,nbls,lnij,lnkl,ngcd),
c    der2_AyCz=der2(21,nbls,lnij,lnkl,ngcd),
c    der2_AzCx=der2(22,nbls,lnij,lnkl,ngcd),
c    der2_AzCy=der2(23,nbls,lnij,lnkl,ngcd),
c    der2_AzCz=der2(24,nbls,lnij,lnkl,ngcd),
c
c    der2_BxBx=der2(25,nbls,lnij,lnkl,ngcd),
c    der2_BxBy=der2(26,nbls,lnij,lnkl,ngcd),
c    der2_BxBz=der2(27,nbls,lnij,lnkl,ngcd),
c    der2_ByBy=der2(28,nbls,lnij,lnkl,ngcd),
c    der2_ByBz=der2(29,nbls,lnij,lnkl,ngcd),
c    der2_BzBz=der2(30,nbls,lnij,lnkl,ngcd),
c
c    der2_BxCx=der2(31,nbls,lnij,lnkl,ngcd),
c    der2_BxCy=der2(32,nbls,lnij,lnkl,ngcd),
c    der2_BxCz=der2(33,nbls,lnij,lnkl,ngcd),
c    der2_ByCx=der2(34,nbls,lnij,lnkl,ngcd),
c    der2_ByCy=der2(35,nbls,lnij,lnkl,ngcd),
c    der2_ByCz=der2(36,nbls,lnij,lnkl,ngcd),
c    der2_BzCx=der2(37,nbls,lnij,lnkl,ngcd),
c    der2_BzCy=der2(38,nbls,lnij,lnkl,ngcd),
c    der2_BzCz=der2(39,nbls,lnij,lnkl,ngcd),
c
c    der2_CxCx=der2(40,nbls,lnij,lnkl,ngcd),
c    der2_CxCy=der2(41,nbls,lnij,lnkl,ngcd),
c    der2_CxCz=der2(42,nbls,lnij,lnkl,ngcd),
c    der2_CyCy=der2(43,nbls,lnij,lnkl,ngcd),
c    der2_CyCz=der2(44,nbls,lnij,lnkl,ngcd),
c    der2_CzCz=der2(45,nbls,lnij,lnkl,ngcd),
c-------------------------------------------------------------------
c Second derivatives:
c 
c Block AA (1-6)
c Block BB (25-30)
c Block CC (40-45)
c
        nder_aa=1
        nder_bb=25
        nder_cc=40
        do icart=1,3
           do jcart=icart,3    
              call block_aa(icart,jcart,
     *                      ngcd,nbls,buf2,lnijr,lnklr,lnij,lnkl,
     *                      nqij,nqkl,
     *                      nder_aa,der2)
!DIR$ NOINLINE
              call block_bb(icart,jcart,
     *                      ngcd,nbls,buf2,lnijr,lnklr,lnij,lnkl,
     *                      nqij,nqkl,
     *                      nder_bb,der2,xab)
              call block_cc(icart,jcart,
     *                      ngcd,nbls,buf2,lnijr,lnklr,lnij,lnkl,
     *                      nqij,nqkl,
     *                      nder_cc,der2)
              nder_aa=nder_aa+1
              nder_bb=nder_bb+1
              nder_cc=nder_cc+1
           enddo
        enddo
c
c Block AB (7-15)
c Block AC (16-24)
c Block BC (31-39)
c
        nder_ab=7
        nder_ac=16
        nder_bc=31
        do icart=1,3
           do jcart=1,3    
              call block_ab(icart,jcart,
     *                      ngcd,nbls,buf2,lnijr,lnklr,lnij,lnkl,
     *                      nqij,nqkl,
     *                      nder_ab,der2,xab)
              call block_ac(icart,jcart,
     *                      ngcd,nbls,buf2,lnijr,lnklr,lnij,lnkl,
     *                      nqij,nqkl,
     *                      nder_ac,der2)
!DIR$ INLINE
              call block_bc(icart,jcart,
     *                      ngcd,nbls,buf2,lnijr,lnklr,lnij,lnkl,
     *                      nqij,nqkl,
     *                      nder_bc,der2,xab)
              nder_ab=nder_ab+1
              nder_ac=nder_ac+1
              nder_bc=nder_bc+1
           enddo
        enddo
c
      end
c=================================================================
      subroutine block_aa(icart,jcart,
     *                    ngcd,nbls,buf2,lnijr,lnklr,lnij,lnkl,
     *                    nqij,nqkl,
     *                    nder_aa,der2)
      implicit real*8 (a-h,o-z)
#include "texas_lpar.fh"      
c
      dimension buf2(nbls,lnijr,lnklr,ngcd,10)
      dimension der2(45,nbls,lnij,lnkl,ngcd)
c
      do 200 kl=nfu(nqkl)+1,lnkl
      do 200 ij=nfu(nqij)+1,lnij
c
      ij_pi=npxyz(icart,ij)
      ij_pj=npxyz(jcart,ij)
      ij_pi_pj=npxyz(jcart,ij_pi)
c----------------------------------
c Three "pawers" are needed :
c     n_ij_00_i=nia(icart,ij)
c     n_ij_pi_j=nia(jcart,ij_pi)
c     n_ij_mi_j=nia(jcart,ij_mi)
c----------------------------------
      n_ij_pi_j=nia(jcart,ij_pi)
c
      n_ij_00_i=nia(icart,ij)
      if(n_ij_00_i.gt.0) then
         ij_mi=nmxyz(icart,ij)
      else
         ij_mi=0
      endif
c
      if(ij_mi.gt.0) then
         n_ij_mi_j=nia(jcart,ij_mi)
      else
         n_ij_mi_j=0
      endif
c
c first add than substract:
      ij_mi_pj=nmxyz(icart,ij_pj)
      ij_pi_mj=nmxyz(jcart,ij_pi)
c
      ij_mi_mj=0
      if(ij_mi.gt.0) ij_mi_mj=nmxyz(jcart,ij_mi)
c                 or ij_mi_mj=nmxyz(icart,ij_mj)
c
      do 225 iqu=1,ngcd
        do 250 ijkl=1,nbls
c
c       four_a_pipj=buf2(8,ijkl,ij_pi_pj, kl,iqu) 
c       two_a_mi_pj=buf2(2,ijkl,ij_mi_pj, kl,iqu) 
c       two_a_pi_mj=buf2(2,ijkl,ij_pi_mj, kl,iqu) 
c       two_0_mi_mj=buf2(1,ijkl,ij_mi_mj, kl,iqu) 
c
      der=buf2(ijkl,ij_pi_pj, kl,iqu,8) 
c
      if(ij_mi_pj.gt.0) then
         der=der-n_ij_00_i*buf2(ijkl,ij_mi_pj, kl,iqu,2) 
      endif
      if(ij_pi_mj.gt.0) then
         der=der-n_ij_pi_j*buf2(ijkl,ij_pi_mj, kl,iqu,2) 
      endif
      if(ij_mi_mj.gt.0) then
         der=der+n_ij_00_i*n_ij_mi_j*buf2(ijkl,ij_mi_mj, kl,iqu,1)
      endif
c
      der2(nder_aa,ijkl,ij,kl,iqu)=der
c
  250   continue
  225 continue
  200 continue
c
      end
c=================================================================
      subroutine block_cc(icart,jcart,
     *                    ngcd,nbls,buf2,lnijr,lnklr,lnij,lnkl,
     *                    nqij,nqkl,
     *                    nder_cc,der2)
      implicit real*8 (a-h,o-z)
#include "texas_lpar.fh"      
c
c2002 dimension buf2(10,nbls,lnijr,lnklr,ngcd)
      dimension buf2(nbls,lnijr,lnklr,ngcd,10)
      dimension der2(45,nbls,lnij,lnkl,ngcd)
c
      do 200 kl=nfu(nqkl)+1,lnkl
      kl_pi=npxyz(icart,kl)
      kl_pj=npxyz(jcart,kl)
      kl_pi_pj=npxyz(jcart,kl_pi)
c
      n_kl_pi_j=nia(jcart,kl_pi)
c
      n_kl_00_i=nia(icart,kl)
      if(n_kl_00_i.gt.0) then
         kl_mi=nmxyz(icart,kl)
      else
         kl_mi=0
      endif
c
      if(kl_mi.gt.0) then
         n_kl_mi_j=nia(jcart,kl_mi)
      else
         n_kl_mi_j=0
      endif
c
c first add than substract:
      kl_mi_pj=nmxyz(icart,kl_pj)
      kl_pi_mj=nmxyz(jcart,kl_pi)
c
      kl_mi_mj=0
      if(kl_mi.gt.0) kl_mi_mj=nmxyz(jcart,kl_mi)
c                 or kl_mi_mj=nmxyz(icart,kl_mj)
c
      do 200 ij=nfu(nqij)+1,lnij
      do 225 iqu=1,ngcd
        do 250 ijkl=1,nbls
c
c
c       four_c_pipj=buf2(10,ijkl,ij,kl_pi_pj,iqu) 
c       two_c_mi_pj=buf2( 4,ijkl,ij,kl_mi_pj,iqu) 
c       two_c_pi_mj=buf2( 4,ijkl,ij,kl_pi_mj,iqu) 
c       two_0_mi_mj=buf2( 1,ijkl,ij,kl_mi_mj,iqu) 
c
      der=buf2(ijkl,ij,kl_pi_pj,iqu,10) 
      if(kl_mi_pj.gt.0) then
        der=der-n_kl_00_i*buf2(ijkl,ij,kl_mi_pj,iqu,4) 
      endif
      if(kl_pi_mj.gt.0) then
        der=der-n_kl_pi_j*buf2(ijkl,ij,kl_pi_mj,iqu,4) 
      endif
      if(kl_mi_mj.gt.0) then
        der=der+n_kl_00_i*n_kl_mi_j*buf2(ijkl,ij,kl_mi_mj,iqu,1) 
      endif
c
      der2(nder_cc,ijkl,ij,kl,iqu)=der
c
  250   continue
  225 continue
  200 continue
c
      end
c=================================================================
      subroutine block_bb(icart,jcart,
     *                    ngcd,nbls,buf2,lnijr,lnklr,lnij,lnkl,
     *                    nqij,nqkl,
     *                    nder_bb,der2,xab)
      implicit real*8 (a-h,o-z)
#include "texas_lpar.fh"      
c
c2002 dimension buf2(10,nbls,lnijr,lnklr,ngcd)
      dimension buf2(nbls,lnijr,lnklr,ngcd,10)
      dimension der2(45,nbls,lnij,lnkl,ngcd)
      dimension xab(nbls,3)
ctest
c     write(6,60) icart,jcart,nder_bb
c 60  format('from BB: icart,jcart=',2i3,' no=',i3)
c
      do 200 kl=nfu(nqkl)+1,lnkl
      do 200 ij=nfu(nqij)+1,lnij
c
      ij_pi=npxyz(icart,ij)
      ij_pj=npxyz(jcart,ij)
      ij_pi_pj=npxyz(jcart,ij_pi)
c
      do 225 iqu=1,ngcd
        do 250 ijkl=1,nbls
c
c       four_b_pipj=buf2(9,ijkl,ij_pi_pj, kl,iqu) 
c       four_b_pi  =buf2(9,ijkl,ij_pi   , kl,iqu) 
c       four_b_pj  =buf2(9,ijkl,ij_pj   , kl,iqu) 
c       four_b_0   =buf2(9,ijkl,ij      , kl,iqu) 
c       two_b_0    =buf2(3,ijkl,ij      , kl,iqu) 
c
      der=                  buf2(ijkl,ij_pi_pj,kl,iqu,9) 
     *              +xab(ijkl,jcart)*buf2(ijkl,ij_pi,kl,iqu,9) 
     *              +xab(ijkl,icart)*buf2(ijkl,ij_pj,kl,iqu,9) 
     *   +xab(ijkl,icart)*xab(ijkl,jcart)*buf2(ijkl,ij,kl,iqu,9) 

      if(jcart.eq.icart) der=der-buf2(ijkl,ij,kl,iqu,3) 
c
      der2(nder_bb,ijkl,ij,kl,iqu)=der
c
c        write(6,66) nder_bb, ij,kl
c  66    format('no=',i3,'ij,kl=',2i4)
c        write(6,67) '(2b)2 int. ++ =',buf2(9,ijkl,ij_pi_pj,kl,iqu)
c        write(6,67) '(2b)2 int. +0 =',buf2(9,ijkl,ij_pi   ,kl,iqu)
c        write(6,67) '(2b)2 int. 0+ =',buf2(9,ijkl,ij_pj   ,kl,iqu)
c        write(6,67) '(2b)2 int. 00 =',buf2(9,ijkl,ij      ,kl,iqu)
c        write(6,67) '(2b)1 int. 00 =',buf2(3,ijkl,ij      ,kl,iqu)
c        write(6,67) '  final deriv =',der
c 67     format(a15,f12.7)
c
  250   continue
  225 continue
  200 continue
c
      end
c=================================================================
      subroutine block_ab(icart,jcart,
     *                    ngcd,nbls,buf2,lnijr,lnklr,lnij,lnkl,
     *                    nqij,nqkl,
     *                    nder_ab,der2,xab)
      implicit real*8 (a-h,o-z)
#include "texas_lpar.fh"      
c
c2002 dimension buf2(10,nbls,lnijr,lnklr,ngcd)
      dimension buf2(nbls,lnijr,lnklr,ngcd,10)
      dimension der2(45,nbls,lnij,lnkl,ngcd)
      dimension xab(nbls,3)
ctest
c     write(6,60) icart,jcart,nder_ab
c 60  format('from AB: icart,jcart=',2i3,' no=',i3)
c
c
      do 200 kl=nfu(nqkl)+1,lnkl
      do 200 ij=nfu(nqij)+1,lnij
c
      ij_pi=npxyz(icart,ij)
      ij_pj=npxyz(jcart,ij)
      n_ij_i=nia(icart,ij)
c
      ij_mi=0
      if(n_ij_i.gt.0) ij_mi=nmxyz(icart,ij)
c
      ij_pi_pj=npxyz(jcart,ij_pi)
c first add than substract:
      ij_mi_pj=nmxyz(icart,ij_pj)
c
c     ij_mi_mj=0
c     if(ij_mi.gt.0) ij_mi_mj=nmxyz(jcart,ij_mi)
c                 or ij_mi_mj=nmxyz(icart,ij_mj)
c
      do 225 iqu=1,ngcd
        do 250 ijkl=1,nbls
c
c       four_ab_pipj=buf2(5,ijkl,ij_pi_pj, kl,iqu) 
c       four_ab_pi  =buf2(5,ijkl,ij_pi   , kl,iqu) 
c
c       two_b_mi_pj=buf2(3,ijkl,ij_mi_pj, kl,iqu) 
c       two_b_mi   =buf2(3,ijkl,ij_mi   , kl,iqu) 
c
      der=                buf2(ijkl,ij_pi_pj, kl,iqu,5) 
     *         +xab(ijkl,jcart)*buf2(ijkl,ij_pi, kl,iqu,5) 
c
      if(n_ij_i.gt.0) then
        if(ij_mi_pj.gt.0) then
          der=der-n_ij_i*buf2(ijkl,ij_mi_pj, kl,iqu,3) 
        endif
        if(ij_mi   .gt.0) then
          der=der -n_ij_i*xab(ijkl,jcart)*buf2(ijkl,ij_mi , kl,iqu,3) 
        endif
      endif
c
c        write(6,66) nder_ab, ij,kl
c  66    format('no=',i3,'ij,kl=',2i4)
c        write(6,67) '(4ab  int. ++ =',buf2(5,ijkl,ij_pi_pj,kl,iqu)
c        write(6,67) '(4ab  int. +  =',buf2(5,ijkl,ij_pi   ,kl,iqu)
c        write(6,67) '(2b   int. -+ =',buf2(3,ijkl,ij_mi_pj,kl,iqu)
c        write(6,67) '(2b   int. -  =',buf2(3,ijkl,ij_mi   ,kl,iqu)
c        write(6,67) '  final deriv =',der
c 67     format(a15,f12.7)
c
c
      der2(nder_ab,ijkl,ij,kl,iqu)=der
c
  250   continue
  225 continue
  200 continue
c
      end
c=================================================================
      subroutine block_ac(icart,jcart,
     *                    ngcd,nbls,buf2,lnijr,lnklr,lnij,lnkl,
     *                    nqij,nqkl,
     *                    nder_ac,der2)
      implicit real*8 (a-h,o-z)
#include "texas_lpar.fh"      
c
c2002 dimension buf2(10,nbls,lnijr,lnklr,ngcd)
      dimension buf2(nbls,lnijr,lnklr,ngcd,10)
      dimension der2(45,nbls,lnij,lnkl,ngcd)
ctest
c     write(6,60) icart,jcart,nder_ac
c 60  format('from AC: icart,jcart=',2i3,' no=',i3)
c
c
      do 200 kl=nfu(nqkl)+1,lnkl
      kl_pj=npxyz(jcart,kl)
      n_kl_j=nia(jcart,kl)
      kl_mj=0
      if(n_kl_j.gt.0) kl_mj=nmxyz(jcart,kl)
c
      do 200 ij=nfu(nqij)+1,lnij
      ij_pi=npxyz(icart,ij)
      n_ij_i=nia(icart,ij)
      ij_mi=0
      if(n_ij_i.gt.0) ij_mi=nmxyz(icart,ij)
c
      do 225 iqu=1,ngcd
        do 250 ijkl=1,nbls
c
c       four_ac_pipj=buf2(6,ijkl,ij_pi,kl_pj,iqu) 
c
c       two_c_mi_pj=buf2(4,ijkl,ij_mi,kl_pj,iqu) 
c       two_a_pi_mj=buf2(2,ijkl,ij_pi,kl_mj,iqu) 
c       two_0_mi_mj=buf2(1,ijkl,ij_mi,kl_mj,iqu) 
c
      der=buf2(ijkl,ij_pi,kl_pj,iqu,6) 
      if(n_ij_i.gt.0 .and. ij_mi.gt.0) then
        der=der-n_ij_i*buf2(ijkl,ij_mi,kl_pj,iqu,4) 
      endif
      if(n_kl_j.gt.0 .and. kl_mj.gt.0) then
        der=der-n_kl_j*buf2(ijkl,ij_pi,kl_mj,iqu,2) 
      endif
      if(n_ij_i.gt.0 .and. n_kl_j.gt.0) then
cccc    if(ij_mi.gt.0 .and. kl_mj.gt.0) then
        der=der+n_ij_i*n_kl_j*buf2(ijkl,ij_mi,kl_mj,iqu,1) 
cccc    endif
      endif
c
      der2(nder_ac,ijkl,ij,kl,iqu)=der
c
c        write(6,66) nder_ac, ij,kl
c        write(6,*)'ij   ,kl   =',ij   ,kl   
c        write(6,*)'ij_pi,kl_pj=',ij_pi,kl_pj
c        write(6,*)'ij_mi,kl_pj=',ij_mi,kl_pj
c        write(6,*)'ij_pi,kl_mj=',ij_pi,kl_mj
c        write(6,*)'ij_mi,kl_mj=',ij_mi,kl_mj
c  66    format('no=',i3,'ij,kl=',2i4)
c        write(6,67) '(4ac  int. ++ =',buf2(6,ijkl,ij_pi,kl_pj,iqu)
c        write(6,67) '(2c   int. -+ =',buf2(4,ijkl,ij_mi,kl_pj,iqu)
c        write(6,67) '(2a   int. +- =',buf2(2,ijkl,ij_pi,kl_mj,iqu)
c        write(6,67) '(     int. -- =',buf2(1,ijkl,ij_mi,kl_mj,iqu)
c        write(6,67) '  final deriv =',der
c 67     format(a15,f12.7)
c
  250   continue
  225 continue
  200 continue
c
      end
c=================================================================
      subroutine block_bc(icart,jcart,
     *                    ngcd,nbls,buf2,lnijr,lnklr,lnij,lnkl,
     *                    nqij,nqkl,
     *                    nder_bc,der2,xab)
      implicit real*8 (a-h,o-z)
#include "texas_lpar.fh"      
c
c2002 dimension buf2(10,nbls,lnijr,lnklr,ngcd)
      dimension buf2(nbls,lnijr,lnklr,ngcd,10)
      dimension der2(45,nbls,lnij,lnkl,ngcd)
      dimension xab(nbls,3)
c
      do 200 kl=nfu(nqkl)+1,lnkl
      kl_pj=npxyz(jcart,kl)
      n_kl_j=nia(jcart,kl)
      kl_mj=0
      if(n_kl_j.gt.0) kl_mj=nmxyz(jcart,kl)
c
      do 200 ij=nfu(nqij)+1,lnij
      ij_pi=npxyz(icart,ij)
c
      do 225 iqu=1,ngcd
        do 250 ijkl=1,nbls
c
c       four_bc_pipj=buf2(7,ijkl,ij_pi,kl_pj,iqu) 
c       four_bc_0ipj=buf2(7,ijkl,ij   ,kl_pj,iqu) 
c
c       two_b_pi_mj=buf2(3,ijkl,ij_pi,kl_mj,iqu) 
c       two_b_0i_mj=buf2(3,ijkl,ij   ,kl_mj,iqu) 
c
c
      der=                buf2(ijkl,ij_pi,kl_pj,iqu,7) 
     *  + xab(ijkl,icart)*buf2(ijkl,ij   ,kl_pj,iqu,7)
c
      if(n_kl_j.gt.0 .and. kl_mj.gt.0) then
          der=der
     *    - n_kl_j*(                  buf2(ijkl,ij_pi,kl_mj,iqu,3) 
     *              + xab(ijkl,icart)*buf2(ijkl,ij   ,kl_mj,iqu,3) )
     *      
      endif
c
      der2(nder_bc,ijkl,ij,kl,iqu)=der
c
  250   continue
  225 continue
  200 continue
c
      end
c=================================================================
