* $Id$
c-----------------------------------------------------------------------
C This files contains a set of routines used to convert pair
c quantities into quartet quantities. These routines used to
c be included in different files, usually where they are called
c from. They are called with ...BL(IADDRESS).. double precision 
c parameters which then are of the INTEGER type in the routines. 
c It causes problems of incompatibile types for some compilers.
c-----------------------------------------------------------------------
c List of subroutines with previous location :
c
c amshift.f:      subroutine convr3(bl,m,nbls,npij,npkl,idx1,idx2,
c 
c derivat.f:      subroutine conv24x(nbls,npij,npkl,idx1,idx2 ,
c derivat.f:      subroutine conv24r(nbls,npij,idx1,xab,xabq)
c 
c spec_calcint.f: subroutine conv1x_1(nbls1,mmax1,npij,lcij, idx1,indx,
c spec_calcint.f: subroutine conv1der_1(nbls1,npij,lci,idx1,indx, aa, 
c spec_calcint.f: subroutine conv1der_2(nbls1,npij,lci,idx1,indx, aa,
c spec_calcint.f: subroutine conv1x_2(nbls1,mmax1,npij,lcij, idx1,indx,
c spec_calcint.f: subroutine conv2x(nbls1,nfumax1,npkl,lckl, idx2,indx,
c-----------------------------------------------------------------------
c=======================================================================
      subroutine convr3(bl,m,nbls,npij,npkl,idx1,idx2,
     *                   xab,xcd, ixabn,ixcdn)
      implicit real*8 (a-h,o-z)
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
c=======================================================================
      subroutine conv24x(nbls,npij,npkl,idx1,idx2 ,
     *                  xab ,xcd, xyab, xycd ,
     *                  xabq,xcdq,xyabq,xycdq )
      implicit real*8 (a-h,o-z)
c
      dimension idx1(nbls),idx2(nbls)
      dimension xab(npij,3) ,xcd(npkl,3) ,xyab(npij,3) ,xycd(npkl,3)
      dimension xabq(nbls,3),xcdq(nbls,3),xyabq(nbls,3),xycdq(nbls,3)
c
      do 100 ijkl=1,nbls
      ijpar=idx1(ijkl)
      klpar=idx2(ijkl)
        do 150 i=1,3
        xabq(ijkl,i)=xab(ijpar,i)
        xcdq(ijkl,i)=xcd(klpar,i)
        xyabq(ijkl,i)=xyab(ijpar,i)
        xycdq(ijkl,i)=xycd(klpar,i)
  150   continue
  100 continue
c
      end
c=======================================================================
      subroutine conv24r(nbls,npij,idx1,xab,xabq)
      implicit real*8 (a-h,o-z)
c
      dimension idx1(nbls)
      dimension xab(npij,3),xabq(nbls,3)
c
      do 100 ijkl=1,nbls
      ijpar=idx1(ijkl)
c     klpar=idx2(ijkl)
        do 150 i=1,3
        xabq(ijkl,i)=xab(ijpar,i)
c       xcdq(ijkl,i)=xcd(klpar,i)
  150   continue
  100 continue
      end
c=======================================================================
      subroutine conv1x_1(nbls1,mmax1,npij,lcij, idx1,indx,
     *                    abnia,xpn,abnix,xpnx )
c-------------------------------------------------------------------
c npij = number of uniqe pairs now
c-------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      dimension idx1(*),indx(*)
      dimension xpn(npij,3,*)
      dimension abnia(npij,mmax1,*)
c
      dimension xpnx(nbls1,3)
      dimension abnix(nbls1,mmax1)
c 
      do 10 i=1,nbls1
      ijkl=indx(i)
      ijpar=idx1(ijkl)
        xpnx(i,1)=xpn(ijpar,1,lcij)
        xpnx(i,2)=xpn(ijpar,2,lcij)
        xpnx(i,3)=xpn(ijpar,3,lcij)
   10 continue
c
      do 20 m=1,mmax1
      do 20 i=1,nbls1
      ijkl=indx(i)
      ijpar=idx1(ijkl)
        abnix(i,m)=abnia(ijpar,m,lcij)
   20 continue
      end
c=======================================================================
      subroutine conv1x_2(nbls1,mmax1,npij,lcij, idx1,indx,xpn,xpnx )
c
c npij = number of uniqe pairs now
c
c
      implicit real*8 (a-h,o-z)
      dimension idx1(*),indx(*)
      dimension xpn(npij,3,*)
      dimension xpnx(nbls1,3)
c 
      do 10 i=1,nbls1
      ijkl=indx(i)
      ijpar=idx1(ijkl)
        xpnx(i,1)=xpn(ijpar,1,lcij)
        xpnx(i,2)=xpn(ijpar,2,lcij)
        xpnx(i,3)=xpn(ijpar,3,lcij)
   10 continue
c
      end
c=======================================================================
      subroutine conv1der_1(nbls1,npij,lci,idx1,indx, aa, aax)
c
c npij = number of uniqe pairs now
c exponents are already rescaled by 2 in precal2a_1
c
      implicit real*8 (a-h,o-z)
      dimension idx1(*),indx(*)
      dimension aa(npij,*)
c output :
      dimension aax(nbls1)
c 
      do 10 i=1,nbls1
      ijkl=indx(i)
      ijpar=idx1(ijkl)
cccc    aax(i)=aa(ijpar,lci)*2.0d0   ! already rescaled
        aax(i)=aa(ijpar,lci)
   10 continue
c
      end
c=======================================================================
      subroutine conv1der_2(nbls1,npij,lci,idx1,indx, aa, aax)
c
c npij = number of uniqe pairs now
c
      implicit real*8 (a-h,o-z)
      dimension idx1(*),indx(*)
      dimension aa(npij,*)
c output :
      dimension aax(nbls1)
c-------------------------------------------------------------------
c this is for iroute=2 only: all exponents in a block are the same:
c exponents are already rescaled by 2 in precal2a_2
c-------------------------------------------------------------------
      aax(1)=aa(1,lci)
c
c     ijkl1 =indx(1)
c     ijpar1=idx1(ijkl1)
c     aax(1)=aa(ijpar1,lci)*2.d0
c----------------------------------------------
c     do 10 i=1,nbls1
c     ijkl=indx(i)
c     ijpar=idx1(ijkl)
c       aax(i)=aa(ijpar,lci)*2.0d0
c  10 continue
c----------------------------------------------
c
      end
c=======================================================================
      subroutine conv2x(nbls1,nfumax1,npkl,lckl, idx2,indx,
     *                  habcd,nfumax, habcdx  )
c
c npkl = number of uniqe pairs now
c
      implicit real*8 (a-h,o-z)
      dimension idx2(*),indx(*)
      dimension habcd(npkl,3,nfumax,*)
      dimension habcdx(nbls1,3,nfumax)
c
         do 32 ifu=1,nfumax1
         do 32 i=1,nbls1
         ijkl=indx(i)
         klpar=idx2(ijkl)
           habcdx(i,1,ifu)=habcd(klpar,1,ifu,lckl)
           habcdx(i,2,ifu)=habcd(klpar,2,ifu,lckl)
           habcdx(i,3,ifu)=habcd(klpar,3,ifu,lckl)
   32    continue
      end
c=======================================================================
