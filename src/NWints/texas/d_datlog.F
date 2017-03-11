c $Id$
C_______________________________________________________
c OCT.,94 KW changed the order of D-orbitals to the consistent
C one : xx xy xz yy yz zz
C
C Corresponding changes have been made in one-electron part
C (intcal1.f in this directory).
C
C The NMR part has NOT been changed .
C_______________________________________________________
      subroutine datlog(inx,ncs,lp1,lp2,lp3,lp4,lp5)
c    *             hnia,ndege,len,lensm,nfu,icoor,icool,
c    *             ifrst,ilast,nia,nmxyz,npxyz)
      implicit real*8 (a-h,o-z)
#include "errquit.fh"
c***********************************************
c* memory for all logical matrices - in commons logicd, logic1-11
c*
c*  denote  : maxtyp1=max(itype1)
c*            itype1=inx(12,ics)
c*            if(itype1>4) itype1=itype1-1
c*            if(itype1>5) itype1=itype1-1
c*
c*  denote  : mmax=4*ndege(maxtype1)-3
c*
c*  Array                  Dimension           Parameter
c*----------------------------------------------------------
c*  ndege, len          -  (maxtype1)            lp4
c*                                      
c*  nfu                 -   (mmax+1)             lp1
c*  lensm               -   ( mmax )             lp5
c
c*  icoor,icool,
c*  ifrst,ilast         -   ( nfu(mmax+1) )      lp2
c*
c*  nia                 -   ( 3, nfu(mmax+1) )   3 * lp2
c*  hnia                -   ( 3, nfu(mmax+1) )   3 * lp2
c*  nmxyz               -   ( 3, nfu(mmax+1) )   3 * lp2
c*  npxyz               -   ( 3, nfu(mmax  ) )   3 * lp3
c
c*
c*  shells maxtyp1 mmax nfu(mmax) nfu(mmax+1) lp1,2, 3,  4,  5
c* -------------------------------------------------------------
c*    ss,ss    1     1      0       1       2,    1,    0,  1,  1
c*    pp,pp    2     5     20      35       6,   35,   20,  2,  5
c*    ll,ll    3     5     20      35       6,   35,   20,  3,  5
c*    dd,dd    4     9    120     165      10,  165,  120,  4,  9 
c*    ff,ff    5    13    364     455      14,  455,  364,  5, 13
c*    gg,gg    6    17    816     969      18,  969,  816,  6, 17 
c*    hh,hh    7    21   1540     1771     22, 1771, 1540,  7, 21
c*    ii,ii    8    25   2600     2925     26, 2925, 2600,  8, 25
c*
c* 
c* --------------------------------------------------------------
c dimensions for logical matrices in TWELINT :
c
c up to ff,ff :
c     parameter (lpar1=14,lpar2= 455,lpar3= 364,lpar4=5,lpar5=13)
c up to gg,gg :
c     parameter (lpar1=18,lpar2= 969,lpar3= 816,lpar4=6,lpar5=17)
c up to hh,hh :
c     parameter (lpar1=22,lpar2=1771,lpar3=1540,lpar4=7,lpar5=21)
c up to ii,ii :
c     parameter (lpar1=26,lpar2=2925,lpar3=2600,lpar4=8,lpar5=25)
c
c---------------------------------------------------------------
      dimension inx(12,*)
cxx
c     dimension  hnia(3,lp2)
c     dimension  ndege(lp4)
c     dimension  len(lp4)
c     dimension  lensm(lp5)
c     dimension  nfu(lp1)
c     dimension  icoor(lp2)
c     dimension  icool(lp2)
c     dimension  ifrst(lp2)
c     dimension  ilast(lp2)
c     dimension  nia(3,lp2)
c     dimension  nmxyz(3,lp2)
c     dimension  npxyz(3,lp3)
cxx
Cedo      parameter (lpar1=30,lpar2=4495,lpar3=4060,lpar4=9,lpar5=29)
#include "texas_lpar.fh"
cxx
c
c     lensm(nsij)=total number of functions up to given nsij
c************************************************************
c
      maxtyp1=0
      do 10 ics=1,ncs
      itype=inx(12,ics)
      itype1=itype
      if(itype.gt.4) itype1=itype-1
      if(itype1.gt.5) itype1=itype1-1
c
      if(itype1.gt.maxtyp1) maxtyp1=itype1
   10 continue
c
c for derivatives of two-el. integ. :
       maxtyp1=maxtyp1+2
c
c
c* check if dimensions of logical arrays are correct :
c
      if(maxtyp1.gt.lp4) then
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        call txs_message
     &      ('datlog :','logic arrays too small,lpar4=',lp4,0)
        call txs_message
     &      ('datlog :','take set of parameters with lpar4=',
     *                maxtyp1,0)
        call txs_message
     &      ('datlog :',' make this change in PREPINT2 ',0,0)
        call txs_message
     &      ('datlog :',' execution stoped ', 0,0)
        call errquit('texas:datlog',0, INT_ERR)
*        stop10
      endif
c
      ndege(1)=1
      ndege(2)=2
      ndege(3)=2
      len(1)=1
      len(2)=3
      len(3)=4
      do 20 ityp1=4,maxtyp1
      ityp11=ityp1-1
      ndege(ityp1)=ityp11
      len(ityp1)=ityp1*ityp11/2
   20 continue
c
      nqi=ndege(maxtyp1)
      mmax=4*nqi-3
      lensm(1)=1
      lensm(2)=4
      nfu(1)=0
      nfu(2)=1
      nfu(3)=4
      do 30 i=3,mmax 
      lensm(i)=lensm(i-1) + i*(i+1)/2
      nfu(i+1)=lensm(i)
   30 continue
c
c nia :
c
      do 40 i=1,10
      nia(1,i)=0
      nia(2,i)=0
      nia(3,i)=0
   40 continue
      nia(1,2)=1
      nia(2,3)=1
      nia(3,4)=1
c d-
      nia(1,5)=2
      nia(2,6)=2
      nia(3,7)=2
      nia(1,8)=1
      nia(2,8)=1
      nia(1,9)=1
      nia(3,9)=1
      nia(2,10)=1
      nia(3,10)=1
c d-new
      nia(1,5)=2
      nia(2,5)=0
      nia(3,5)=0
c
      nia(1,6)=1
      nia(2,6)=1
      nia(3,6)=0
c
      nia(1,7)=1
      nia(2,7)=0
      nia(3,7)=1
c
      nia(1,8)=0
      nia(2,8)=2
      nia(3,8)=0
c
      nia(1,9)=0 
      nia(2,9)=1
      nia(3,9)=1
c
      nia(1,10)=0
      nia(2,10)=0
      nia(3,10)=2
c f-
      ijk=10
      do 43 i=1,3
      do 43 j=i,3
      do 43 k=j,3
      ijk=ijk+1
      ix=0
      iy=0
      iz=0
        if(i.eq.1) ix=ix+1
        if(j.eq.1) ix=ix+1
        if(k.eq.1) ix=ix+1
c
        if(i.eq.2) iy=iy+1
        if(j.eq.2) iy=iy+1
        if(k.eq.2) iy=iy+1
c
        if(i.eq.3) iz=iz+1
        if(j.eq.3) iz=iz+1
        if(k.eq.3) iz=iz+1
c
        nia(1,ijk)=ix
        nia(2,ijk)=iy
        nia(3,ijk)=iz
c
   43 continue
c
c g- and higher
c
      do 50 nq=5,mmax
      nq1=nq-1
      nful=nfu(nq)-nfu(nq1)
cxxx  nfuc=nfu(nq+1)-nfu(nq)
      iful=nfu(nq1)
      ifuc=nfu(nq)
         do 51 i=1,nful
         iful=iful+1
         ifuc=ifuc+1
         nia(1,ifuc)=nia(1,iful)+1
         nia(2,ifuc)=nia(2,iful)
         nia(3,ifuc)=nia(3,iful)
   51    continue
c
         do 52 i=1,nq
         i1=i-1
         ifuc=ifuc+1
         nia(1,ifuc)=0
         nia(2,ifuc)=nq1-i1
         nia(3,ifuc)=i1
   52    continue
   50 continue
c
c* total number of functions
c
      nfun=nfu(mmax+1)
      nfu1=nfu(mmax)
c
c* hnia matrix :
c
      do 55 i=1,nfun
      hnia(1,i)=0.5d0*dble(nia(1,i))
      hnia(2,i)=0.5d0*dble(nia(2,i))
      hnia(3,i)=0.5d0*dble(nia(3,i))
   55 continue
c
c* nmxyz and npxyz matrices :
c
      do 60 nf=1,nfun
      ix=nia(1,nf)
      iy=nia(2,nf)
      iz=nia(3,nf)
          do 65 nf1=1,nfun
          ix1=nia(1,nf1)
          iy1=nia(2,nf1)
          iz1=nia(3,nf1)
      if(ix1.eq.ix-1.and.iy1.eq.iy.and.iz1.eq.iz) nmxyz(1,nf)=nf1
      if(ix1.eq.ix.and.iy1.eq.iy-1.and.iz1.eq.iz) nmxyz(2,nf)=nf1
      if(ix1.eq.ix.and.iy1.eq.iy.and.iz1.eq.iz-1) nmxyz(3,nf)=nf1
ccc
      if(nf.le.nfu1) then
      if(ix1.eq.ix+1.and.iy1.eq.iy.and.iz1.eq.iz) npxyz(1,nf)=nf1
      if(ix1.eq.ix.and.iy1.eq.iy+1.and.iz1.eq.iz) npxyz(2,nf)=nf1
      if(ix1.eq.ix.and.iy1.eq.iy.and.iz1.eq.iz+1) npxyz(3,nf)=nf1
      endif
   65     continue
   60 continue
c
c* icoor, icool and ifrst, ilast marices :
c
ctest icool(1)=0
ctest icoor(1)=0
      icool(1)=1
      icoor(1)=1
      ifrst(1)=1
      ilast(1)=1
      do 70 nf=2,nfun
      ix=nia(1,nf)
      iy=nia(2,nf)
      iz=nia(3,nf)
         if(ix.ne.0) then
           icool(nf)=1
         else if(iy.ne.0) then
           icool(nf)=2
         else
           icool(nf)=3
         endif
cc
         if(iz.ne.0) then
           icoor(nf)=3
         else if(iy.ne.0) then
           icoor(nf)=2
         else
           icoor(nf)=1
         endif
c
      ilast(nf)=nmxyz(icool(nf),nf)
      ifrst(nf)=nmxyz(icoor(nf),nf)
c
  70  continue
c
      return
      end
ccccccc
