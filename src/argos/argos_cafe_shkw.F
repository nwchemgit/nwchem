      subroutine argos_cafe_shkw(nwmloc,iwgm,iwdt,xw,yw,mwb,
     + nwb,nbp,ibnd,bnd,iwat,wgt,rwx,rwi1,rwi2,rwi6,f,swt)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      integer nwmloc,iwgm(mwm),iwdt(mwm)
      real*8 xw(mwm,3,mwa),yw(mwm,3,mwa)
      integer mwb,nwb,nbp
      integer ibnd(mwb,3),iwat(mwa)
      real*8 bnd(mwb,nbp,6),wgt(mat,mset)
      real*8 rwx(mscr,3),rwi1(mscr),rwi2(mscr),rwi6(mscr),f(mscr)
      real*8 swt(mwm,3,mwb)
c
      logical ready
c
      integer ix,iwb,iwm,iwa,jwa,nshake,ndtl
      real*8 bw,tw,w,dif,tisum,cdwb4,ydif1,ydif2,ydif3
      real*8 ydist,ep2sum,cdwb2,ep3sum,cdwb3,wiwa,wjwa
c
#include "bitops.fh"
c
      if(ithint.and.ith(7)) then
      do 6 ix=1,3
      do 5 iwb=1,nwb
      do 4 iwm=1,nwmloc
      swt(iwm,ix,iwb)=zero
    4 continue
    5 continue
    6 continue
      endif
c
      niterw=0
    1 continue
      ready=.true.
      do 2 iwb=1,nwb
      if(iand(ibnd(iwb,3),icnstr).eq.0) goto 2
      iwa=ibnd(iwb,1)
      jwa=ibnd(iwb,2)
      wiwa=one/wgt(iwat(iwa),iset)
      wjwa=one/wgt(iwat(jwa),iset)
      bw=bnd(iwb,1,iset)*bnd(iwb,1,iset)
      tw=bw*tolsha
      w=wiwa+wjwa
      nshake=0
      do 3 iwm=1,nwmloc
      rwx(iwm,1)=xw(iwm,1,iwa)-xw(iwm,1,jwa)
      rwx(iwm,2)=xw(iwm,2,iwa)-xw(iwm,2,jwa)
      rwx(iwm,3)=xw(iwm,3,iwa)-xw(iwm,3,jwa)
      rwi1(iwm)=bw-rwx(iwm,1)**2-rwx(iwm,2)**2-rwx(iwm,3)**2
      if(abs(rwi1(iwm)).gt.tw) nshake=nshake+1
    3 continue
      if(nshake.gt.0) then
      ready=.false.
      do 7 iwm=1,nwmloc
      dif=yw(iwm,1,iwa)-yw(iwm,1,jwa)
      rwi2(iwm)=dif*rwx(iwm,1)
      rwx(iwm,1)=dif
      dif=yw(iwm,2,iwa)-yw(iwm,2,jwa)
      rwi2(iwm)=rwi2(iwm)+dif*rwx(iwm,2)
      rwx(iwm,2)=dif
      dif=yw(iwm,3,iwa)-yw(iwm,3,jwa)
      rwi2(iwm)=rwi2(iwm)+dif*rwx(iwm,3)
      rwx(iwm,3)=dif
    7 continue
      ndtl=0
      do 8 iwm=1,nwmloc
      if(rwi2(iwm).lt.small) ndtl=iwgm(iwm)
    8 continue
      if(ndtl.gt.0)
     + call md_abort('Deviation too large for solvent',ndtl)
      do 9 iwm=1,nwmloc
      rwi6(iwm)=half*rwi1(iwm)/(rwi2(iwm)*w)
    9 continue
c
      if(ithint.and.ith(7)) then
      do 12 ix=1,3
      do 11 iwm=1,nwmloc
      swt(iwm,ix,iwb)=swt(iwm,ix,iwb)+rwi6(iwm)*rwx(iwm,ix)
   11 continue
   12 continue
      endif
c
      do 10 iwm=1,nwmloc
      f(iwm)=rwi6(iwm)*rwx(iwm,1)
      xw(iwm,1,iwa)=xw(iwm,1,iwa)+f(iwm)*wiwa
      xw(iwm,1,jwa)=xw(iwm,1,jwa)-f(iwm)*wjwa
      f(iwm)=rwi6(iwm)*rwx(iwm,2)
      xw(iwm,2,iwa)=xw(iwm,2,iwa)+f(iwm)*wiwa
      xw(iwm,2,jwa)=xw(iwm,2,jwa)-f(iwm)*wjwa
      f(iwm)=rwi6(iwm)*rwx(iwm,3)
      xw(iwm,3,iwa)=xw(iwm,3,iwa)+f(iwm)*wiwa
      xw(iwm,3,jwa)=xw(iwm,3,jwa)-f(iwm)*wjwa
   10 continue
      endif
    2 continue
      niterw=niterw+1
      if(niterw.gt.mshitw)
     + call md_abort('Too many shake iterations',niterw)
      if(.not.ready) goto 1
c
      do 24 iwm=1,nwmloc
      if(iand(iwdt(iwm),mfixed).eq.lfixed) then
      do 25 iwa=1,mwa
      xw(iwm,1,iwa)=yw(iwm,1,iwa)
      xw(iwm,2,iwa)=yw(iwm,2,iwa)
      xw(iwm,3,iwa)=yw(iwm,3,iwa)
   25 continue
      if(ithint.and.ith(7)) then
      do 26 iwb=1,nwb
      swt(iwm,1,iwb)=zero
      swt(iwm,2,iwb)=zero
      swt(iwm,3,iwb)=zero
   26 continue
      endif
      endif
   24 continue
c
c      if(nwr.gt.0) then
c      if(iropt.eq.2.or.(iropt.eq.1.and.forest.gt.tiny)) then
c      do 21 iwr=1,nwr
c      cmx1=zero
c      cmx2=zero
c      cmx3=zero
c      do 22 iwa=1,mwa
c      cmx1=cmx1+xw(idwr(iwr),1,iwa)*wwaf(iwa)
c      cmx2=cmx2+xw(idwr(iwr),2,iwa)*wwaf(iwa)
c      cmx3=cmx3+xw(idwr(iwr),3,iwa)*wwaf(iwa)
c   22 continue
c      do 23 iwa=1,mwa
c      xw(idwr(iwr),1,iwa)=xw(idwr(iwr),1,iwa)+xwr(iwr,1)-cmx1
c      xw(idwr(iwr),2,iwa)=xw(idwr(iwr),2,iwa)+xwr(iwr,2)-cmx2
c      xw(idwr(iwr),3,iwa)=xw(idwr(iwr),3,iwa)+xwr(iwr,3)-cmx3
c   23 continue
c   21 continue
c      endif
c      endif
c
      if(ithint.and.ith(7)) then
      tisum=zero
      do 14 iwb=1,nwb
      iwa=ibnd(iwb,1)
      jwa=ibnd(iwb,2)
      cdwb4=bnd(iwb,1,4)
      do 15 iwm=1,nwmloc
      ydif1=yw(iwm,1,iwa)-yw(iwm,1,jwa)
      ydif2=yw(iwm,2,iwa)-yw(iwm,2,jwa)
      ydif3=yw(iwm,3,iwa)-yw(iwm,3,jwa)
      ydist=sqrt(ydif1**2+ydif2**2+ydif3**2)
      tisum=tisum+(swt(iwm,1,iwb)*ydif1+swt(iwm,2,iwb)*ydif2+
     + swt(iwm,3,iwb)*ydif3)*cdwb4/ydist
   15 continue
   14 continue
      deriv(7,1)=tisum*tstepi*tstepi
      endif
      if(ip2(5)) then
      ep2sum=zero
      do 16 iwb=1,nwb
      iwa=ibnd(iwb,1)
      jwa=ibnd(iwb,2)
      cdwb2=bnd(iwb,1,2)-bnd(iwb,1,1)
      do 17 iwm=1,nwmloc
      ydif1=yw(iwm,1,iwa)-yw(iwm,1,jwa)
      ydif2=yw(iwm,2,iwa)-yw(iwm,2,jwa)
      ydif3=yw(iwm,3,iwa)-yw(iwm,3,jwa)
      ydist=sqrt(ydif1**2+ydif2**2+ydif3**2)
      ep2sum=ep2sum+(swt(iwm,1,iwb)*ydif1+swt(iwm,2,iwb)*ydif2+
     + swt(iwm,3,iwb)*ydif3)*cdwb2/ydist
   17 continue
   16 continue
      ep2(1)=ep2(1)+ep2sum*tstepi*tstepi
      endif
      if(ip3(5)) then
      ep3sum=zero
      do 18 iwb=1,nwb
      iwa=ibnd(iwb,1)
      jwa=ibnd(iwb,2)
      cdwb3=bnd(iwb,1,3)-bnd(iwb,1,1)
      do 19 iwm=1,nwmloc
      ydif1=yw(iwm,1,iwa)-yw(iwm,1,jwa)
      ydif2=yw(iwm,2,iwa)-yw(iwm,2,jwa)
      ydif3=yw(iwm,3,iwa)-yw(iwm,3,jwa)
      ydist=sqrt(ydif1**2+ydif2**2+ydif3**2)
      ep3sum=ep3sum+(swt(iwm,1,iwb)*ydif1+swt(iwm,2,iwb)*ydif2+
     + swt(iwm,3,iwb)*ydif3)*cdwb3/ydist
   19 continue
   18 continue
      ep3(1)=ep3(1)+ep3sum*tstepi*tstepi
      endif
c
      return
      end
c $Id$
