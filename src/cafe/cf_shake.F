      subroutine cf_shakew(xw,yw,iwgm,iwdt,numw)
c $Id: cf_shake.F,v 1.10 2000-02-18 07:02:30 d3j191 Exp $
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer numw
      real*8 xw(mwm,3,mwa),yw(mwm,3,mwa)
      integer iwgm(mwm),iwdt(mwm)
c
      call cf_shkw(numw,iwgm,iwdt,xw,yw,
     + mbt(1),numb(1),mbp(1),int_mb(i_ibnd(1)),dbl_mb(i_bnd(1)),
     + int_mb(i_iwa),dbl_mb(i_mas),
     + dbl_mb(i_s3r1),dbl_mb(i_s1r1),dbl_mb(i_s1r2),dbl_mb(i_s1r3),
     + dbl_mb(i_s1r4),dbl_mb(i_smr1))
c
      return
      end
      subroutine cf_shakes(lhandl,xs,ys,igan,isat,isdt,isfr,isto)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer lhandl,isfr,isto
      real*8 xs(msa,3),ys(msa,3)
      integer igan(msa),isat(msa),isdt(msa)
c
      integer lptr,nbonds
c
      lptr=int_mb(i_list+lhandl+4)
      nbonds=int_mb(lptr)
      if(nbonds.gt.0) then
      call cf_shks(nbonds,int_mb(lptr+1),mbt(2),mbp(2),
     + int_mb(i_ibnd(2)),dbl_mb(i_bnd(2)),dbl_mb(i_rbnd(2)),
     + isto,msa,igan,isat,isdt,xs,ys,dbl_mb(i_mas))
      endif
c
      return
      end
      logical function cf_shakep(xs,ys,igan,isat,isdt,nums)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer nums
      real*8 xs(msa,3),ys(msa,3)
      integer igan(msa),isat(msa),isdt(msa)
      integer nitall
c
      cf_shakep=.true.
c
      if(lpmf.and.lpmfc) then
      call cf_shkp(nums,xs,ys,igan,isat,int_mb(i_ipmf),int_mb(i_jpmf),
     + dbl_mb(i_rpmf),dbl_mb(i_xpmf),dbl_mb(i_ypmf),dbl_mb(i_mas))
      nitall=niters-1
      call ga_igop(mcf_64,nitall,1,'+')
      cf_shakep=nitall.le.0
      endif
c
      return
      end
      subroutine cf_shkp(nums,xs,ys,isgan,isat,
     + ipmf,jpmf,rpmf,xpmf,ypmf,wgt)
c
      implicit none
c
#include "cf_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer nums
      real*8 xs(msa,3),ys(msa,3)
      integer isgan(msa),isat(msa)
      integer ipmf(numpmf,5),jpmf(numpmf,4,npmfa)
      real*8 rpmf(numpmf,2,6),xpmf(numpmf,4,4),ypmf(numpmf,4,3)
      real*8 wgt(mat,6),ft(4,3),rv(2),w(4)
      integer ip(3),it(2),ign(4),iml(4),idt(4),ich(4)
c
      integer i,j,k,l
c
      epmf=zero
c
      do 1 i=1,numpmf
      do 2 j=1,4
      xpmf(i,j,1)=zero
      xpmf(i,j,2)=zero
      xpmf(i,j,3)=zero
      xpmf(i,j,4)=zero
      ypmf(i,j,1)=zero
      ypmf(i,j,2)=zero
      ypmf(i,j,3)=zero
    2 continue
      do 3 j=1,ipmf(i,1)
      do 4 k=1,ipmf(i,j+1)
      do 5 l=1,nums
      if(isgan(l).eq.jpmf(i,j,k)) then
      xpmf(i,j,1)=xpmf(i,j,1)+xs(l,1)
      xpmf(i,j,2)=xpmf(i,j,2)+xs(l,2)
      xpmf(i,j,3)=xpmf(i,j,3)+xs(l,3)
      xpmf(i,j,4)=xpmf(i,j,4)+wgt(isat(l),iset)
      ypmf(i,j,1)=ypmf(i,j,1)+ys(l,1)
      ypmf(i,j,2)=ypmf(i,j,2)+ys(l,2)
      ypmf(i,j,3)=ypmf(i,j,3)+ys(l,3)
      endif
    5 continue
    4 continue
      xpmf(i,j,1)=xpmf(i,j,1)/dble(ipmf(i,j+1))
      xpmf(i,j,2)=xpmf(i,j,2)/dble(ipmf(i,j+1))
      xpmf(i,j,3)=xpmf(i,j,3)/dble(ipmf(i,j+1))
      ypmf(i,j,1)=ypmf(i,j,1)/dble(ipmf(i,j+1))
      ypmf(i,j,2)=ypmf(i,j,2)/dble(ipmf(i,j+1))
      ypmf(i,j,3)=ypmf(i,j,3)/dble(ipmf(i,j+1))
      do 14 k=1,ipmf(i,j+1)
      do 15 l=1,nums
      if(isgan(l).eq.jpmf(i,j,k)) then
      xs(l,1)=xs(l,1)-xpmf(i,j,1)
      xs(l,2)=xs(l,2)-xpmf(i,j,2)
      xs(l,3)=xs(l,3)-xpmf(i,j,3)
      endif
   15 continue
   14 continue
    3 continue
    1 continue
c
      call ga_dgop(mcf_59,xpmf,12*numpmf,'+')
      call ga_dgop(mcf_63,ypmf,12*numpmf,'+')
c
      do 6 i=1,numpmf
c
      if(ipmf(i,1).eq.2) then
      do 7 j=1,2
      w(j)=xpmf(i,j,4)
      ign(j)=j
      iml(j)=j
      idt(j)=0
      ich(j)=1
      do 8 k=1,3
      ft(j,k)=zero
    8 continue
    7 continue
      it(1)=1
      ip(1)=1
      ip(2)=2
      ip(3)=1
      call cf_shks(1,it,1,2,ip,rpmf(i,1,1),rv,2,4,ign,iml,idt,
     + xpmf(i,1,1),ypmf(i,1,1),w)
      do 9 j=1,ipmf(i,1)
      do 10 k=1,ipmf(i,j+1)
      do 11 l=1,nums
      if(isgan(l).eq.jpmf(i,j,k)) then
      xs(l,1)=xs(l,1)+xpmf(i,j,1)
      xs(l,2)=xs(l,2)+xpmf(i,j,2)
      xs(l,3)=xs(l,3)+xpmf(i,j,3)
      endif
   11 continue
   10 continue
    9 continue
c
      endif
c      if(ipmf(i,1).eq.3) then
c      call errquit('Angle pmf not implemented',0)
c      endif
c      if(ipmf(i,1).eq.4) then
c      call errquit('Torsion pmf not implemented',0)
c      endif
    6 continue
c
      return
      end
      subroutine cf_shkw(nwmloc,iwgm,iwdt,xw,yw,mwb,nwb,nbp,ibnd,bnd,
     + iwat,wgt,rwx,rwi1,rwi2,rwi6,f,swt)
c
      implicit none
c
#include "cf_common.fh"
c
      integer nwmloc,iwgm(mwm),iwdt(mwm)
      real*8 xw(mwm,3,mwa),yw(mwm,3,mwa)
      integer mwb,nwb,nbp
      integer ibnd(mwb,3),iwat(mwa)
      real*8 bnd(mwb,nbp,6),wgt(mat,6)
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
     + call errquit('Deviation too large for solvent',ndtl)
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
     + call errquit('Too many shake iterations',niterw)
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
      subroutine cf_shks(nbonds,indexl,msb,msp,ibnd,bnd,rbnd,
     + natoms,matoms,igan,isat,isdt,xs,ys,wgt)
c
      implicit none
c
#include "cf_common.fh"
c
      integer msb,msp
      integer ibnd(msb,3)
      real*8 bnd(msb,msp,6),rbnd(msb,2)
      integer nbonds
      integer indexl(nbonds)
      integer matoms,natoms
      integer igan(matoms),isat(matoms),isdt(matoms)
      real*8 xs(matoms,3),ys(matoms,3)
      real*8 wgt(mat,6)
c
      real*8 wsai,wsaj
      logical ready,dtl
      integer i,j,isb,isa,jsa,isbmax,ndtl
      integer isaglo,jsaglo
      real*8 dmax,bs,ts,w,rsx1,rsx2,rsx3,rsi1,dif1,dif2,dif3,rsi2
      real*8 rsi6,ff1,ff2,ff3
c
#include "bitops.fh"
c
      niters=0
c
c     main iterative loop
c
    1 continue
      ready=.true.
      dtl=.false.
      ndtl=0
      dmax=zero
      isbmax=0
c
      do 2 i=1,nbonds
c
      isb=indexl(i)
      if(iand(ibnd(isb,3),icnstr).ne.0) then
c
      isa=0
      jsa=0
      do 3 j=1,natoms
      if(ibnd(isb,1).eq.igan(j)) isa=j
      if(ibnd(isb,2).eq.igan(j)) jsa=j
    3 continue
c
      rsx1=xs(isa,1)-xs(jsa,1)
      rsx2=xs(isa,2)-xs(jsa,2)
      rsx3=xs(isa,3)-xs(jsa,3)
      rbnd(isb,1)=sqrt(rsx1*rsx1+rsx2*rsx2+rsx3*rsx3)
      rbnd(isb,2)=zero
c
      if(iand(isdt(isa),mfixed).ne.lfixed.and.
     + iand(isdt(jsa),mfixed).ne.lfixed) then
c
      bs=bnd(isb,1,iset)*bnd(isb,1,iset)
      ts=bs*tolsha
      wsai=one/wgt(isat(isa),iset)
      wsaj=one/wgt(isat(jsa),iset)
      w=wsai+wsaj
      rsi1=bs-rsx1**2-rsx2**2-rsx3**2
      if(abs(rsi1).gt.dmax) then
      dmax=abs(rsi1)
      isbmax=isb
      endif
      if(abs(rsi1).gt.ts) then
      ready=.false.
      dif1=ys(isa,1)-ys(jsa,1)
      dif2=ys(isa,2)-ys(jsa,2)
      dif3=ys(isa,3)-ys(jsa,3)
      rsi2=dif1*rsx1+dif2*rsx2+dif3*rsx3
      rsx1=dif1
      rsx2=dif2
      rsx3=dif3
      if(rsi2.lt.small) then
      rsi2=small
      dtl=.true.
      ndtl=isb
      endif
      rsi6=half*rsi1/(rsi2*w)
      ff1=rsi6*rsx1
      ff2=rsi6*rsx2
      ff3=rsi6*rsx3
      if(iand(isdt(isa),mdynam).eq.ldynam.or.
     + iand(isdt(isa),mrestr).eq.lrestr) then
      xs(isa,1)=xs(isa,1)+ff1*wsai
      xs(isa,2)=xs(isa,2)+ff2*wsai
      xs(isa,3)=xs(isa,3)+ff3*wsai
      endif
      if(iand(isdt(jsa),mdynam).eq.ldynam.or.
     + iand(isdt(jsa),mrestr).eq.lrestr) then
      xs(jsa,1)=xs(jsa,1)-ff1*wsaj
      xs(jsa,2)=xs(jsa,2)-ff2*wsaj
      xs(jsa,3)=xs(jsa,3)-ff3*wsaj
      endif
      if(ithint) then
      deriv(19,1)=deriv(19,1)+tstepi*tstepi*
     + (rsx1*ff1+rsx2*ff2+rsx3*ff3)*bnd(isb,1,4)/
     + sqrt(rsx1*rsx1+rsx2*rsx2+rsx3*rsx3)
      endif
      endif
      endif
      endif
    2 continue
c
      niters=niters+1
c
      if(niters.gt.mshits) call errquit('Too many iterations',0)
c
      if(dtl) then
      isaglo=ibnd(ndtl,1)
      jsaglo=ibnd(ndtl,2)
      call errquit('Deviation too large solute',isaglo)
      endif
c
      if(.not.ready) goto 1
c
      return
      end

