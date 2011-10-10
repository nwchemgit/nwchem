      subroutine argos_cafe_shks(nbonds,indexl,msb,msp,ibnd,bnd,rbnd,
     + natoms,matoms,igan,isat,isgm,isdt,ishp,xs,ys,wgt,lseq)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      integer msb,msp
      integer ibnd(msb,3)
      real*8 bnd(msb,msp,6),rbnd(msb,2)
      integer nbonds
      integer indexl(nbonds)
      integer matoms,natoms
      integer igan(matoms),isat(matoms),isdt(matoms),ishp(matoms)
      integer isgm(matoms),lseq(mseq)
      real*8 xs(matoms,3),ys(matoms,3)
      real*8 wgt(mat,mset)
c
      real*8 wsai,wsaj
      logical ready,dtl
      integer i,j,k,isb,isa,jsa,isbmax,ndtl,icset
      integer isaglo,jsaglo
      real*8 dmax,bs,ts,w,rsx1,rsx2,rsx3,rsi1,dif1,dif2,dif3,rsi2
      real*8 rsi6,ff1,ff2,ff3
c
      character*15 filerr
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
c
      isa=0
      jsa=0
      do 3 j=1,natoms
      if(ibnd(isb,1).eq.igan(j)) isa=j
      if(ibnd(isb,2).eq.igan(j)) jsa=j
    3 continue
c
      if(nfhop.eq.0) then
      icset=iset
      else
      icset=lseq(isgm(isa))
      endif
c
      if(iand(ibnd(isb,3),icnstr).ne.0.and.
     + bnd(isb,2,icset).gt.zero) then
c
      rsx1=xs(isa,1)-xs(jsa,1)
      rsx2=xs(isa,2)-xs(jsa,2)
      rsx3=xs(isa,3)-xs(jsa,3)
      rbnd(isb,1)=sqrt(rsx1*rsx1+rsx2*rsx2+rsx3*rsx3)
      rbnd(isb,2)=zero
c      write(*,'(5i5,f12.6)') isb,isa,jsa,ishp(isa),ishp(jsa),rbnd(isb,1)
c
      if(iand(isdt(isa),mfixed).ne.lfixed.or.
     + iand(isdt(jsa),mfixed).ne.lfixed) then
c
      bs=bnd(isb,1,icset)*bnd(isb,1,icset)
      ts=bs*tolsha
      wsai=one/wgt(isat(isa),icset)
      wsaj=one/wgt(isat(jsa),icset)
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
      if(iand(ishp(jsa),1).ne.1) then
      if(iand(isdt(isa),mdynam).eq.ldynam.or.
     + iand(isdt(isa),mrestr).eq.lrestr) then
      xs(isa,1)=xs(isa,1)+ff1*wsai
      xs(isa,2)=xs(isa,2)+ff2*wsai
      xs(isa,3)=xs(isa,3)+ff3*wsai
      endif
      endif
      if(iand(ishp(isa),1).ne.1) then
      if(iand(isdt(jsa),mdynam).eq.ldynam.or.
     + iand(isdt(jsa),mrestr).eq.lrestr) then
      xs(jsa,1)=xs(jsa,1)-ff1*wsaj
      xs(jsa,2)=xs(jsa,2)-ff2*wsaj
      xs(jsa,3)=xs(jsa,3)-ff3*wsaj
      endif
      endif
      if(ithint) then
      deriv(19,1)=deriv(19,1)+tstepi*tstepi*
     + (rsx1*ff1+rsx2*ff2+rsx3*ff3)*bnd(isb,1,4)/
     + sqrt(rsx1*rsx1+rsx2*rsx2+rsx3*rsx3)
      endif
      endif
      endif
c
c     place dummy hydrogens (acceptor sites) onto bound heavy atoms
c
      if(iand(ishp(isa),1).eq.1.and.iand(ishp(jsa),1).ne.1) then
      xs(isa,1)=xs(jsa,1)
      xs(isa,2)=xs(jsa,2)
      xs(isa,3)=xs(jsa,3)
      endif
      if(iand(ishp(isa),1).ne.1.and.iand(ishp(jsa),1).eq.1) then
      xs(jsa,1)=xs(isa,1)
      xs(jsa,2)=xs(isa,2)
      xs(jsa,3)=xs(isa,3)
      endif
      endif
c
    2 continue
c
      niters=niters+1
c
      if(niters.gt.mshits) call md_abort('Too many iterations',0)
c
      if(dtl) then
      isaglo=ibnd(ndtl,1)
      jsaglo=ibnd(ndtl,2)
      write(filerr,1000) me
 1000 format('shake_',i3.3,'.error')
      open(unit=16,file=filerr)
      do 4 i=1,nbonds
      isb=indexl(i)
      if(iand(ibnd(isb,3),icnstr).ne.0) then
      isa=0
      jsa=0
      do 5 j=1,natoms
      if(ibnd(isb,1).eq.igan(j)) isa=j
      if(ibnd(isb,2).eq.igan(j)) jsa=j
    5 continue
      if(isb.eq.ndtl) then
      write(16,1001) igan(isa),igan(jsa),
     + (ys(isa,k),k=1,3),(ys(jsa,k),k=1,3),
     + (xs(isa,k),k=1,3),(xs(jsa,k),k=1,3),' < '
 1001 format(2i7,12f8.3,a)
      else
      write(16,1001) igan(isa),igan(jsa),
     + (ys(isa,k),k=1,3),(ys(jsa,k),k=1,3),
     + (xs(isa,k),k=1,3),(xs(jsa,k),k=1,3)
      endif
      endif
    4 continue
      close(unit=16)
      call md_abort('Deviation too large solute',isaglo)
      endif
c
      if(.not.ready) goto 1
c
      return
      end
 
c $Id$
