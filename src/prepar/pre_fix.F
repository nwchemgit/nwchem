      logical function pre_fix(lfnout,ifix,jfix,rfix,cfix,
     + iwmr,xw,mwm,mwa,nwm,nwa,iwmrc,xwc,mwmc,nwmc,
     + isar,isgm,csa,xs,msa,nsa)
c
c $Id: pre_fix.F,v 1.1 1998-07-09 20:20:57 d3j191 Exp $
c
      implicit none
c
#include "util.fh"
c
      integer lfnout
      integer mwm,mwa,nwm,nwa,msa,nsa,mwmc,nwmc
      integer ifix,jfix
      real*8 rfix
      character*80 cfix
c
      real*8 xw(3,mwa,mwm),xs(3,msa),xwc(3,mwa,mwmc)
      integer iwmrc(mwmc),iwmr(mwm),isar(msa),isgm(msa)
      character*16 csa(msa)
c
      integer iatom,isgmnt,jsgmnt,ndx,ndy
      integer i,j,ifrst,ilast,jfrst,jlast
      character*6 aname
      logical lsin
      real*8 dx,dy,dz,d2,rfix2
c
      rfix2=rfix*rfix
c
c     determine center
c
      isgmnt=0
      aname='      '
      ndx=index(cfix,':')
      ndy=index(cfix,' ')
      if(ndx.gt.0) then
      read(cfix(1:ndx-1),*) isgmnt
      if(ndy-ndx.gt.5)
     + call errquit('pre_fix: atom name problem',9999)
      aname(1:ndy-ndx)=cfix(ndx+1:ndy)
      else
      if(ndy.gt.5) call errquit('pre_modify: atom name problem',9999)
      aname(1:ndy-1)=cfix(1:ndy-1)
      endif
      if(aname(1:1).eq.'_') aname(1:1)=' '
c
      if(util_print('where',print_high)) then
      write(lfnout,1000) ifix,jfix,rfix,isgmnt,aname
 1000 format('FIX ',2i5,f12.6,1x,i5,':',a)
      endif
c
      iatom=0
      ifrst=0
      ilast=0
      do 1 i=1,nsa
      if(isgm(i).eq.isgmnt) then
      if(csa(i)(11:16).eq.aname) iatom=i
      if(ifrst.eq.0) ifrst=i
      ilast=i
      endif
    1 continue
      if(iatom.eq.0.or.ifrst.eq.0.or.ilast.eq.0)
     + call errquit('pre_fix: center not found',9999)
c
      jsgmnt=isgm(1)
      jfrst=1
      jlast=1
      lsin=.false.
      do 3 j=1,nsa
      if(jsgmnt.ne.isgm(j).or.j.eq.nsa) then
      if((ifix.gt.0.and..not.lsin).or.(ifix.lt.0.and.lsin)) then
      do 5 i=jfrst,jlast
      isar(i)=1
    5 continue
      endif
      jsgmnt=isgm(j)
      jfrst=j
      lsin=.false.
      else
      jlast=j
      endif
      if(.not.lsin) then
      do 6 i=ifrst,ilast
      dx=abs(xs(1,i)-xs(1,j))
      if(dx.le.rfix) then
      lsin=.true.
      else
      dy=abs(xs(1,i)-xs(1,j))
      if(dx.le.rfix) then
      lsin=.true.
      else
      dz=abs(xs(1,i)-xs(1,j))
      if(dx.le.rfix) then
      lsin=.true.
      else
      d2=dx*dx+dy*dy+dz*dz
      if(d2.le.rfix2) lsin=.true.
      endif
      endif
      endif
    6 continue
      endif
    3 continue
c
      pre_fix=.true.
      return
c
 9999 continue
      pre_fix=.false.
      return
      end
