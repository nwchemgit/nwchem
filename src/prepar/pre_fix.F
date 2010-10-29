      logical function pre_fix(lfnout,ifix,jfix,rfix,cfix,
     + iwmr,xw,mwm,mwa,nwm,nwa,iwmrc,xwc,mwmc,nwmc,
     + isar,isgm,csa,xs,msa,nsa)
c
c $Id$
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
      integer i,j,ifrst,ilast,jfrst
      character*6 aname
      logical lsin
      real*8 dx,dy,dz,d2,rfix2
      character*4 csgmnt,string
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
      csgmnt='    '
      read(cfix(1:ndx-1),*,err=99) isgmnt
      goto 999
   99 continue
      isgmnt=0
      if(ndx.le.4) then
      csgmnt(1:ndx-1)=cfix(1:ndx-1)
      else
      csgmnt=cfix(1:4)
      endif
  999 continue
      if(ndy-ndx.gt.5)
     + call md_abort('pre_fix: atom name problem',9999)
      aname(1:ndy-ndx)=cfix(ndx+1:ndy)
      else
      if(ndy.gt.5) call md_abort('pre_modify: atom name problem',9999)
      aname(1:ndy-1)=cfix(1:ndy-1)
      endif
      do 1 i=1,index(aname,' ')-1
      if(aname(i:i).eq.'_') aname(i:i)=' '
    1 continue
c
  100 continue
c
      if(csgmnt.ne.'    ') then
      do 101 i=1,nsa
      string=csa(i)(1:4)
      if(csgmnt.eq.string.and.isgmnt.lt.isgm(i)) then
      isgmnt=isgm(i)
      goto 102
      endif
  101 continue
      goto 199
      endif
  102 continue
c
c     determine iatom : atom number of specified atom
c               ifrst : number first atom of the segment
c               ilast : number last atom of the segment
c
      iatom=0
      ifrst=0
      ilast=0
      do 2 i=1,nsa
      if(isgm(i).eq.isgmnt) then
      if(csa(i)(11:16).eq.aname) iatom=i
      if(ifrst.eq.0) ifrst=i
      ilast=i
      endif
    2 continue
c
      if(iatom.eq.0.or.ifrst.eq.0.or.ilast.eq.0) then
      if(util_print('restart',print_low)) then
      write(lfnout,'(a,i8,a,a)') 'Unable to find ',isgmnt,':',aname
      endif
      call md_abort('pre_fix: center not found',9999)
      endif
c
      if(util_print('restart',print_high)) then
      write(lfnout,1000) ifix,jfix,rfix,isgmnt,aname,ifrst,ilast
 1000 format(' fixed ',2i5,f12.6,1x,i5,':',a,5x,2i7)
      endif
c
c     fix solvent molecules
c     ---------------------
c
      do 3 i=1,nwm
      do 4 j=1,nwa
      dx=abs(xw(1,j,i)-xs(1,iatom))
      dy=abs(xw(2,j,i)-xs(2,iatom))
      dz=abs(xw(3,j,i)-xs(3,iatom))
      d2=dx*dx+dy*dy+dz*dz
      if(ifix.eq.1) then
      if(d2.ge.rfix2) iwmr(i)=2
      else
      if(d2.le.rfix2) iwmr(i)=2
      endif
    4 continue
    3 continue
      do 5 i=1,nwmc
      do 6 j=1,nwa
      dx=abs(xwc(1,j,i)-xs(1,iatom))
      dy=abs(xwc(2,j,i)-xs(2,iatom))
      dz=abs(xwc(3,j,i)-xs(3,iatom))
      d2=dx*dx+dy*dy+dz*dz
      if(ifix.eq.1) then
      if(d2.ge.rfix2) iwmrc(i)=2
      else
      if(d2.le.rfix2) iwmrc(i)=2
      endif
    6 continue
    5 continue
c
c     fix atoms
c     ---------
c
      if(jfix.eq.1) then
      do 7 i=1,nsa
      dx=abs(xs(1,i)-xs(1,iatom))
      dy=abs(xs(2,i)-xs(2,iatom))
      dz=abs(xs(3,i)-xs(3,iatom))
      d2=dx*dx+dy*dy+dz*dz
      if(ifix.eq.1) then
      if(d2.ge.rfix2) isar(i)=2
      else
      if(d2.le.rfix2) isar(i)=2
      endif
    7 continue
      else
c
c     fix segments
c     ------------
c
      jfrst=1
      do 8 i=1,nsa
      if(i.eq.nsa) then
      jsgmnt=0
      else
      jsgmnt=isgm(i+1)
      endif
c
      if(jsgmnt.ne.isgm(i)) then
      lsin=.false.
      do 9 j=jfrst,i
      dx=abs(xs(1,j)-xs(1,iatom))
      dy=abs(xs(2,j)-xs(2,iatom))
      dz=abs(xs(3,j)-xs(3,iatom))
      d2=dx*dx+dy*dy+dz*dz
      if(d2.le.rfix2) lsin=.true.
    9 continue
      if((lsin.and.ifix.eq.-1).or.(.not.lsin.and.ifix.eq.1)) then
      do 10 j=jfrst,i
      isar(j)=2
   10 continue
      endif
      jfrst=i+1
      endif
    8 continue
c
      endif
c
      if(csgmnt.ne.'    ') goto 100
  199 continue
c
      pre_fix=.true.
      return
c
 9999 continue
      pre_fix=.false.
      return
      end
