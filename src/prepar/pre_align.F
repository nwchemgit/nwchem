      logical function pre_align(card,lfnout,
     + csa,isat,isgm,xs,vs,msa,nsa,xw,vw,mwm,mwa,nwm,nwa,
     + xwc,vwc,mwmc,nwmc)
c
c $Id$
c
      implicit none
c
#include "util.fh"
c
      integer msa,nsa,mwm,nwm,mwa,nwa,mwmc,nwmc,lfnout
      integer isat(msa),isgm(msa)
      character*16 csa(msa)
      character*80 card
      real*8 xs(3,msa),vs(3,msa),xw(3,mwa,mwm),vw(3,mwa,mwm)
      real*8 xwc(3,mwa,mwmc),vwc(3,mwa,mwmc)
c
      integer i,j,k,iatom,jatom,katom,iseg,jseg,kseg
      real*8 xx,rv(3),o(3),angle
      character*4 atomi,atomj,atomk
c
      read(card(8:36),'(3(i5,a4))') iseg,atomi,jseg,atomj,kseg,atomk
c
      iatom=0
      jatom=0
      katom=0
      do 1 i=1,nsa
      if(isgm(i).eq.iseg.and.csa(i)(11:14).eq.atomi) iatom=i
      if(isgm(i).eq.jseg.and.csa(i)(11:14).eq.atomj) jatom=i
      if(isgm(i).eq.kseg.and.csa(i)(11:14).eq.atomk) katom=i
    1 continue
c
c     put iatom on the origin
c
      do 2 k=1,3
      xx=xs(k,iatom)
      do 3 i=1,nsa
      xs(k,i)=xs(k,i)-xx
    3 continue
      do 4 j=1,nwa
      do 5 i=1,nwm
      xw(k,j,i)=xw(k,j,i)-xx
    5 continue
      do 6 i=1,nwmc
      xwc(k,j,i)=xwc(k,j,i)-xx
    6 continue
    4 continue
    2 continue
c
c     rotation vector
c
      rv(1)=xs(2,jatom)
      rv(2)=-(xs(1,jatom))
      rv(3)=0.0d0
      o(1)=0.0d0
      o(2)=0.0d0
      o(3)=0.0d0
c
c     rotation angle
c
      angle=acos(xs(3,jatom)/sqrt(xs(1,jatom)**2+xs(2,jatom)**2+
     + xs(3,jatom)**2))
c
c     put jatom on z axis
c
      do 7 i=1,nsa
      call rotate(o,rv,angle,xs(1,i),xs(1,i))
    7 continue
      do 8 j=1,nwa
      do 9 i=1,nwm
      call rotate(o,rv,angle,xw(1,j,i),xw(1,j,i))
    9 continue
      do 10 i=1,nwm
      call rotate(o,rv,angle,xwc(1,j,i),xwc(1,j,i))
   10 continue
    8 continue
c
      if(katom.gt.0) then
c
c     rotation vector
c
      rv(1)=0.0d0
      rv(2)=0.0d0
      rv(3)=1.0d0
      o(1)=0.0d0
      o(2)=0.0d0
      o(3)=0.0d0
c
c     rotation angle
c
      angle=-acos((xs(1,katom)+xs(2,katom))/
     + sqrt(2*(xs(1,katom)**2+xs(2,katom)**2)))
c
c     put k atom in x-y plane
c
c
c     put jatom on z axis
c
      do 11 i=1,nsa
      call rotate(o,rv,angle,xs(1,i),xs(1,i))
   11 continue
      do 12 j=1,nwa
      do 13 i=1,nwm
      call rotate(o,rv,angle,xw(1,j,i),xw(1,j,i))
   13 continue
      do 14 i=1,nwm
      call rotate(o,rv,angle,xwc(1,j,i),xwc(1,j,i))
   14 continue
   12 continue
      endif
c
      if(util_print('restart',print_medium)) then
      write(lfnout,1000) iseg,atomi,jseg,atomj,kseg,atomk
 1000 format(' Aligned atoms ',i5,':',a,i5,':',a,i5,':',a)
      endif
c
      pre_align=.true.
      return
c
 9999 continue
      pre_align=.false.
      return
      end
