      subroutine argos_diana_center_fold(xs,xw,imol,nw)
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "util.fh"
c
      real*8 xs(msa,3),xw(mwm,mwa,3)
      integer imol(msa),nw
c
      integer i,j,k,nmol,ncount
      real*8 xc(3),dx,factor
c
      real*8 boxh(3)
c
      boxh(1)=half*box(1)
      boxh(2)=half*box(2)
      boxh(3)=half*box(3)
c
      xc(1)=zero
      xc(2)=zero
      xc(3)=zero
      ncount=0
c
      nmol=0
      do 1 i=1,nsa
      nmol=max(nmol,imol(i))
      if(imol(i).ge.icent.and.imol(i).le.jcent) then
      ncount=ncount+1
      xc(1)=xc(1)+xs(i,1)
      xc(2)=xc(2)+xs(i,2)
      xc(3)=xc(3)+xs(i,3)
      endif
    1 continue
      if(ncount.eq.0) call md_abort('No atom in molecule',icent)
      xc(1)=xc(1)/dble(ncount)
      xc(2)=xc(2)/dble(ncount)
      xc(3)=xc(3)/dble(ncount)
c
      do 2 i=1,nsa
      xs(i,1)=xs(i,1)-xc(1)
      xs(i,2)=xs(i,2)-xc(2)
      xs(i,3)=xs(i,3)-xc(3)
    2 continue
      do 3 j=1,nwa
      do 4 i=1,nw
      xw(i,j,1)=xw(i,j,1)-xc(1)
      xw(i,j,2)=xw(i,j,2)-xc(2)
      xw(i,j,3)=xw(i,j,3)-xc(3)
    4 continue
    3 continue
c
      do 5 j=1,nmol
      xc(1)=zero
      xc(2)=zero
      xc(3)=zero
      ncount=0
      do 6 i=1,nsa
      nmol=max(nmol,imol(i))
      if(imol(i).eq.j) then
      ncount=ncount+1
      xc(1)=xc(1)+xs(i,1)
      xc(2)=xc(2)+xs(i,2)
      xc(3)=xc(3)+xs(i,3)
      endif
    6 continue
      if(ncount.eq.0) goto 5
      xc(1)=xc(1)/dble(ncount)
      xc(2)=xc(2)/dble(ncount)
      xc(3)=xc(3)/dble(ncount)
      do 7 k=1,3
      dx=zero
      if(xc(k).lt.-boxh(k)) dx=box(k)
      if(xc(k).gt.boxh(k)) dx=-box(k)
      do 8 i=1,nsa
      if(imol(i).eq.j) xs(i,k)=xs(i,k)+dx
    8 continue
    7 continue
    5 continue
      if(nwa.gt.0.and.nw.gt.0) then
      factor=one/dble(nwa)
      do 9 i=1,nw
      xc(1)=zero
      xc(2)=zero
      xc(3)=zero
      do 10 j=1,nwa
      xc(1)=xc(1)+xw(i,j,1)
      xc(2)=xc(2)+xw(i,j,2)
      xc(3)=xc(3)+xw(i,j,3)
   10 continue
      xc(1)=xc(1)*factor
      xc(2)=xc(2)*factor
      xc(3)=xc(3)*factor
      do 11 k=1,3
      dx=zero
      if(xc(k).lt.-boxh(k)) dx=box(k)
      if(xc(k).gt.boxh(k)) dx=-box(k)
      do 12 j=1,nwa
      xw(i,j,k)=xw(i,j,k)+dx
   12 continue
   11 continue
    9 continue
      endif
c
      return
      end
c $Id$
