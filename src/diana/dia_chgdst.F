      subroutine dia_chgdst(hist,xs,qs)
c
      implicit none
c
#include "dia_common.fh"
c
      real*8 hist(2,ldist),xs(msa,3),qs(msa)
c
      integer i,j
      real*8 range,drange
c
      do 1 i=1,ldist
      hist(1,i)=0.0d0
      hist(2,i)=0.0d0
    1 continue
c
      range=xsmax(3)-xsmin(3)
      do 2 i=1,nsa
      j=int((dble(ldist)*(xs(i,3)-xsmin(3)))/range)+1
      if(j.gt.0.and.j.le.ldist) then
      if(qs(i).lt.0.0d0) then
      hist(1,j)=hist(1,j)+qs(i)
      else
      hist(2,j)=hist(2,j)+qs(i)
      endif
      endif
    2 continue
c
      drange=range/dble(ldist)
      if(me.eq.0) then
      write(lfnhis,1000) time,ldist
 1000 format(f12.6,i5)
      write(lfnhis,1001) (xsmin(3)+(dble(i)-0.5d0)*drange,
     + hist(1,i),hist(2,i),hist(1,i)+hist(2,i),i=1,ldist)
 1001 format(4f12.6)
      endif
c
      return
      end
c $Id$
