      subroutine ana_rama(x,iram)
c
c $Id: ana_rama.F,v 1.1 2000-09-14 21:38:22 d3j191 Exp $
c
      implicit none
c
#include "ana_params.fh"
#include "ana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      real*8 ana_torsion
      external ana_torsion
c
      real*8 x(msa,3)
      integer iram(msgm,7)
c
      real*8 phi,psi,rhb1,rhb2,dist
      integer i,j,ihb1,ihb2,k
      logical lphi,lpsi
c
      if(me.eq.0) then
      write(lfnram,1000) time
 1000 format(f12.6)
      endif
c
      rhb1=0.0d0
      rhb2=0.0d0
c
      do 1 i=1,nsgm
      lphi=iram(i,1).gt.0
      lpsi=iram(i,5).gt.0
      phi=0.0d0
      psi=0.0d0
      if(lphi) phi=ana_torsion(dbl_mb(i_xdat),msa,
     + dbl_mb(i_wdat),mwm,mwa,
     + iram(i,1),0,iram(i,2),0,iram(i,3),0,iram(i,4),0,)
      if(lpsi) psi=ana_torsion(dbl_mb(i_xdat),msa,
     + dbl_mb(i_wdat),mwm,mwa,
     + iram(i,2),0,iram(i,3),0,iram(i,4),0,iram(i,5),0,)
      ihb1=0
      ihb2=0
c
      if(iram(i,6).gt.0) then
      do 2 j=1,nsgm
      if(i.ne.j.and.iram(j,7).gt.0) then
      dist=(x(iram(i,6),1)-x(iram(j,7),1))**2+
     + (x(iram(i,6),2)-x(iram(j,7),2))**2+
     + (x(iram(i,6),3)-x(iram(j,7),3))**2
      if(ihb1.eq.0.or.dist.lt.rhb1) then
      ihb1=j
      rhb1=dist
      endif
      endif
    2 continue
      endif
c
      if(iram(i,7).gt.0) then
      do 3 j=1,nsgm
      if(i.ne.j.and.iram(j,6).gt.0) then
      dist=(x(iram(i,7),1)-x(iram(j,6),1))**2+
     + (x(iram(i,7),2)-x(iram(j,6),2))**2+
     + (x(iram(i,7),3)-x(iram(j,6),3))**2
      if(ihb2.eq.0.or.dist.lt.rhb2) then
      ihb2=j
      rhb2=dist
      endif
      endif
    3 continue
      endif
c
      if(ihb1.gt.0) then
      rhb1=sqrt(rhb1)
      if(rhb1.lt.rhbmin.or.rhb1.gt.rhbmax) then
      ihb1=0
      rhb1=0.0d0
      endif
      endif
c
      if(ihb2.gt.0) then
      rhb2=sqrt(rhb2)
      if(rhb2.lt.rhbmin.or.rhb2.gt.rhbmax) then
      ihb2=0
      rhb2=0.0d0
      endif
      endif
c
      if(me.eq.0) then
      if(lphi.or.lpsi) then
      write(lfnram,1001) i,lphi,lpsi,phi,psi,ihb1,ihb2,rhb1,rhb2
 1001 format(i7,1x,2l1,2f12.6,2i7,2f12.6)
      endif
      endif
    1 continue
c
      if(me.eq.0) then
      write(lfnram,1002) 0
 1002 format(i7)
      endif
c
      return
      end
