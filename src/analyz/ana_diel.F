      subroutine ana_diel(xs,qs,isel,xw,qw,iwsel)
c
c $Id$
c
      implicit none
c
#include "ana_params.fh"
#include "ana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      real*8 xs(msa,3),qs(msa),xw(mwm,mwa,3),qw(mwa)
      integer isel(msa),iwsel(mwm)
c
      integer i,j,nxp,nxm,nxpw,nxmw
      real*8 xp(3),xm(3),qp,qm
c
      do 1 i=1,3
      xp(i)=0.0d0
      xm(i)=0.0d0
    1 continue
      qp=0.0d0
      qm=0.0d0
      nxp=0
      nxm=0
      nxpw=0
      nxmw=0
c
      do 2 i=1,nsa
      if(isel(i).gt.0) then
      if(qs(i).lt.0.0d0) then
      nxm=nxm+1
      xm(1)=xm(1)+qs(i)*xs(i,1)
      xm(2)=xm(2)+qs(i)*xs(i,2)
      xm(3)=xm(3)+qs(i)*xs(i,3)
      qm=qm+qs(i)
      else
      nxp=nxp+1
      xp(1)=xp(1)+qs(i)*xs(i,1)
      xp(2)=xp(2)+qs(i)*xs(i,2)
      xp(3)=xp(3)+qs(i)*xs(i,3)
      qp=qp+qs(i)
      endif
      endif
    2 continue
c
      do 3 i=1,nwm
      do 4 j=1,nwa
      if(iwsel(i).gt.0) then
      if(qw(j).lt.0.0d0) then
      nxmw=nxmw+1
      xm(1)=xm(1)+qw(j)*xw(i,j,1)
      xm(2)=xm(2)+qw(j)*xw(i,j,2)
      xm(3)=xm(3)+qw(j)*xw(i,j,3)
      qm=qm+qw(j)
      else
      nxpw=nxpw+1
      xp(1)=xp(1)+qw(j)*xw(i,j,1)
      xp(2)=xp(2)+qw(j)*xw(i,j,2)
      xp(3)=xp(3)+qw(j)*xw(i,j,3)
      qp=qp+qw(j)
      endif
      endif
    4 continue
    3 continue
c
      write(lfndie,1000) nxp,nxm,nxpw,nxmw,qp,qm,(xp(i)+xm(i),i=1,3)
 1000 format(4i5,2f10.3,6f12.3)
c
      return
      end
      subroutine ana_diehdr(sgmnam)
c
      implicit none
c
#include "ana_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      character*16 sgmnam(nsa)
c
      return
      end
