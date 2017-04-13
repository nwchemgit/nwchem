      logical function argos_prepare_solvat(islv,lfnout,boxs,xslv,vslv,
     + mwms,mwa,nwms,nwa,touch,box,iwmr,xw,vw,
     + mwm,nwm,xs,msa,nsa,xwc,vwc,mwmc,nwmc)
c
c $Id$
c
      implicit none
c
#include "util.fh"
c
      integer islv,mwms,mwa,nwms,nwa,mwm,nwm,msa,nsa,lfnout
      integer nwmc,mwmc
      integer iwmr(mwm)
      real*8 xw(3,mwa,mwm),xs(3,nsa),xslv(3,mwa,mwms)
      real*8 xwc(3,mwa,mwmc),vwc(3,mwa,mwmc)
      real*8 vw(3,mwa,mwm),vslv(3,mwa,mwms)
      real*8 box(3),boxs(3),touch,dx,dy,dz
c
      integer i,j,k,is,ix,iy,iz,nx,ny,nz
      real*8 boxh(3),xc(3),touch2,x,y,z,pbc
      logical ltry
c
      if(islv.eq.0.or.nwms.le.0) goto 9999
c
      nwm=0
      touch2=touch*touch
c
      boxh(1)=0.5d0*box(1)
      boxh(2)=0.5d0*box(2)
      boxh(3)=0.5d0*box(3)
c
      nx=int(box(1)/boxs(1))+2
      ny=int(box(2)/boxs(2))+2
      nz=int(box(3)/boxs(3))+2
c
      do 4 is=1,nwms
      do 5 k=1,3
      x=0.0d0
      do 6 j=1,nwa
      x=x+xslv(k,j,is)
    6 continue
      xc(k)=x/dble(nwa)
    5 continue
c
      do 7 ix=-nx,nx
      do 8 iy=-ny,ny
      do 9 iz=-nz,nz
c
      x=xc(1)+dble(ix)*boxs(1)
      y=xc(2)+dble(iy)*boxs(2)
      z=xc(3)+dble(iz)*boxs(3)
c
      ltry=.false.
      if(islv.eq.1.or.islv.eq.3) then
      ltry=abs(x).lt.boxh(1).and.abs(y).lt.boxh(2).and.abs(z).lt.boxh(3)
      elseif(islv.eq.2) then
      ltry=(x*x+y*y+z*z).le.boxh(1)*boxh(1)
      endif
c
      if(ltry) then
      do 10 i=1,nsa
      dx=x-xs(1,i)
      dy=y-xs(2,i)
      dz=z-xs(3,i)
      dx=dx-nint(dx/box(1))*box(1)
      dy=dy-nint(dy/box(2))*box(2)
      dz=dz-nint(dz/box(3))*box(3)
      if(dx*dx+dy*dy+dz*dz.lt.touch2) goto 9
   10 continue
      if(nwmc.gt.0) then
      do 11 i=1,nwmc
      do 12 j=1,nwa
      dx=x-xwc(1,j,i)
      dy=y-xwc(2,j,i)
      dz=z-xwc(3,j,i)
      dx=dx-nint(dx/box(1))*box(1)
      dy=dy-nint(dy/box(2))*box(2)
      dz=dz-nint(dz/box(3))*box(3)
      if(dx*dx+dy*dy+dz*dz.lt.touch2) goto 9
   12 continue
   11 continue
      endif
      if(islv.eq.3) then
      pbc=aint(4.0d0*
     + (abs(x)/box(1)+abs(y)/box(2)+abs(z)/box(3))/3.0d0)
      if(pbc.ne.0.0d0) goto 9
      endif
      if(nwm.ge.mwm) 
     + call md_abort('argos_prepare_solvat: increase mwm',9999)
      nwm=nwm+1
      do 13 j=1,nwa
      xw(1,j,nwm)=xslv(1,j,is)+dble(ix)*boxs(1)
      xw(2,j,nwm)=xslv(2,j,is)+dble(iy)*boxs(2)
      xw(3,j,nwm)=xslv(3,j,is)+dble(iz)*boxs(3)
      iwmr(nwm)=0
      do 14 k=1,3
      vw(k,j,nwm)=vslv(k,j,is)
   14 continue
   13 continue
      endif
c
    9 continue
    8 continue
    7 continue
    4 continue
c
      if(util_print('restart',print_medium)) then
      if(islv.eq.1) write(lfnout,1000) box
 1000 format(/,' Solvated in box',t40,3f12.6)
      if(islv.eq.2) write(lfnout,1001) 0.5d0*box(1)
 1001 format(/,' Solvated in sphere with radius',t40,f12.6)
      if(islv.eq.3) write(lfnout,1002) box
 1002 format(/,' Solvated in truncated octahedron with size',t40,3f12.6)
      endif
c
      argos_prepare_solvat=.true.
      return
c
 9999 argos_prepare_solvat=.false.
      return
      end
 
