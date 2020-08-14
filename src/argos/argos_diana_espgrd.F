      subroutine argos_diana_espgrd(x,q,t,g,lprint,isel)
c
      implicit none
c
#include "argos_diana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      real*8 x(nsa,3),q(nsa,nq),g(ngx,ngy,ngz)
      character*2 t(nsa)
      logical lprint
      integer isel(nsa)
c
      integer ix,iy,iz,i,j,k,nxfr,nxto
      real*8 xg,yg,zg,value
      real*8 dx,dy,dz
      real*8 boxh(3)
c
      mgrd=ngx*ngy*ngz
      lplt=1
c
      if(np.eq.1) then
      nxfr=1
      nxto=ngx
      else
      nxfr=me*(ngx/np)+1
      nxto=(me+1)*(ngx/np)
      if(me.eq.np-1) nxto=ngx
      if(nxto.gt.ngx) nxto=ngx
      endif
c
      boxh(1)=0.5d0*box(1)
      boxh(2)=0.5d0*box(2)
      boxh(3)=0.5d0*box(3)
c
      do 1 i=1,ngx
      do 2 j=1,ngy
      do 3 k=1,ngz
      g(i,j,k)=0.0d0
    3 continue
    2 continue
    1 continue
c
      do 4 ix=nxfr,nxto
      xg=dble(ix-1)*spac(1)+xmin(1)
      do 5 iy=1,ngy
      yg=dble(iy-1)*spac(2)+xmin(2)
      do 6 iz=1,ngz
      zg=dble(iz-1)*spac(3)+xmin(3)
      do 7 i=1,nsa
      if(isel(i).gt.0) then
      dx=xg-x(i,1)
      dy=yg-x(i,2)
      dz=zg-x(i,3)
      if(iesppb.eq.2.or.iesppb.eq.3) then
      if(abs(dx).gt.boxh(1)) then
      dx=dx-anint(dx/box(1))*box(1)
      endif
      if(abs(dy).gt.boxh(2)) then
      dy=dy-anint(dy/box(2))*box(2)
      endif
      if(iesppb.eq.1.or.iesppb.eq.3) then
      if(abs(dz).gt.boxh(3)) then
      dz=dz-anint(dz/box(3))*box(3)
      endif
      endif
      endif
      g(ix,iy,iz)=g(ix,iy,iz)+138.9354d0*
     + q(i,lplt)/sqrt(dx**2+dy**2+dz**2)
      endif
    7 continue
    6 continue
    5 continue
    4 continue
c
c     global sum
c
      if(np.gt.1) then
      call ga_dgop(mag_d08,g,mgrd,'+')
      endif
c
c     write to plt file
c
      if(me.eq.0) then
      open(unit=lfnplt,file=filplt(1:index(filplt,' ')-1),
     + form='formatted',status='unknown')
      rewind(lfnplt)
      write(lfnplt,1000) 3,3
 1000 format(2i7)
      write(lfnplt,1001) ngz,ngy,ngx
 1001 format(3i7)
      write(lfnplt,1002)
     + 10.0d0*xmin(3),10.0d0*(xmin(3)+dble(ngz-1)*spac(3)),
     + 10.0d0*xmin(2),10.0d0*(xmin(2)+dble(ngy-1)*spac(2)),
     + 10.0d0*xmin(1),10.0d0*(xmin(1)+dble(ngx-1)*spac(1))
 1002 format(6e12.5)
      write(lfnplt,1003) (((g(ix,iy,iz),ix=1,ngx),iy=1,ngy),iz=1,ngz)
 1003 format(2e20.12)
      close(unit=lfnplt)
      if(lprint) write(*,1004) filplt(1:index(filplt,' ')-1)
 1004 format(' Grid data written to ',a)
c
      open(unit=lfnepz,file=filepz(1:index(filepz,' ')-1),
     + form='formatted',status='unknown')
      rewind(lfnepz)
      do 8 iz=1,ngz
      value=0.0d0
      do 9 iy=1,ngy
      do 10 ix=1,ngx
      value=value+g(ix,iy,iz)
   10 continue
    9 continue
      write(lfnepz,1005) xmin(3)+dble(iz)*spac(3),value/dble(ngx*ngy)
 1005 format(2f12.6)
    8 continue
      close(unit=lfnepz)
      endif
c
      return
      end
c $Id$
