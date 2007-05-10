      subroutine dia_espplt()
c
c $Id: dia_espplt.F,v 1.12 2007-05-10 21:17:40 d3j191 Exp $
c
      implicit none
c
#include "dia_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
#include "util.fh"
c
      character*255 string
      integer i
      integer ilt,iht,jlt,jht,ndxt,ldt
c
      string=card(8:80)
      read(string(1:4),'(i4)') ngrid
      read(string(5:8),'(i4)') iesppb
      read(string(9:20),'(f12.6)') rcut
      strplt=string(21:80)
      strepz=string(21:80)
c
      call ga_distribution(ga_trj,me,il,ih,jl,jh)
      call ga_access(ga_trj,il,ih,jl,jh,ndx,ld)
      call ga_distribution(ga_trt,me,ilt,iht,jlt,jht)
      call ga_access(ga_trt,ilt,iht,jlt,jht,ndxt,ldt)
c
      call dia_pltgrd(int_mb(i_isel),dbl_mb(i_qdat),
     + dbl_mb(ndx),dbl_mb(ndxt))
c
      return
      end
      subroutine dia_pltgrd(isel,qs,xtrj,trt)
c
      implicit none
c
#include "dia_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer isel(nsa)
      real*8 qs(nsa),xtrj(ntrj,3,nfrdim),trt(6,nfrdim)
c
      integer i,j,k
      real*8 d
c
c     allocate memory for grid
c
      do 12 j=1,3
      if(iesppb.eq.0.or.rcut.lt.0.0d0.or.
     + (iesppb.eq.1.and.j.ne.3).or.
     + (iesppb.eq.2.and.j.eq.3)) then
      d=xsmax(j)-xsmin(j)+2.0d0*abs(rcut)
      else
      d=box(j)
      endif
      ngrd(j)=ngrid*int(d)
   12 continue
c
      ngx=ngrd(1)
      ngy=ngrd(2)
      ngz=ngrd(3)
      mgrd=ngx*ngy*ngz
      if(.not.ma_push_get(mt_dbl,mgrd,'grid',l_grid,i_grid))
     + call md_abort('Failed to allocate memory for grid',me)
c
      do 1 k=1,nfrme
c
      write(filplt,1000) strplt(1:index(strplt,' ')-1),
     + me*nfrdim+k,'.plt '
      write(filepz,1000) strepz(1:index(strepz,' ')-1),
     + me*nfrdim+k,'.epz '
 1000 format(a,i5.5,a)
c
c     calculate the electrostatic potential on the grid
c
      call dia_espgrd(isel,xtrj(1,1,k),trt(1,k),qs,dbl_mb(i_grid))
c
    1 continue
c
c     deallocate memory
c
      if(.not.ma_pop_stack(l_grid))
     + call md_abort('Failed to deallocate grid',0)
c
      return
 9999 continue
      call md_abort('dia_pltgrd: error',0)
      return
      end
      subroutine dia_espgrd(isel,x,t,q,g)
c
      implicit none
c
#include "dia_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
#if defined(CRAY_T3D) || defined(CRAY_T3E)
      real*8 cf_anint
      external cf_anint
#endif
c
      real*8 x(nsa,3),q(nsa),g(ngx,ngy,ngz),t(6)
      integer isel(nsa)
c
      integer ix,iy,iz,i,j,k,nxfr,nxto
      real*8 xg,yg,zg,value
      real*8 dx,dy,dz
      real*8 boxh(3),xd(3)
c
      mgrd=ngx*ngy*ngz
      lplt=1
c
      box(1)=t(4)
      box(2)=t(5)
      box(3)=t(6)
      boxh(1)=0.5d0*box(1)
      boxh(2)=0.5d0*box(2)
      boxh(3)=0.5d0*box(3)
c
      xmin(1)=x(1,1)
      xmin(2)=x(1,2)
      xmin(3)=x(1,3)
      xmax(1)=x(1,1)
      xmax(2)=x(1,2)
      xmax(3)=x(1,3)
      do 11 i=1,nsa
      xmin(1)=min(xmin(1),x(i,1))
      xmin(2)=min(xmin(2),x(i,2))
      xmin(3)=min(xmin(3),x(i,3))
      xmax(1)=max(xmax(1),x(i,1))
      xmax(2)=max(xmax(2),x(i,2))
      xmax(3)=max(xmax(3),x(i,3))
   11 continue
c
      do 12 j=1,3
      if(iesppb.eq.0.or.rcut.lt.0.0d0.or.
     + (iesppb.eq.1.and.j.ne.3).or.
     + (iesppb.eq.2.and.j.eq.3)) then
      xd(j)=xmax(j)-xmin(j)+2.0d0*abs(rcut)
      xmin(j)=xmin(j)-abs(rcut)
      else
      xd(j)=box(j)
      xmin(j)=0.5d0*(xmin(j)+xmax(j)-box(j))
      endif
   12 continue
c
      do 13 j=1,3
      spac(j)=xd(j)/(ngrd(j)-1)
   13 continue
c
      do 1 i=1,ngx
      do 2 j=1,ngy
      do 3 k=1,ngz
      g(i,j,k)=0.0d0
    3 continue
    2 continue
    1 continue
c
      do 4 ix=1,ngx
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
#if defined(CRAY_T3D) || defined(CRAY_T3E)
      dx=dx-cf_anint(dx/box(1))*box(1)
#else
      dx=dx-anint(dx/box(1))*box(1)
#endif
      endif
      if(abs(dy).gt.boxh(2)) then
#if defined(CRAY_T3D) || defined(CRAY_T3E)
      dy=dy-cf_anint(dy/box(2))*box(2)
#else
      dy=dy-anint(dy/box(2))*box(2)
#endif
      endif
      if(iesppb.eq.1.or.iesppb.eq.3) then
      if(abs(dz).gt.boxh(3)) then
#if defined(CRAY_T3D) || defined(CRAY_T3E)
      dz=dz-cf_anint(dz/box(3))*box(3)
#else
      dz=dz-anint(dz/box(3))*box(3)
#endif
      endif
      endif
      endif
      g(ix,iy,iz)=g(ix,iy,iz)+138.9354d0*
     + q(i)/sqrt(dx**2+dy**2+dz**2)
      endif
    7 continue
    6 continue
    5 continue
    4 continue
c
c     write to plt file
c
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
      write(lfnepz,1005) xmin(3)+dble(iz-1)*spac(3),value/dble(ngx*ngy)
 1005 format(2f12.6)
    8 continue
      close(unit=lfnepz)
c
      return
      end
