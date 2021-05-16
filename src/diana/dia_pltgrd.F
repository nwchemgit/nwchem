      subroutine dia_espplt()
c
c $Id$
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
c
      string=card(8:80)
      read(string(1:4),'(i4)') ngrid
      read(string(5:8),'(i4)') iesppb
      read(string(9:20),'(f12.6)') rcut
      string=card(28:80)
      i=index(string,' ')
      filchg=string(1:i)
      filplt=string(i+1:80)
      filepz=filplt(1:index(filplt,'.'))//'epz '
c
      if(filchg(1:1).ne.' ') call dia_pltgrd(.true.,int_mb(i_isel))
c
      return
      end
      subroutine dia_pltgrd(lprint,isel)
c
c $Id$
c
      implicit none
c
#include "dia_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      logical lprint
      integer isel(nsa)
c
      ngrd(1)=ngrid
      ngrd(2)=ngrid
      ngrd(3)=ngrid
c
c
c     open the charge file
c
      if(me.eq.0) then
      open(unit=lfnchg,file=filchg(1:index(filchg,' ')-1),
     + form='formatted',status='old',err=9999)
c
      read(lfnchg,1000) nsa,nq,box
 1000 format(i7,i3,3f12.6)
      endif
      if(np.gt.1) then
      call ga_brdcst(mag_d01,nsa,ma_sizeof(mt_int,1,mt_byte),0)
      call ga_brdcst(mag_d01,nq,ma_sizeof(mt_int,1,mt_byte),0)
      call ga_brdcst(mag_d01,box,3*ma_sizeof(mt_dbl,1,mt_byte),0)
      endif
c
c     allocate memory for coordinates and charges
c
      if(.not.ma_push_get(mt_dbl,3*nsa,'x',l_x,i_x))
     + call md_abort('Failed to allocate memory for x',me)
      if(.not.ma_push_get(mt_dbl,nq*nsa,'q',l_q,i_q))
     + call md_abort('Failed to allocate memory for q',me)
      if(.not.ma_push_get(mt_byte,2*nsa,'t',l_t,i_t))
     + call md_abort('Failed to allocate memory for t',me)
c
c     read coordinates and charges
c
      call dia_rdchg(dbl_mb(i_x),dbl_mb(i_q),byte_mb(i_t))
c
c     allocate memory for the grid
c
      ngx=ngrd(1)
      ngy=ngrd(2)
      ngz=ngrd(3)
      mgrd=ngx*ngy*ngz
c
      if(me.eq.0) then
      if(lprint) write(*,2000) ngrd,rcut,
     + filchg(1:index(filchg,' ')-1),
     + filplt(1:index(filplt,' ')-1)
 2000 format(' Electrostatic potential on ',i5,'x',i5,'x',i5,
     + ' grid extending',
     + f12.6,' nm',/,' Charges from ',a,' written to ',a)
      if(lprint.and.iesppb.gt.0) write(*,2001)
 2001 format(' Periodic boundary conditions')
c
      endif
c
      if(.not.ma_push_get(mt_dbl,mgrd,'grid',l_grid,i_grid)) then
      call md_abort('Failed to allocate memory for grid',me)
      endif
c
c     calculate the electrostatic potential on the grid
c
      call dia_espgrd(dbl_mb(i_x),dbl_mb(i_q),byte_mb(i_t),
     + dbl_mb(i_grid),lprint,isel)
c
c     deallocate memory
c
      if(.not.ma_pop_stack(l_grid))
     + call md_abort('Failed to deallocate grid',0)
      if(.not.ma_pop_stack(l_t))
     + call md_abort('Failed to deallocate t',0)
      if(.not.ma_pop_stack(l_q))
     + call md_abort('Failed to deallocate q',0)
      if(.not.ma_pop_stack(l_x))
     + call md_abort('Failed to deallocate x',0)
c
      close(unit=lfnchg)
      return
 9999 continue
      call md_abort('dia_pltgrd: error',0)
      return
      end
      subroutine dia_rdchg(x,q,t)
c
      implicit none
c
#include "dia_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      real*8 x(nsa,3),q(nsa,nq)
      character*2 t(nsa)
c
      integer i,j
      real*8 xd(3)
c
      if(me.eq.0) then
c
      do 1 i=1,nsa
      read(lfnchg,1000) t(i),(x(i,j),j=1,3),(q(i,j),j=1,nq)
 1000 format(a2,3f12.6,4f12.6)
      x(i,1)=0.1d0*x(i,1)
      x(i,2)=0.1d0*x(i,2)
      x(i,3)=0.1d0*x(i,3)
      if(i.eq.1) then
      do 2 j=1,3
      xmin(j)=x(i,j)
      xmax(j)=x(i,j)
    2 continue
      else
      do 3 j=1,3
      if(xmin(j).gt.x(i,j)) xmin(j)=x(i,j)
      if(xmax(j).lt.x(i,j)) xmax(j)=x(i,j)
    3 continue
      endif
    1 continue
c
      do 4 j=1,3
      if(iesppb.eq.0.or.rcut.lt.0.0d0.or.
     + (iesppb.eq.1.and.j.ne.3).or.
     + (iesppb.eq.2.and.j.eq.3)) then
      xd(j)=xmax(j)-xmin(j)+2.0d0*abs(rcut)
      xmin(j)=xmin(j)-abs(rcut)
      else
      xd(j)=box(j)
      xmin(j)=0.5d0*(xmin(j)+xmax(j)-box(j))
      endif
    4 continue
c
      do 5 j=1,3
      spac(j)=1.0d0/(ngrid-1)
    5 continue
c
      do 6 j=1,3
      ngrd(j)=xd(j)/spac(j)+1
    6 continue
c
      do 7 j=1,3
      spac(j)=xd(j)/dble(ngrd(j))
    7 continue
c
      endif
c
c     broadcast to all nodes
c
      if(np.gt.1) then
      call ga_brdcst(mag_d01,x,3*nsa*ma_sizeof(mt_dbl,1,mt_byte),0)
      call ga_brdcst(mag_d02,xmin,3*ma_sizeof(mt_dbl,1,mt_byte),0)
      call ga_brdcst(mag_d03,ngrd,3*ma_sizeof(mt_int,1,mt_byte),0)
      call ga_brdcst(mag_d04,spac,3*ma_sizeof(mt_dbl,1,mt_byte),0)
      call ga_brdcst(mag_d05,q,nsa*nq*ma_sizeof(mt_dbl,1,mt_byte),0)
      endif
c
      return
      end
      subroutine dia_espgrd(x,q,t,g,lprint,isel)
c
      implicit none
c
#include "dia_common.fh"
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
