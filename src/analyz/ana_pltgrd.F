      subroutine ana_pltgrd(card)
c
c $Id: ana_pltgrd.F,v 1.5 2000-01-05 00:45:02 d3j191 Exp $
c
      implicit none
c
#include "ana_common.fh"
#include "mafdecls.fh"
c
      integer mgrd,i
      character*80 card,string
c
      print*,card
      string=card(8:80)
c
      read(string(1:4),'(i4)') mgrd
      read(string(5:16),'(f12.6)') rcut
      spac=rcut/dble(mgrd-1)
      string=card(24:80)
      i=index(string,' ')
      filchg=string(1:i)
      filplt=string(i+1:80)
      print*,rcut,spac
c
      ngrd(1)=mgrd
      ngrd(2)=mgrd
      ngrd(3)=mgrd
c
      write(*,2000) mgrd,rcut,
     + filchg(1:index(filchg,' ')-1),
     + filplt(1:index(filplt,' ')-1)
 2000 format(' Electrostatic potential on ',i4,' cubed grid extending',
     + f12.6,' nm for charges from ',a,' written to ',a)
c
c     open the charge file
c
      open(unit=lfnchg,file=filchg(1:index(filchg,' ')-1),
     + form='formatted',status='old',err=9999)
c
      read(lfnchg,1000) nsa,nq
 1000 format(i7,i3)
c
c     allocate memory for coordinates and charges
c
      if(.not.ma_push_get(mt_dbl,3*nsa,'x',l_x,i_x))
     + call errquit('Failed to allocate memory for x',me)
      if(.not.ma_push_get(mt_dbl,nq*nsa,'q',l_q,i_q))
     + call errquit('Failed to allocate memory for q',me)
      if(.not.ma_push_get(mt_byte,2*nsa,'t',l_t,i_t))
     + call errquit('Failed to allocate memory for t',me)
c
c     read coordinates and charges
c
      call ana_rdchg(dbl_mb(i_x),dbl_mb(i_q),byte_mb(i_t))
c
c     allocate memory for the grid
c
      ngx=ngrd(1)
      ngy=ngrd(2)
      ngz=ngrd(3)
      mgrd=ngx*ngy*ngz
c
      if(.not.ma_push_get(mt_dbl,mgrd,'grid',l_grid,i_grid)) then
      call errquit('Failed to allocate memory for grid',me)
      endif
c
c     calculate the electrostatic potential on the grid
c
      call ana_espgrd(dbl_mb(i_x),dbl_mb(i_q),byte_mb(i_t),
     + dbl_mb(i_grid))
c
c     deallocate memory
c
      if(.not.ma_pop_stack(l_grid))
     + call errquit('Failed to deallocate grid',0)
      if(.not.ma_pop_stack(l_t))
     + call errquit('Failed to deallocate t',0)
      if(.not.ma_pop_stack(l_q))
     + call errquit('Failed to deallocate q',0)
      if(.not.ma_pop_stack(l_x))
     + call errquit('Failed to deallocate x',0)
c
      close(unit=lfnchg)
      return
 9999 continue
      call errquit('ana_pltgrd: error',0)
      return
      end
      subroutine ana_rdchg(x,q,t)
c
      implicit none
c
#include "ana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      real*8 x(nsa,3),q(nsa,nq)
      character*2 t(nsa)
c
      integer i,j
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
      ngrd(j)=(xmax(j)-xmin(j)+2.0d0*rcut)/spac+1
      xmin(j)=xmin(j)-rcut
    4 continue
c
      endif
c
c     broadcast to all nodes
c
      if(np.gt.1) then
      call ga_brdcst(mag_d01,x,3*nsa*ma_sizeof(mt_dbl,1,mt_byte),0)
      call ga_brdcst(mag_d02,xmin,3*ma_sizeof(mt_dbl,1,mt_byte),0)
      call ga_brdcst(mag_d03,ngrd,3*ma_sizeof(mt_int,1,mt_byte),0)
      endif
c
      return
      end
      subroutine ana_espgrd(x,q,t,g)
c
      implicit none
c
#include "ana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      real*8 x(nsa,3),q(nsa,nq),g(ngx,ngy,ngz)
      character*2 t(nsa)
c
      integer ix,iy,iz,i,j,k,mgrd,nxfr,nxto
      real*8 xg,yg,zg
c
      mgrd=ngx*ngy*ngz
      print*,ngx,ngy,ngz
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
      do 1 i=1,ngx
      do 2 j=1,ngy
      do 3 k=1,ngz
      g(i,j,k)=0.0d0
    3 continue
    2 continue
    1 continue
c
      do 4 ix=nxfr,nxto
      xg=dble(ix-1)*spac+xmin(1)
      do 5 iy=1,ngy
      yg=dble(iy-1)*spac+xmin(2)
      do 6 iz=1,ngz
      zg=dble(iz-1)*spac+xmin(3)
      do 7 i=1,nsa
      g(ix,iy,iz)=g(ix,iy,iz)+138.9354d0*
     + q(i,lplt)/sqrt((xg-x(i,1))**2+(yg-x(i,2))**2+(zg-x(i,3))**2)
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
     + 10.0d0*xmin(3),10.0d0*(xmin(3)+dble(ngz-1)*spac),
     + 10.0d0*xmin(2),10.0d0*(xmin(2)+dble(ngy-1)*spac),
     + 10.0d0*xmin(1),10.0d0*(xmin(1)+dble(ngx-1)*spac)
 1002 format(6e12.5)
      write(lfnplt,1003) (((g(ix,iy,iz),ix=1,ngx),iy=1,ngy),iz=1,ngz)
 1003 format(2e20.12)
      close(unit=lfnplt)
c
c     write coordinates in xyz format
c
      open(unit=lfnxyz,file=filxyz(1:index(filxyz,' ')-1),
     + form='formatted',status='unknown')
      rewind(lfnxyz)
      write(lfnxyz,2000) nsa
 2000 format('* Coordinates in XYZ format',/,i5)
      do 8 i=1,nsa
      write(lfnxyz,2001) t(i),
     + (10.0d0*x(i,j),j=1,3)
 2001 format(a2,1x,3f12.6)
    8 continue
      close(lfnxyz)
      endif
c
      return
      end
