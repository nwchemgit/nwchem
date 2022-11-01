      subroutine esp_plt(x,q,g,t)
c
c $Id$
c
      implicit none
c
#include "esp_params.fh"
#include "esp_common.fh"
#include "global.fh"
#include "msgids.fh"
#include "geom.fh"
#include "util.fh"
#include "errquit.fh"
c
      real*8 x(3,natoms),q(ndim,4)
      real*8 g(*)
      character*16 t(natoms)
c
      integer i,ix,iy,iz,j,jf,jl,num
      real*8 xc,yc,zc,eq
      character*10 today,now
c
      if(mgrid(1)+2.gt.3*mxgrid) call md_abort('Memory problem',0)
c
      if(me.eq.0) then
      call util_file_name('plt',.false.,.false.,pltfil)
      open(unit=lfnplt,file=pltfil,form='formatted',status='unknown')
      write(lfnplt,1000) 3,3
 1000 format(2i7)
      write(lfnplt,1001) mgrid(3)+2,mgrid(2)+2,mgrid(1)+2
 1001 format(3i7)
      write(lfnplt,1002)
     + 10.0d0*cau2nm*(xmin(3)-rcut),
     + 10.0d0*cau2nm*(xmin(3)-rcut+dble(mgrid(3)+1)*spac),
     + 10.0d0*cau2nm*(xmin(2)-rcut),
     + 10.0d0*cau2nm*(xmin(2)-rcut+dble(mgrid(2)+1)*spac),
     + 10.0d0*cau2nm*(xmin(1)-rcut),
     + 10.0d0*cau2nm*(xmin(1)-rcut+dble(mgrid(1)+1)*spac)
 1002 format(6e12.5)
      endif
c
      if(np.gt.1) then
      num=natoms/np+1
      jf=me*num
      jl=(me+1)*num-1
      if(jl.gt.natoms) jl=natoms
      else
      jf=1
      jl=natoms
      endif
c
      do 1 iz=0,mgrid(3)+1
      zc=xmin(3)-rcut+dble(iz)*spac
      do 2 iy=0,mgrid(2)+1
      yc=xmin(2)-rcut+dble(iy)*spac
      do 3 ix=0,mgrid(1)+1
c      if(me.eq.0) then
c      call swatch(today,now)
c      write(*,1005) today,now,ix,iy,iz
c 1005 format(2a10,': plot point ',3i10)
c      endif
      xc=xmin(1)-rcut+dble(ix)*spac
      eq=zero
      if(jf.le.jl) then
      do 4 j=jf,jl
      eq=eq+q(j,nf)/dsqrt((xc-x(1,j))*(xc-x(1,j))+
     + (yc-x(2,j))*(yc-x(2,j))+(zc-x(3,j))*(zc-x(3,j)))
    4 continue
      endif
      g(ix+1)=eq
    3 continue
      if(np.gt.1) call ga_dgop(mre_003,g,mgrid(1)+2,'+')
      if(me.eq.0) then
      call swatch(today,now)
      write(lfnplt,1003) (g(ix),ix=1,mgrid(1)+2)
 1003 format(2e20.12)
c      write(*,1004) today,now,iy,iz
c 1004 format(2a10,': plot grid column ',2i10)
      endif
    2 continue
    1 continue
c
      if(me.eq.0) then
      close(lfnplt)
c
      call util_file_name('xyz',.false.,.false.,xyzfil)
      open(unit=lfnxyz,file=xyzfil,form='formatted',status='unknown')
#if 0      
      write(lfnxyz,2000) natoms
 2000 format('* Coordinates in XYZ format',/,i5)
      do 5 i=1,natoms
      write(lfnxyz,2001) t(i)(1:index(t(i),' ')-1),
     + (10.0d0*cau2nm*x(j,i),j=1,3),q(i,4)
 2001 format(a,t3,4f12.6)
    5 continue
#else
      if(.not.geom_print_xyz(igeom, lfnxyz))
     A     call errquit(' failed geom_print_xyz',
     E     0, GEOM_ERR)
#endif      
      close(lfnxyz)
      endif
c
      return
      end
