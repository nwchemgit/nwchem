      subroutine esp_grdcnt(x,q)
c
c $Id$
c
      implicit none
c
#include "esp_params.fh"
#include "esp_common.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      real*8 x(3,natoms),q(natoms),g(3)
c
      integer i,ia,ix,iy,iz
      real*8 d,dmin,rad2,rcut2
c
      mxgrid=0
      if(me.eq.0) then
      rcut2=rcut*rcut
      mxgrid=1
      do 1 iz=0,mgrid(3)+1
      do 2 iy=0,mgrid(2)+1
      do 3 ix=0,mgrid(1)+1
      g(1)=xmin(1)-rcut+dble(ix)*spac
      g(2)=xmin(2)-rcut+dble(iy)*spac
      g(3)=xmin(3)-rcut+dble(iz)*spac
      dmin=rcut2
      ia=1
      do 4 i=1,natoms
      rad2=(cnm2au*radius(int(q(i))))**2
      d=(g(1)-x(1,i))*(g(1)-x(1,i))
     + +(g(2)-x(2,i))*(g(2)-x(2,i))
     + +(g(3)-x(3,i))*(g(3)-x(3,i))
c      if(d.lt.rad2) goto 3
       if(rad2-d.gt.small) goto 3
c      if(d.lt.dmin) dmin=d
      if(dmin-d.gt.small) dmin=d
    4 continue
c      if(dmin.lt.rcut2) then
      if(rcut2-dmin.gt.small) then
      mxgrid=mxgrid+1
      endif
    3 continue
    2 continue
    1 continue
      endif
c
      if(np.gt.0) then
      call ga_brdcst(mre_001,mxgrid,ma_sizeof(mt_int,1,mt_byte),0)
      endif
c
      return
      end
