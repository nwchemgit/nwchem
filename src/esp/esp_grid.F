      subroutine esp_grid(x,q,g)
c
c $Id: esp_grid.F,v 1.1 1997-10-06 23:00:34 d3j191 Exp $
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
      real*8 x(3,natoms),q(natoms),g(3,mxgrid)
      real*8 d,dmin,rcut2,rad2
      integer i,ix,iy,iz
c
c     determine actual number of grid points
c
      if(me.eq.0) then
      rcut2=rcut*rcut
      ngrid=1
      do 1 iz=1,mgrid(3)
      do 2 iy=1,mgrid(2)
      do 3 ix=1,mgrid(1)
      g(1,ngrid)=xmin(1)-rcut+dble(ix)*spac
      g(2,ngrid)=xmin(2)-rcut+dble(iy)*spac
      g(3,ngrid)=xmin(3)-rcut+dble(iz)*spac
      dmin=rcut2
      do 4 i=1,natoms
      rad2=(cnm2au*radius(int(q(i))))**2
      d=(g(1,ngrid)-x(1,i))*(g(1,ngrid)-x(1,i))
     + +(g(2,ngrid)-x(2,i))*(g(2,ngrid)-x(2,i))
     + +(g(3,ngrid)-x(3,i))*(g(3,ngrid)-x(3,i))
      if(d.lt.rad2) goto 3
      if(d.lt.dmin) dmin=d
    4 continue
      if(dmin.lt.rcut2) ngrid=ngrid+1
    3 continue
    2 continue
    1 continue
      ngrid=ngrid-1
      endif
c
c     broadcast grid to other nodes
c
      if(np.gt.0) then
      call ga_brdcst(mre_001,ngrid,ma_sizeof(mt_int,1,mt_byte),0)
      call ga_brdcst(mre_002,g,3*ngrid*ma_sizeof(mt_dbl,1,mt_byte),0)
      endif
c
      return
      end
