      subroutine esp_getx(t,x,q,n)
c
c $Id: esp_getx.F,v 1.4 2003-10-19 06:07:31 d3j191 Exp $
c
      implicit none
c
#include "esp_params.fh"
#include "esp_common.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
c
      character*16 t(natoms)
      integer n(natoms)
      real*8 x(3,natoms),q(natoms)
      integer i,j
c
      if(.not.geom_cart_get(igeom,natoms,t,x,q))
     + call md_abort('esp: geom_cart_get failed',me)
c
c     convert coordinates to nm and find extrema
c
      do 1 j=1,3
      xmin(j)=x(j,1)
      xmax(j)=x(j,1)
      do 2 i=1,natoms
      if(x(j,i).lt.xmin(j)) xmin(j)=x(j,i)
      if(x(j,i).gt.xmax(j)) xmax(j)=x(j,i)
    2 continue
    1 continue
      do 3 i=1,natoms
      n(i)=int(q(i))
    3 continue
c
c     calculate grid size
c
      do 4 i=1,3
      mgrid(i)=(xmax(i)-xmin(i)+2.0d0*rcut)/spac+3
    4 continue
      mxgrid=mgrid(1)*mgrid(2)*mgrid(3)
c
      return
      end
