      subroutine esp_getx(t,x,q,n)
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
c
      character*16 t(*)
      integer n(*)
      real*8 x(3,*),q(*)
      integer i,j
      integer natoms_in
c
      
      if(.not.geom_cart_get(igeom,natoms_in,t,x,q))
     +     call md_abort('esp: geom_cart_get failed',me)
c     get rid of zero charges
      natoms=0
      do i=1,natoms_in
         if(abs(q(i)).gt.1d-9) then
            natoms=natoms+1
            t(natoms)=t(i)
            q(natoms)=q(i)
            x(3,natoms)=x(3,i)
         endif
      enddo
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
