      subroutine esp_rms(x,qf,g,v)
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
      real*8 x(3,natoms),qf(ndim,6)
      real*8 g(3,mxgrid),v(mxgrid)
      double precision v2(6)
c
      integer i,j,k
      real*8 eq
c
      do 1 k=1,nf
      rms(1,k)=zero
      rms(2,k)=zero
      v2(k)=0d0
      do 2 i=1,ngrid
      eq=zero
      do 3 j=1,natoms
      eq=eq+qf(j,k)/dsqrt((g(1,i)-x(1,j))*(g(1,i)-x(1,j))+
     + (g(2,i)-x(2,j))*(g(2,i)-x(2,j))+(g(3,i)-x(3,j))*(g(3,i)-x(3,j)))
    3 continue
      rms(1,k)=rms(1,k)+(v(i)-eq)*(v(i)-eq)
      rms(2,k)=rms(2,k)+(v(i)-eq)*(v(i)-eq)
      v2(k)=v2(k)+v(i)*v(i)
    2 continue
    1 continue
      do 4 k=1,nf
      rms(1,k)=cau2kj*dsqrt(rms(1,k))/dble(ngrid)
      rms(2,k)=100d0*sqrt(rms(2,k)/v2(k))
    4 continue
c
      return
      end
