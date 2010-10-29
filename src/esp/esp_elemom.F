      subroutine esp_elemom(x,qf)
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
      real*8 x(3,natoms),qf(ndim,6)
c
      real*8 xc(3),dm,ch,xq(3),rq
      integer i,ix,j,jx,k,l
c
c     dipole moment
c
      do 1 j=1,nf
      dipmom(j)=zero
      do 2 ix=1,3
      xc(ix)=zero
      ch=zero
      do 3 i=1,natoms
      ch=ch+abs(qf(i,j))
      xc(ix)=xc(ix)+x(ix,i)*abs(qf(i,j))
    3 continue
      if(ch.gt.zero) then
      xc(ix)=xc(ix)/ch
      dm=zero
      do 4 i=1,natoms
      dm=dm+qf(i,j)*(x(ix,i)-xc(ix))
    4 continue
      dipmom(j)=dipmom(j)+dm*dm
      endif
    2 continue
      dipmom(j)=dsqrt(dipmom(j))
    1 continue
c
c     quadrupole moment
c
      do 5 k=1,nf
      do 6 i=1,3
      do 7 j=1,3
      quamom(i,j,k)=zero
    7 continue
    6 continue
      do 8 l=1,natoms
      rq=zero
      do 9 ix=1,3
      xq(ix)=(x(ix,l)-xc(ix))
      rq=rq+xq(ix)*xq(ix)
    9 continue
      do 10 ix=1,3
      quamom(ix,ix,k)=quamom(ix,ix,k)+qf(l,k)*(three*xq(ix)*xq(ix)-rq)
      do 11 jx=ix+1,3
      quamom(ix,jx,k)=quamom(ix,jx,k)+qf(l,k)*(three*xq(ix)*xq(jx))
   11 continue
   10 continue
    8 continue
      do 12 ix=1,2
      do 13 jx=ix+1,3
      quamom(jx,ix,k)=quamom(ix,jx,k)
   13 continue
   12 continue
      call esp_jacobi(quamom(1,1,k))
    5 continue
c
      return
      end
