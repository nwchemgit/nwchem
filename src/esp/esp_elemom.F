      subroutine esp_elemom(x,qf)
c
c $Id: esp_elemom.F,v 1.3 1997-10-07 23:41:40 d3j191 Exp $
c
      implicit none
c
#include "esp_params.fh"
#include "esp_common.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
c
      real*8 x(3,natoms),qf(ndim,4)
c
      integer i,ix,j
      real*8 qm,qp,xm(3),xp(3)
c
      do 3 j=1,nf
      qm=zero
      qp=zero
      do 1 ix=1,3
      xm(ix)=zero
      xp(ix)=zero
    1 continue
      do 2 i=1,natoms
      if(qf(i,j).lt.zero) then
      qm=qm+qf(i,j)
      xm(1)=xm(1)+qf(i,j)*x(1,i)
      xm(2)=xm(2)+qf(i,j)*x(2,i)
      xm(3)=xm(3)+qf(i,j)*x(3,i)
      else
      qp=qp+qf(i,j)
      xp(1)=xp(1)+qf(i,j)*x(1,i)
      xp(2)=xp(2)+qf(i,j)*x(2,i)
      xp(3)=xp(3)+qf(i,j)*x(3,i)
      endif
    2 continue
      dipmom(j)=half*(qp-qm)*
     + dsqrt((xp(1)/qp-xm(1)/qm)*(xp(1)/qp-xm(1)/qm)+
     +       (xp(2)/qp-xm(2)/qm)*(xp(2)/qp-xm(2)/qm)+
     +       (xp(3)/qp-xm(3)/qm)*(xp(3)/qp-xm(3)/qm))
    3 continue
c
      return
      end
