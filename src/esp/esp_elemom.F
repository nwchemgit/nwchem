      subroutine esp_elemom(x,q)
c
c $Id: esp_elemom.F,v 1.1 1997-10-06 23:00:32 d3j191 Exp $
c
      implicit none
c
#include "esp_params.fh"
#include "esp_common.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
c
      real*8 x(3,natoms),q(natoms)
c
      integer i,ix
      real*8 qm,qp,xm(3),xp(3)
c
      qm=0.0d0
      qp=0.0d0
      do 1 ix=1,3
      xm(ix)=0.0d0
      xp(ix)=0.0d0
    1 continue
      do 2 i=1,natoms
      if(q(i).lt.0.0d0) then
      qm=qm+q(i)
      xm(1)=xm(1)+q(i)*x(1,i)
      xm(2)=xm(2)+q(i)*x(2,i)
      xm(3)=xm(3)+q(i)*x(3,i)
      else
      qp=qp+q(i)
      xp(1)=xp(1)+q(i)*x(1,i)
      xp(2)=xp(2)+q(i)*x(2,i)
      xp(3)=xp(3)+q(i)*x(3,i)
      endif
    2 continue
      dipmom=0.5d0*(qp-qm)*
     + dsqrt((xp(1)/qp-xm(1)/qm)*(xp(1)/qp-xm(1)/qm)+
     +       (xp(2)/qp-xm(2)/qm)*(xp(2)/qp-xm(2)/qm)+
     +       (xp(3)/qp-xm(3)/qm)*(xp(3)/qp-xm(3)/qm))
c
      return
      end
