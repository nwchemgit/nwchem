      real*8 function ana_bond(x,nx,i,j)
c
c $Id: ana_bond.F,v 1.1 2000-03-16 20:27:40 d3j191 Exp $
c
      implicit none
c
      integer nx,i,j
      real*8 x(nx,3)
c
c     determine the bond length
c     -------------------------
c
      ana_bond=sqrt((x(i,1)-x(j,1))**2+(x(i,2)-x(j,2))**2+
     + (x(i,3)-x(j,3))**2)
c
      return
      end
