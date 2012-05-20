      SUBROUTINE tool_volme(latt,vol)

      implicit none

      real*8 x,y,z,latt,vol

      dimension latt(3,3)

      x=latt(2,2)*latt(3,3)-latt(2,3)*latt(2,3)
      y=latt(3,2)*latt(1,3)-latt(1,2)*latt(3,3)
      z=latt(1,2)*latt(2,3)-latt(2,2)*latt(1,3)

      vol=abs(latt(1,1)*x+latt(2,1)*y+latt(3,1)*z)

      return

      END
c $Id$
